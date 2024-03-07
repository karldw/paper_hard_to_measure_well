# Limit the number of threads OpenBLAS uses, since we're using snakemake to run
# in parallel at a higher level. (Only affects OpenBLAS; must be run before
# we import numpy)
# Also guard against MKL bugs, in case using MKL.
# See https://www.alexpghayes.com/post/2022-10-18-intel-mkl-data-race/
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["MKL_NUM_THREADS"] = "1"


# conda:
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import scipy.sparse
import scipy.stats
import pyprojroot
import ipopt  # installed from conda as cyipopt, except on windows


# standard lib
from pathlib import Path
from collections import namedtuple
from functools import partial
import json
import datetime
import logging
import re
import warnings
import contextlib
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

# The math details (FOC etc) are in outcomes_analysis_derivative_functions.py
# The logic of using those FOC are in this file.
from outcomes_analysis_helpers import (
    ProblemSetup,
    dwl_per_well,
    calc_well_audit_prob_uniform,
    prob_leak_with_policy,
    is_probability,
    read_constants,
    abatement_cost_per_pad,
    AuditInfo,
    DataParam,
    OutcomesOneDraw,
    FEE_PERCENTILES_TO_REPORT,
)

const = read_constants()
TAU_LEVELS = const["TAU_LEVELS"]
T_LEVELS = const["T_LEVELS"]
AUDIT_COST_TO_COMPARE_FIXED_VS_OPTIM_DWL = const[
    "AUDIT_COST_TO_COMPARE_FIXED_VS_OPTIM_DWL"
]
SOCIAL_COST_METHANE_PER_KG = const["SOCIAL_COST_METHANE_PER_KG"]
MODEL_NAMES = const["MODEL_NAMES"]
METHANE_GWP = const["METHANE_GWP"]


def identity(x):
    return x


def _helper_common_cost_param(b, Y, X):
    resid = Y - X @ np.atleast_2d(b).T
    return resid[:, 0]


def calc_common_cost_param(prob_leak, e_size_expect, price_H):
    """
    Use least squares to solve for the A and α coefficients of the cost function
    for this MCMC draw.
    (Doing this here because I'm having a hard time getting it to work in Stan)
    """
    N = len(prob_leak)
    assert N == len(e_size_expect) == len(price_H) > 1
    e_p = e_size_expect * price_H
    Y = np.expand_dims(np.log(e_p), axis=1)
    X = np.column_stack((np.ones_like(prob_leak), np.log(prob_leak)))
    # These are the max acceptable values (with some fudge factor) to make the
    # cost functions work out.
    # Here we estimate with constrained least squares.
    max_acceptable_A_log = np.log(np.min(e_p) * 0.999)
    max_acceptable_α = -1.001
    soln = scipy.optimize.least_squares(
        fun=_helper_common_cost_param,
        x0=np.array((max_acceptable_A_log, max_acceptable_α)),
        bounds=(
            (-np.inf, -np.inf),  # lower bounds
            (max_acceptable_A_log, max_acceptable_α),  # upper bounds
        ),
        args=(Y, X),
        method="dogbox",
    )
    coefs = soln.x
    A = np.exp(coefs[0])
    α = coefs[1]
    assert A <= np.min(e_p)
    assert α <= max_acceptable_α
    coefs_mat = np.full((N, 2), np.array([[A, α]]))
    return coefs_mat


def read_generated_iter(input, varnames):
    """Read the parquet files from the Stan generated quantities
    Provides a generator over MCMC draws to operate on each individually

    param: input is snakemake.input or other dictionary-like with file paths
    varnames is the variable name to read. Variables are stored one per file.
    """
    assert len(varnames) >= 1
    pq_files = [pq.ParquetFile(input[v], memory_map=True) for v in varnames]
    ncol = pq_files[0].metadata.num_columns
    for draw_id in range(ncol):
        draw_id = str(draw_id + 1)  # convert from R index to python zero-based
        draws = [
            f.read(draw_id).column(0).to_pandas(zero_copy_only=True).to_numpy()
            for f in pq_files
        ]
        yield draws


def read_generated_col(draw_id, input, varnames):
    """Read the parquet files from the Stan generated quantities
    Provides a generator over MCMC draws to operate on each individually

    param: input is snakemake.input or other dictionary-like with file paths
    varnames is the variable name to read. Variables are stored one per file.
    """
    assert len(varnames) >= 1
    draw_id = str(draw_id + 1)  # convert from R index to python zero-based
    pq_files = [pq.ParquetFile(input[v], memory_map=True) for v in varnames]
    draws = [
        f.read(draw_id).column(0).to_pandas(zero_copy_only=True).to_numpy()
        for f in pq_files
    ]
    return draws


def read_sdata(json_file):
    """Read the data that was input to Stan"""

    with open(json_file, "rt") as f:
        sdata = json.load(f)
    df = pd.DataFrame(sdata["X"], columns=sdata["X_varnames"])
    # Note that "Y", "noise", "price", and "gas_avg_mcfd" are not columns in X.
    # See model_data_prep.R
    for c in {"Y", "noise", "price", "gas_avg_mcfd"}:
        df[c] = sdata[c]
    return df


def run_ipopt(data_param, audit_info, r_guess):
    N = len(data_param.e_size_expect)
    if audit_info.audit_frac > 0:
        budget_const = [N * audit_info.audit_frac]
    else:
        budget_const = []
    num_constr = len(budget_const)  # (variable bounds don't count)
    r_bound_low = np.zeros_like(r_guess)
    r_bound_high = np.ones_like(r_guess)
    prob = ipopt.problem(
        n=len(r_guess),  # N variables to optimize
        m=num_constr,
        problem_obj=ProblemSetup(data_param, audit_info),
        lb=r_bound_low,
        ub=r_bound_high,
        cu=budget_const,  # no lower bound
    )
    # print_level 0 still prints some stuff, but at least we'll save a little
    # time printing progress reports into the void.
    prob.addOption("print_level", 0)
    prob.addOption(
        "acceptable_iter", 30
    )  # How long to stay at an okay but not great solution? (Higher than default of 15)
    with contextlib.redirect_stdout(None):
        r, info = prob.solve(r_guess)
    # Pull the shadow value out
    assert len(info["mult_g"]) == num_constr
    if num_constr == 1:
        λ = -info["mult_g"][0]
    elif num_constr == 0:
        λ = -audit_info.audit_cost
    else:
        raise ValueError("unexpected number of shadow values")
    successful_solve = info["status"] == 0
    if not successful_solve:
        if info["status"] == 1:
            warnings.warn(
                "Solving to desired tolerances failed, but reached 'acceptable' tolerances",
                RuntimeWarning,
            )
        else:
            raise RuntimeError(f"Solving failed, status message {info['status_msg']}")
    return r, λ


def calc_outcomes_once(audit_info, data_param, r_guess=None):
    """
    Calculate the outcomes of one policy, for one MCMC draw.

    Output:
    Instance of OutcomesOneDraw named tuple, and r

    """
    N = len(data_param.e_size_expect)
    if r_guess is None:
        if audit_info.audit_rule == "uniform":
            r_guess = np.array([audit_info.audit_frac])
        elif audit_info.audit_rule == "target_e" and audit_info.detect_threshold > 0:
            # 2 probs per well -- see math writeup.
            # (Code should also accept shape (2,))
            r_guess = np.full((2 * N,), audit_info.audit_frac)
        else:
            r_guess = np.full((N,), audit_info.audit_frac)

    # Use a bunch of if-else here (rather than the cleaner dict of fns) to
    # get better tracebacks.
    if audit_info.audit_rule == "none":
        r = np.zeros_like(data_param.e_size_expect)
        λ = 0.0
        expected_fee_charged = audit_info.τT * r
    elif audit_info.audit_rule == "remote":
        r = np.zeros_like(data_param.e_size_expect)
        λ = 0.0
        # Wells with e < audit_info.detect_threshold face no price
        expected_fee_charged = audit_info.τT * (
            data_param.e_size_draw >= audit_info.detect_threshold
        )
    elif audit_info.audit_rule == "uniform" and audit_info.audit_frac > 0:
        r, λ = calc_well_audit_prob_uniform(data_param, audit_info, r_guess)
        expected_fee_charged = audit_info.τT * r
    elif audit_info.audit_rule == "target_e" and audit_info.detect_threshold > 0:
        r_big_r_small, λ = run_ipopt(data_param, audit_info, r_guess)
        # The well knows their size, so the effective fee per kg is now r_big if
        # the e is big and r_small if the e is small. This is different than the
        # probability-mix we had in the regulator's problem.
        # See the math notes for more.
        # r_big_r_small has shape (2 * N,), with r for big-e wells and r for
        # small-e wells, but appended together into one long vector because
        # that's how ipopt needs things to be. Use np.where to pick.
        r = np.where(
            data_param.e_size_draw > audit_info.detect_threshold,
            r_big_r_small[:N],
            r_big_r_small[N:],
        )
        expected_fee_charged = audit_info.τT * r
        # In small tests, it's faster to use the uniform than to start from the last point
    else:
        r, λ = run_ipopt(data_param, audit_info, r_guess)
        expected_fee_charged = audit_info.τT * r

    # Need to have the right shape here, so special-case target_e_high
    if audit_info.audit_rule == "target_e" and audit_info.detect_threshold > 0:
        r_guess_next = r_big_r_small
    else:
        r_guess_next = r

    # Score the r chosen above, now using e_size_draw instead of e_size_expect
    # Note that dwl here will include the audit cost when audit_cost > 0, but
    # not when audit_frac > 0, so we add it in.
    # (Easier to do it here, before the calc_relative_outcomes function)
    dwl = dwl_per_well(expected_fee_charged, data_param)
    if audit_info.audit_frac > 0 and audit_info.audit_cost == 0:
        additional_audit_cost_mean = (
            audit_info.audit_frac * AUDIT_COST_TO_COMPARE_FIXED_VS_OPTIM_DWL
        )
    else:
        additional_audit_cost_mean = 0.0
    additional_audit_cost_tot = additional_audit_cost_mean * N
    new_prob_leak = prob_leak_with_policy(
        expected_fee_charged,
        data_param.cost_coef,
        data_param.e_size_draw,
        data_param.price_H,
    )
    bau_prob_leak = prob_leak_with_policy(
        np.zeros_like(expected_fee_charged),
        data_param.cost_coef,
        data_param.e_size_draw,
        data_param.price_H,
    )
    emis = data_param.e_size_draw * data_param.time_H * new_prob_leak
    # Note: this is a little redundant to calculate BAU values here and
    # separately in calc_all_outcomes_per_draw, but I can live with it. I want
    # it available here, before aggregating across wells, but I don't want to
    # change the rest of the code that works with aggregated data.

    emis_bau = data_param.e_size_draw * data_param.time_H * bau_prob_leak
    # Note: "expected_fee_charged" is the well's expectation of the fee they'll
    # face per unit of leak rate (i.e. tau * T * r_i = total fee / e_i). The
    # units are dollars-hours per kg. That works for thinking about incentives,
    # but we also want to make comparisons with the social cost, so provide a
    # number that's just dollars per kilogram (using the true kg of emissions,
    # according to the model).
    # Total expected fee: tau * T * r_i * e_i * (1 - q_i)
    # Total expected emissions: e_i * (1 - q_i) * H
    # Expected fee per unit of emissions: tau * T * r_i / H
    expected_fee_per_kg = expected_fee_charged / data_param.time_H
    fee_quantiles = np.quantile(
        expected_fee_per_kg, tuple(FEE_PERCENTILES_TO_REPORT.values())
    )
    # Note that tot_cost_per_kg is the expected cost to the operator of having
    # one more kg emissions: the expected fee plus the commodity value lost.
    # It does *not* include other costs, like abatement or audits.
    gas_price_per_kg_ch4 = data_param.price_H / data_param.time_H

    tot_cost_per_kg = expected_fee_per_kg + gas_price_per_kg_ch4

    # Now, we also want to compare to the gas price, so convert back from kg CH4
    # to mcf natural gas (see notes in model_data_prep.r and
    # match_jpl_measurements.R)
    methane_kg_per_mcf = 18.8916  # Note: this is for pure CH4, not nat gas
    approx_ch4_fraction = 0.95
    gas_price_per_mcf_gas = (
        gas_price_per_kg_ch4 * methane_kg_per_mcf * approx_ch4_fraction
    )

    # direct_private_cost_per_pad answers the question "how much do the well
    # pad's costs go up per period H?" This includes fees and abatement costs.
    # It does notinclude changes in revenue (from additional gas captured).
    direct_private_cost_per_pad = (
        abatement_cost_per_pad(new_prob_leak, data_param)
        - abatement_cost_per_pad(bau_prob_leak, data_param)
        + expected_fee_per_kg * emis
    )
    # net private cost = abatement cost + expected fee - additional revenue
    gas_prod_mcf_per_H = data_param.gas_avg_mcfd / 24.0 * data_param.time_H
    # We want the fractional change (fraction of price, rather than $) because
    # we'll use this number in an elasticity calculation later.
    net_private_cost_per_mcf_pct_price = (
        100
        * (direct_private_cost_per_pad - gas_price_per_kg_ch4 * (emis_bau - emis))
        / gas_prod_mcf_per_H
        / gas_price_per_mcf_gas
    )
    net_private_cost_per_mcf_pct_price_weighted = np.average(
        net_private_cost_per_mcf_pct_price, weights=data_param.gas_avg_mcfd
    )

    outcomes = OutcomesOneDraw(
        dwl_mean=np.mean(dwl) + additional_audit_cost_mean,
        dwl_tot=np.sum(dwl) + additional_audit_cost_tot,
        emis_mean=np.mean(emis),
        emis_tot=np.sum(emis),
        tot_cost_per_kg_mean=np.mean(tot_cost_per_kg),
        fee_per_kg_mean=np.mean(expected_fee_per_kg),
        fee_per_kg_p01=fee_quantiles[0],
        fee_per_kg_p05=fee_quantiles[1],
        fee_per_kg_p10=fee_quantiles[2],
        fee_per_kg_p20=fee_quantiles[3],
        fee_per_kg_p25=fee_quantiles[4],
        fee_per_kg_p30=fee_quantiles[5],
        fee_per_kg_p40=fee_quantiles[6],
        fee_per_kg_med=fee_quantiles[7],
        fee_per_kg_p60=fee_quantiles[8],
        fee_per_kg_p70=fee_quantiles[9],
        fee_per_kg_p75=fee_quantiles[9],
        fee_per_kg_p80=fee_quantiles[11],
        fee_per_kg_p90=fee_quantiles[12],
        fee_per_kg_p95=fee_quantiles[13],
        fee_per_kg_p99=fee_quantiles[14],
        net_private_cost_per_mcf_pct_price=net_private_cost_per_mcf_pct_price_weighted,
        shadow_price=λ,
        audit_rule=audit_info.audit_rule,
        audit_frac=audit_info.audit_frac,
        τT=audit_info.τT,
        detect_threshold=audit_info.detect_threshold,
        audit_cost=audit_info.audit_cost,
    )
    # returning r here is only useful to use as r_guess for the next iter
    # (and in testing)
    return outcomes, r_guess_next


def check_outcomes(df):
    """Run some tests on the outcomes. Return the original df"""

    # 1. Uniform should have uniform fee
    df_uni = df.query("audit_rule == 'uniform'")
    assert (df_uni.fee_per_kg_mean == df_uni.fee_per_kg_med).all()

    return df


def get_input_files(snakemake):
    """
    All model-generate produce leak_size_expect.parquet, stan_data.json,
    and leak_size_draw.parquet.
    The code also depends on either prob_leak.parquet or cost_param_A.parquet
    and cost_param_alpha.parquet in the same folder.
    Read those in based on the model name.
    """
    input_dir = Path(snakemake.input["stan_data_json"]).resolve().parent
    assert input_dir.is_dir()
    model_name = snakemake.wildcards["model_name"]
    input_files = {
        "stan_data_json": input_dir / "stan_data.json",
        "leak_size_expect": input_dir / "leak_size_expect.parquet",
        "leak_size_draw": input_dir / "leak_size_draw.parquet",
        "prob_size_above_threshold": input_dir / "prob_size_above_threshold.parquet",
    }
    if model_name in MODEL_NAMES["cost_coef_models"]:
        input_files["cost_param_A"] = input_dir / "cost_param_A.parquet"
        input_files["cost_param_alpha"] = input_dir / "cost_param_alpha.parquet"
    else:
        input_files["prob_leak"] = input_dir / "prob_leak.parquet"
    return input_files


def calc_all_outcomes_per_draw(
    draw_id, input_files, audit_info, price_H, time_H, gas_avg_mcfd, r_guess=None
):
    # Read and unpack common values from parquet files.
    data_to_read = (
        "leak_size_expect",
        "leak_size_draw",
        "prob_size_above_threshold",
    )
    data_list = read_generated_col(draw_id, input_files, data_to_read)
    e_size_expect = data_list[0]
    e_size_draw = data_list[1]
    prob_size_above_threshold = data_list[2]

    if "cost_param_alpha" in input_files.keys():
        cost_coef = np.column_stack(
            read_generated_col(
                draw_id, input_files, ("cost_param_A", "cost_param_alpha")
            )
        )
    else:
        cost_coef = calc_common_cost_param(
            prob_leak=read_generated_col(draw_id, input_files, ["prob_leak"]),
            e_size_expect=e_size_expect,
            price_H=price_H,
        )
    data_param = DataParam(
        price_H=price_H,
        e_size_expect=e_size_expect,
        e_size_draw=e_size_draw,
        cost_coef=cost_coef,
        time_H=time_H,
        prob_is_large=prob_size_above_threshold,
        gas_avg_mcfd=gas_avg_mcfd,
    )

    try:
        outcomes_policy, _ = calc_outcomes_once(audit_info, data_param, r_guess)
    except RuntimeError:
        raise RuntimeError(f"Failed to converge for draw {draw_id}, audit {audit_info}")
    # Also do BAU and optimal. These are faster, since they're analytical
    audit_bau = AuditInfo(
        audit_rule="none", audit_frac=0.0, τT=0.0, detect_threshold=0.0, audit_cost=0.0
    )
    audit_optimal = AuditInfo(
        audit_rule="remote",
        audit_frac=0.0,
        τT=SOCIAL_COST_METHANE_PER_KG * time_H,
        detect_threshold=0.0,
        audit_cost=0.0,
    )

    outcome_bau, _ = calc_outcomes_once(audit_bau, data_param)
    outcome_optimal, _ = calc_outcomes_once(audit_optimal, data_param)
    return outcomes_policy, outcome_bau, outcome_optimal


def calc_all_outcomes_all_draws(snakemake):
    """
    Calculate all outcomes for all MCMC draws.
    (this takes a while, depending on the number of cases considered)
    """
    input_files = get_input_files(snakemake)
    sdata = read_sdata(input_files["stan_data_json"])
    price = sdata["price"]  # gas_price_per_kg_ch4
    gas_avg_mcfd = sdata["gas_avg_mcfd"]
    time_H = parse_period_wildcard(snakemake.wildcards["time_period"])
    price_H_dollar_hr_per_kg = price * time_H
    audit_info = parse_audit_info(snakemake.wildcards)
    logging.info(audit_info)
    num_draws = 4000
    # Here I designed the loop to allow you to provide last round's r as a guess
    # for this round's r, thinking that would allow ipopt to converge faster.
    # It's actually slower, so we'll just provide r_guess = None for every iter.
    r_guess = None
    outcome_policy_list = []
    outcome_bau_list = []
    outcome_optimal_list = []

    # Windows will throw an error if >61
    # if snakemake.threads is 1, this will still launch a sub-process.
    mp_threads = max(min(snakemake.threads - 1, 61), 1)
    mp_context = multiprocessing.get_context(method="spawn")
    mp_function = partial(
        calc_all_outcomes_per_draw,
        input_files=input_files,
        audit_info=audit_info,
        price_H=price_H_dollar_hr_per_kg,
        time_H=time_H,
        gas_avg_mcfd=gas_avg_mcfd,
        r_guess=r_guess,
    )
    # For the slow ones, process in chunks. For the fast ones, just do it all in one.
    if audit_info.audit_rule in {"target_e", "target_x"}:
        mp_chunksize = 50
    else:
        mp_chunksize = num_draws
    with ProcessPoolExecutor(max_workers=mp_threads, mp_context=mp_context) as executor:
        parallel_results = executor.map(mp_function, range(num_draws), chunksize=mp_chunksize)

    for res in parallel_results:
        outcomes_policy, outcome_bau, outcome_optimal = res
        outcome_policy_list.append(outcomes_policy)
        outcome_bau_list.append(outcome_bau)
        outcome_optimal_list.append(outcome_optimal)

    df_policy = pd.DataFrame.from_records(
        outcome_policy_list, columns=OutcomesOneDraw._fields
    )
    df_bau = pd.DataFrame.from_records(
        outcome_bau_list, columns=OutcomesOneDraw._fields
    )
    df_optimal = pd.DataFrame.from_records(
        outcome_optimal_list, columns=OutcomesOneDraw._fields
    )
    df_policy["draw_id"] = range(num_draws)
    df_bau["draw_id"] = range(num_draws)
    df_optimal["draw_id"] = range(num_draws)

    df_all = calc_relative_outcomes(
        df_policy, df_bau, df_optimal, time_H, gas_avg_mcfd
    ).pipe(check_outcomes)
    return df_all


def calc_relative_outcomes(df_policy, df_bau, df_optimal, time_H, gas_avg_mcfd):
    """
    Calculate outcomes relative to the first best outcome for
    variables dwl_mean, dwl_tot, emis_mean, and emis_tot.
    1. Relative outcomes, e.g. emis_mean_rel_pct, are on a scale from 0 (BAU) to
    100 (first best). Note that first-best emissions won't be zero kg.
    2. Difference in outcomes, e.g. emis_reduce_mean, from BAU.
    """
    # Don't construct these from strings so I can grep for them later.
    var_to_normalize = {
        "dwl_mean": "dwl_mean_rel_pct",
        "dwl_tot": "dwl_tot_rel_pct",
        "emis_mean": "emis_mean_rel_pct",
        "emis_tot": "emis_tot_rel_pct",
    }
    # Isolate the BAU outcomes (one row per draw) and select the variables we
    # want to calculate from, plus the draw_id
    cols_to_keep = list(var_to_normalize.keys())
    cols_to_keep.append("draw_id")

    expected_shapes = (df_optimal.shape == df_bau.shape == df_policy.shape) and (
        df_bau.shape[0] > 0
    )
    if not expected_shapes:
        raise ValueError(
            "Bad shapes for first best / BAU outcomes. See output in log   ."
        )

    df = df_policy.merge(
        df_bau,
        how="inner",
        on="draw_id",
        validate="1:1",
        copy=True,
        suffixes=(None, "_bau"),
    ).merge(
        df_optimal,
        how="inner",
        on="draw_id",
        validate="1:1",
        copy=True,
        suffixes=(None, "_best"),
    )
    vars_to_drop = []
    # Calculate the relative change:
    for current, new_var in var_to_normalize.items():
        bau = current + "_bau"
        best = current + "_best"
        # For both DWL and emiss: 0 <= best <= current <= BAU
        # (Because audit costs are included, it's possible current > BAU for
        # very ineffective policies. For numerical issues, we can end up with
        # very small negative nuumbers as well.)
        df[new_var] = 100 * (df[bau] - df[current]) / (df[bau] - df[best])
        vars_to_drop.append(bau)
        vars_to_drop.append(best)
        assert df[new_var].between(-1, 110, inclusive="both").all()
    var_to_diff = {
        "emis_mean": "emis_reduce_mean",
        "emis_tot": "emis_reduce_tot",
    }
    # Check that the bau_outcomes df has the keys we need. (<= is subset)
    assert set(var_to_diff.keys()) <= set(cols_to_keep)
    # Now also calculate the change in levels:
    for current, new_var in var_to_diff.items():
        bau = current + "_bau"
        df[new_var] = df[bau] - df[current]
    # Clean up the *_bau and *_best variables
    df.drop(columns=vars_to_drop, inplace=True)
    return df


def conf_low(x):
    """Return the 2.5% quantile of x"""
    return x.quantile(0.025)


def conf_high(x):
    """Return the 97.5% quantile of x"""
    return x.quantile(0.975)


def summarize_outcomes(df):
    """
    For every outcome, calculate the mean and 95% CI.
    """
    assert df.notna().all().all()
    agg_fn = ["mean", conf_low, conf_high]
    summ = (
        df.drop(columns="draw_id")
        .groupby(list(AuditInfo._fields))
        .aggregate(agg_fn)
        .reset_index()
    )
    # Here pandas has made the columns into a MultiIndex...
    # summ.rename fails here because it won't handle the column name pairs
    new_columns = []
    for c in summ.columns:
        assert isinstance(c, tuple)
        assert len(c) == 2
        if c[1] == "":
            # keep original
            new_columns.append(c[0])
        else:
            new_columns.append("_".join(c))
    summ.columns = new_columns
    summ.rename(columns={"τT": "tau_T"}, inplace=True)
    assert "shadow_price_mean" in summ.columns
    return summ


def set_memory_limit(mem_mb):
    """
    Limit the available memory to `mem_mb`. Only works on unix systems.
    "available memory" includes memory_map files.
    """
    import resource

    mem_bytes = mem_mb * 1024**2
    _, lim_hard = resource.getrlimit(resource.RLIMIT_AS)
    new_lim = (mem_bytes, lim_hard)
    resource.setrlimit(resource.RLIMIT_AS, new_lim)


def extract_regex_match(str_to_search, regex):
    """Convenience wrapper to pull out one regex group"""
    match = re.search(regex, str_to_search)
    if not match:
        raise ValueError(f"Could not match {regex} in string {str_to_search}")
    return match.group(1)


def parse_period_wildcard(wildcard_str):
    """Parse the 'time_period' Snakemake wildcard.

    If the wildcard is empty (""), return 8760.0.
    """
    if wildcard_str == "":
        return 8760.0
    time_period_hr = float(extract_regex_match(wildcard_str, r"-period_(\d+)_hours"))
    return time_period_hr


def parse_τT_wildcard(wildcard_str):
    # numbers (strict: no commas or scientific notation, decimals require 0)
    τT_number_matches = re.search(r"^\-?(\d+\.)?\d+$", wildcard_str)
    if τT_number_matches:
        return float(wildcard_str)

    τT_regex = f"({'|'.join(TAU_LEVELS.keys())})-({'|'.join(T_LEVELS.keys())})"
    τT_regex_match = re.search(τT_regex, wildcard_str)
    if not τT_regex_match:
        raise ValueError(f"Failed to parse τT from {wildcard_str}")

    τ_str = τT_regex_match.group(1)
    T_str = τT_regex_match.group(2)
    τT = TAU_LEVELS[τ_str] * T_LEVELS[T_str]
    return τT


def parse_audit_info(wildcards):
    """
    Here we parse some values out of a filename like this:
    "{model_name}{prior_only}{bootstrap}{time_period}" /
    "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.parquet"
    into an AuditInfo tuple
    """
    # Here we translate what it means for detect to be "high" or "low" into numbers.
    const = read_constants()
    detect_opts = {
        "low": const["POLICY_DETECT_THRESHOLD_LOW"],
        "high": const["POLICY_DETECT_THRESHOLD_HIGH"],
    }
    # audit_frac_opts = {"low": 0.0, "med": 0.01, "high": 0.10}
    # audit_cost_opts = {"low": 0.0, "med": 100.0, "high": 600.0}
    # Parse (none|uniform|remote_low|remote_high|target_x|target_e_low|target_e_high)
    # into an audit rule (e.g. "target_e") and detection threshold (e.g. "low")
    match_audit = re.search(
        "^(none|uniform|remote|target_x|target_e)_?(low|high)?$",
        wildcards["audit_rule"],
    )
    if not match_audit:
        raise ValueError("Failed to match audit pattern")
    audit_rule = match_audit.group(1)
    match_detect_threshold = match_audit.group(2)
    if not match_detect_threshold:
        # Note: detect_threshold isn't relevant to some rules, but still need
        # to fill in a value
        match_detect_threshold = "low"
    detect_threshold = detect_opts[match_detect_threshold]

    # parse "audit_amount", which has a wildcard constraint of
    # "(0pct|1pct|10pct|optimal-100usd|optimal-600usd)"
    if wildcards["audit_amount"].endswith("pct"):
        audit_frac = extract_regex_match(wildcards["audit_amount"], r"^(0?\.?\d+)pct$")
        audit_frac = float(audit_frac) / 100.0
        audit_cost = 0.0
    elif wildcards["audit_amount"].startswith("optimal-"):
        audit_frac = 0.0
        audit_cost = float(
            extract_regex_match(wildcards["audit_amount"], r"^optimal-(\d+)usd$")
        )
    else:
        raise ValueError(f"Failed to parse {wildcards['audit_amount']}")
    τT = parse_τT_wildcard(wildcards["audit_tauT"])
    audit_info = AuditInfo(
        audit_rule=audit_rule,
        audit_frac=audit_frac,
        τT=τT,
        detect_threshold=detect_threshold,
        audit_cost=audit_cost,
    )

    # Run a few checks to make sure we're not requsting nonsense audit policy
    if (
        audit_frac < 0
        or audit_frac > 1
        or audit_cost < 0
        or detect_threshold < 0
        or τT < 0
        or audit_rule not in {"none", "uniform", "target_x", "target_e", "remote"}
    ):
        logging.error(wildcards)
        logging.error(audit_info)
        raise ValueError("Impossible audit values")
    if (audit_frac == 0 and audit_cost == 0) or (audit_frac > 0 and audit_cost > 0):
        if audit_rule in {"uniform", "target_x", "target_e"}:
            logging.error(wildcards)
            raise ValueError(
                f"Bad combo of audit_frac ({audit_frac}) and audit_cost ({audit_cost}) for rule {audit_rule}"
            )
    if (
        audit_rule in ("none", "uniform", "target_x")
        and detect_threshold != detect_opts["low"]
    ):
        logging.error(wildcards)
        raise ValueError(
            f"Audit rules that don't depend on detect_threshold should set to {detect_opts['low']}"
        )
    if audit_rule in ("none", "remote") and (audit_cost != 0 or audit_frac != 0):
        logging.error(wildcards)
        raise ValueError(
            "Remote or no audits should set audit_cost = 0 and audit_frac = 0"
        )
    if audit_cost > 0 and detect_threshold > 0:
        logging.error(audit_info)
        raise NotImplementedError(
            "audit_cost > 0 and detect_threshold > 0 isn't implemented."
        )
    return audit_info


class Timer(object):
    """
    Context manager to time an expression (once; not for benchmarking)
    https://stackoverflow.com/a/5849861
    https://stackoverflow.com/a/3427051
    """

    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = datetime.datetime.now().replace(microsecond=0)

    def __exit__(self, type, value, traceback):
        if self.name:
            logging.info("[%s]" % self.name)
        elapsed = datetime.datetime.now().replace(microsecond=0) - self.tstart
        logging.info(f"Elapsed: {elapsed}")


def make_placeholder_snakemake():
    logging.warn("Using placeholder snakemake")
    wildcards = {
        "model_name": "08_twopart_lognormal_heterog_alpha",
        # "model_name": "01_twopart_lognormal",
        "bootstrap": "-bootstrap",
        "time_period": "-period_8760_hours",
        # "time_period": "",
        "prior_only": "",
    }
    policy_outcomes = (
        Path(pyprojroot.here("data/generated"))
        / f"policy_outcomes/{wildcards['model_name']}{wildcards['prior_only']}{wildcards['bootstrap']}{wildcards['time_period']}"
    )
    stan_fits = (
        Path(pyprojroot.here("scratch"))
        / f"stan_fits/{wildcards['model_name']}{wildcards['prior_only']}{wildcards['bootstrap']}{wildcards['time_period']}"
    )
    SnakemakePlaceholder = namedtuple(
        "SnakemakePlaceholder", ["input", "output", "threads", "resources", "wildcards"]
    )
    snakemake = SnakemakePlaceholder(
        input={
            "leak_size_draw": stan_fits / "leak_size_draw.parquet",
            "leak_size_expect": stan_fits / "leak_size_expect.parquet",
            "stan_data_json": stan_fits / "stan_data.json",
        },
        output={"results_summary": policy_outcomes / "audit_outcome_summary.parquet"},
        threads=4,
        resources={"mem_mb": 7000},
        wildcards=wildcards,
    )
    return snakemake


if __name__ == "__main__":
    # Important: do not assume 'snakemake' will be available as a global.
    # Instead, pass the snakemake object around between functions.
    if not "snakemake" in globals():
        snakemake = make_placeholder_snakemake()
    if snakemake.wildcards["model_name"] not in MODEL_NAMES["all"]:
        raise ValueError(f"Unknown model {snakemake.wildcards['model_name']}")
    logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, encoding="utf-8")
    try:
        set_memory_limit(snakemake.resources["mem_mb"])
    except:
        logging.warning("Note: failed to set memory limit.")

    with Timer(snakemake.output["results_summary"]):
        results = calc_all_outcomes_all_draws(snakemake).pipe(summarize_outcomes)
        results.to_parquet(snakemake.output["results_summary"])
