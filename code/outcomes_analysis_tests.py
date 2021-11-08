import numpy as np
import scipy.stats
import ipopt  # installed from conda as cyipopt, except on windows

from outcomes_analysis_helpers import (
    is_probability,
    DataParam,
    AuditInfo,
    ProblemSetup,
    read_constants,
)
from outcomes_analysis import (
    parse_period_wildcard,
    calc_common_cost_param,
    calc_outcomes_once,
    parse_audit_info,
)


def test_parse_period_wildcard():
    """Fool me once..."""
    assert parse_period_wildcard("-period_24_hours") == 24
    assert parse_period_wildcard("-period_2_hours") == 2
    assert parse_period_wildcard("") == 8760


def test_derivatives():
    """Use Ipopt's built-in derivative checker"""
    N = 50
    audit_frac = 0.01
    M = audit_frac * N
    # Scaled down version so the derivatives aren't as large (only tests for abs difference, not relative difference)
    τT = 4.0
    H = 10.0
    rng = np.random.default_rng(0)

    e_size = rng.lognormal(3, 2, N)
    recip_α = rng.uniform(-1, 0, N)
    recip_A = rng.uniform(0, 3, N)
    price = rng.uniform(0.5, 3, N)
    r_initial_guess = np.full(N, 0.01)
    r_bound_low = np.zeros_like(e_size)
    r_bound_high = np.ones_like(e_size)
    constr_bound_low = np.array([-np.inf])
    constr_bound_high = np.array([M])
    data_param = DataParam(
        e_size_expect=e_size,
        e_size_draw=e_size,
        time_H=H,
        price_H=price * H,
        cost_coef=np.hstack((1.0 / recip_A, 1.0 / recip_α)),
        # Note: this might not be right, but it doesn't matter:
        prob_is_large=1 - scipy.stats.lognorm.cdf(e_size, 3, 2),
        gas_avg_mcfd=e_size,  # also doesn't matter
    )
    for audit_cost in (0.0, 100.0):
        if audit_cost > 0:
            num_constr = 0
        else:
            num_constr = 1
        for audit_frac in (0.0, 0.01):
            if (audit_cost > 0) + (audit_frac > 0) != 1:
                # Exactly one of these should be positive
                continue
            for audit_rule in ("uniform", "target_x", "target_e"):
                for detect_threshold in (0.0, 100.0):
                    if audit_rule == "uniform" and audit_frac > 0:
                        # Use the analytical solution.
                        continue
                    if (
                        audit_rule == "target_e"
                        and detect_threshold > 0
                        and audit_cost > 0
                    ):
                        # NOTE this case isn't working yet.
                        continue

                    if detect_threshold > 0 and audit_rule != "target_e":
                        continue
                    if audit_rule == "uniform":
                        num_vars = 1
                    elif detect_threshold > 0 and audit_rule == "target_e":
                        num_vars = N * 2
                    else:
                        num_vars = N
                    r_initial_guess = np.full(num_vars, 0.01)
                    r_bound_low = np.zeros_like(r_initial_guess)
                    r_bound_high = np.ones_like(r_initial_guess)
                    audit_info = AuditInfo(
                        audit_rule=audit_rule,
                        audit_frac=audit_frac,
                        τT=τT,
                        detect_threshold=detect_threshold,
                        audit_cost=audit_cost,
                    )
                    try_to_print(audit_info)
                    prob = ipopt.problem(
                        n=num_vars,
                        m=num_constr,
                        problem_obj=ProblemSetup(data_param, audit_info),
                        lb=r_bound_low,
                        ub=r_bound_high,
                        cl=constr_bound_low,
                        cu=constr_bound_high,
                    )
                    prob.addOption("derivative_test", "second-order")
                    # Note: chose derivative_test_perturbation to be large enough that
                    # none of the numerical derivatives are zero when they shouldn't be
                    prob.addOption("derivative_test_perturbation", 1e-5)
                    prob.addOption("derivative_test_tol", 5e-3)
                    prob.addOption("print_level", 4)
                    # Note: I would like to capture and parse the result here, but
                    # I'm having trouble capturing it.
                    # So just read the logs, I guess.
                    _, info = prob.solve(r_initial_guess)
                    if info["status"] != 0:
                        try_to_print(audit_info)
                        raise RuntimeError(
                            f"Solving failed, status message {info['status_msg']}"
                        )


def try_to_print(x):
    """Wrapper for print to deal with unicode issues on Windows"""
    try:
        print(x)
    except:
        pass


def test_audit_assignment():
    e_size1 = 10.0 + np.array(range(10))
    e_size2 = 10.0 + (np.array(range(10)) / 10.0)
    H = 8760.0
    price_H = np.full_like(e_size1, 0.5) * H
    prob_leak = np.full_like(e_size1, 0.015) - (np.array(range(10)) / 2000)
    cost_coef = calc_common_cost_param(
        prob_leak=prob_leak, e_size_expect=e_size1, price_H=price_H
    )
    for audit_info in make_audit_possibilities_for_testing(time_H=H):
        rule = audit_info.audit_rule
        if rule in {"uniform", "none", "remote"}:
            continue
        data_param1 = DataParam(
            e_size_expect=e_size1,
            e_size_draw=e_size1,
            time_H=H,
            price_H=price_H,
            cost_coef=cost_coef,
            prob_is_large=e_size2 - 10,
            gas_avg_mcfd=e_size1,  # doesn't matter
        )
        data_param2 = DataParam(
            e_size_expect=e_size2,
            e_size_draw=e_size2,
            time_H=H,
            price_H=price_H,
            cost_coef=cost_coef,
            prob_is_large=e_size2 - 10,
            gas_avg_mcfd=e_size2,  # doesn't matter
        )
        # Tests:
        out1, r1 = calc_outcomes_once(audit_info, data_param1)
        out2, r2 = calc_outcomes_once(audit_info, data_param2)
        assert len(r1) == len(r2) == len(e_size1)
        assert is_probability(r1)
        assert is_probability(r2)
        # ptp == max - min
        all_ones = np.ones_like(r1)
        all_zero = np.zeros_like(r1)
        range_greater = (
            (np.ptp(r1) >= np.ptp(r2))
            or np.allclose(r1, all_ones)
            or np.allclose(r1, all_zero)
            or (rule == "target_e" and audit_info.detect_threshold > 0)
        )
        if not range_greater:
            try_to_print(audit_info)
            print(r1, r2)
            raise ValueError("range_greater test failed")
        if rule == "target_e":
            assert is_sorted(r1) or is_sorted(-r1)
            assert is_sorted(r2) or is_sorted(-r2)
        else:
            assert is_sorted(r1)
            assert is_sorted(r2)


def is_sorted(x):
    return np.all(x[:-1] <= x[1:])


def compare_audit_info(x1, x2):
    any_mismatch = False
    for v1, v2, field in zip(x1, x2, AuditInfo._fields):
        if v1 != v2:
            print(f"Mismatch for {field}: {v1} != {v2}")
            any_mismatch = True
    if any_mismatch:
        raise AssertionError("Values do not match. See above.")


def test_parse_audit_info():
    wildcard = {
        "audit_rule": "target_x",
        "audit_amount": "1pct",
        "audit_tauT": "med-1week",
    }
    expected = AuditInfo(
        audit_rule="target_x",
        audit_frac=0.01,
        τT=2 * 168,
        detect_threshold=0,
        audit_cost=0,
    )
    compare_audit_info(parse_audit_info(wildcard), expected)
    wildcard = {
        "audit_rule": "target_e_high",
        "audit_amount": "10pct",
        "audit_tauT": "low-3month",
    }
    expected = AuditInfo(
        audit_rule="target_e",
        audit_frac=0.1,
        τT=(5 / 1000 * METHANE_GWP) * 2190.0,
        detect_threshold=100,
        audit_cost=0,
    )
    compare_audit_info(parse_audit_info(wildcard), expected)

    wildcard = {
        "audit_rule": "remote_low",
        "audit_amount": "0pct",
        "audit_tauT": "0.38",
    }
    expected = AuditInfo(
        audit_rule="remote", audit_frac=0, τT=0.38, detect_threshold=0, audit_cost=0
    )
    compare_audit_info(parse_audit_info(wildcard), expected)

    wildcard = {
        "audit_rule": "none",
        "audit_amount": "0pct",
        "audit_tauT": "0",
    }
    expected = AuditInfo(
        audit_rule="none", audit_frac=0, τT=0, detect_threshold=0, audit_cost=0
    )
    compare_audit_info(parse_audit_info(wildcard), expected)


def benchmark_reading(input_files):
    """
    Note: reading the 4000 columns and calculating cost_coef for each takes ~9s.
    That's not a ridicuous amount of overhead to run separately for each scenario
    """
    sdata = read_sdata(input_files["stan_data_json"])
    price = sdata["price"]
    start_time = time.perf_counter()
    for e_size_expect, prob_leak in read_generated(input_files):
        cost_coef = calc_common_cost_param(
            prob_leak=prob_leak, e_size_expect=e_size_expect, price=price
        )
    elapsed = time.perf_counter() - start_time
    print(elapsed)


def make_audit_possibilities_for_testing(time_H):
    """
    Define a list of audit possibilities
    (A list of AuditInfo named tuples)
    """
    # Note: all of these could easily be extended.
    # If changing tau_opts, also change TAU_LEVELS in policy_*.R
    tau_opts = (
        2.0 * SOCIAL_COST_METHANE_PER_KG,  # high
        1.0 * SOCIAL_COST_METHANE_PER_KG,  # med
        5.0 / 1000.0 * 29.8,  # low = $5/ton CO2e
    )
    # Note that this is a small set of T values for development. We probably eventually want more.
    T_opts = (
        # 24.0,  # 1 day
        168.0,  # 1 week
        # 730.0,  # 1 month
        2190.0,  # 3 months
        # 8760.0,  # 1 year
    )
    detect_opts = (0.0, 100.0)  # low and high detection thresholds
    # Note that this is a slim set of audit_frac options for speed during
    # development. We may want more values, particularly for the shadow price plot
    audit_frac_opts = (0.0, 0.01, 0.10)
    audit_cost_opts = (0.0, 100.0, 600.0)
    # Note:
    # we use audit_frac_opts when we want to consider allocating a given audit budget.
    # we use audit_cost_opts when we want to consider the optimal number of audits, given a cost.
    # At most one of these should be non-zero (remote and none have both zero)

    # Create a list of possibilities
    # - Only 1 for "none"
    # - only detect_threshold and tau varies for "remote"
    # etc

    # audit_frac value doesn't matter for "none" and "remote", but it's helpful
    # for it to be in the set of values we consider.
    # τT doesn't matter for "none", but it's again it should be in the set.
    # possibilities is a set to avoid worrying about adding duplicates.
    possibilities = {
        # BAU, no policy:
        AuditInfo(
            audit_rule="none",
            audit_frac=0.0,
            τT=tau_opts[0] * T_opts[0],
            detect_threshold=detect_opts[0],
            audit_cost=0.0,
        )
    }

    # Add a range of remote fees. This accomplishes two goals:
    # - Make sure we have the first-best Pigouvian tax for comparison
    # - Have a range of intermediate taxes to trace out the abatement curve,
    #   re-using the existing machinery.
    abatement_curve_τT = np.linspace(
        start=0, stop=SOCIAL_COST_METHANE_PER_KG * time_H * 1.5, num=100, endpoint=True
    )
    for τT in abatement_curve_τT:
        possibilities.add(
            # First best,
            AuditInfo(
                audit_rule="remote",
                audit_frac=0.0,
                τT=τT,
                detect_threshold=0.0,
                audit_cost=0.0,
            )
        )

    # Remote:
    for d in detect_opts:
        for τ in tau_opts:
            for T in T_opts:
                τT = τ * T
                if τT > SOCIAL_COST_METHANE_PER_KG * time_H:
                    continue
                possibilities.add(
                    AuditInfo(
                        audit_rule="remote",
                        audit_frac=0.0,
                        τT=τT,
                        detect_threshold=d,
                        audit_cost=0.0,
                    )
                )
    # Uniform, Target-X, and Target-e:
    for τ in tau_opts:
        for T in T_opts:
            for frac in audit_frac_opts:
                for cost in audit_cost_opts:
                    if (frac == 0 and cost == 0) or (frac > 0 and cost > 0):
                        # We're interested in either fixed budget or cost-based
                        # The two don't mix.
                        continue
                    for audit_rule in ("uniform", "target_x"):
                        possibilities.add(
                            AuditInfo(
                                audit_rule=audit_rule,
                                audit_frac=frac,
                                τT=τ * T,
                                detect_threshold=0.0,
                                audit_cost=cost,
                            )
                        )
                    for detect in detect_opts:
                        if cost > 0 and detect > 0:
                            # NOTE: this case isn't working yet.
                            continue
                        possibilities.add(
                            AuditInfo(
                                audit_rule="target_e",
                                audit_frac=frac,
                                τT=τ * T,
                                detect_threshold=detect,
                                audit_cost=cost,
                            )
                        )
    return possibilities


def all_tests():
    # Note: test_derivatives requires the user to check the output (see note
    # above. It's working for now, but if you make changes, read through the output)
    test_parse_period_wildcard()
    test_parse_audit_info()
    test_derivatives()
    test_audit_assignment()


if __name__ == "__main__":
    const = read_constants()
    TAU_LEVELS = const["TAU_LEVELS"]
    T_LEVELS = const["T_LEVELS"]
    SOCIAL_COST_METHANE_PER_KG = const["SOCIAL_COST_METHANE_PER_KG"]
    METHANE_GWP = const["METHANE_GWP"]
    all_tests()
