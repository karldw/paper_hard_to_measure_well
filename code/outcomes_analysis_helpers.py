# outcomes_analysis_derivative_functions.py
# (major revision #2)
# Imported by outcomes_analysis.py

import numpy as np
import scipy.sparse

from collections import namedtuple
import json


AuditInfo = namedtuple(
    "AuditInfo",
    ["audit_rule", "audit_frac", "τT", "detect_threshold", "audit_cost"]
    # Expected types: string,   float,        float,  float
)

# Container to pass data around
DataParam = namedtuple(
    "DataParam",
    [
        "price_H",  # price * H, units of dollar-hour per kg. Array with shape (N,)
        "e_size_expect",  # expected leak size, units of kg per hr. Array with shape (N,)
        "e_size_draw",  # draw of leak size, units of kg per hr. Array with shape (N,)
        "cost_coef",  # expected leak size, units of kg per hr. Array with shape (N, 2)
        "time_H",  # Time period H the well uses. Float.
        "prob_is_large",  # probability leak size is over detection threshold. Shape (N,)
        "gas_avg_mcfd",  # well pad's average daily production for the month observed. Shape (N,)
    ],
)

OutcomesOneDraw = namedtuple(
    "OutcomesOneDraw",
    [
        "dwl_mean",
        "dwl_tot",
        "emis_mean",
        "emis_tot",
        "tot_cost_per_kg_mean",
        "fee_per_kg_mean",
        "fee_per_kg_med",
        "fee_per_kg_p10",
        "fee_per_kg_p90",
        "net_private_cost_per_mcf_pct_price",
        "shadow_price",
        "audit_rule",
        "audit_frac",
        "τT",
        "detect_threshold",
        "audit_cost",
    ],
)


def dwl_per_well(t_effective, data_param):
    """
    Deadweight loss from setting t < δ

    Note that this function is used for scoring, but not for targeting audits,
    so we use draw_param.e_size_draw.
    """
    A = data_param.cost_coef[..., 0]
    α = data_param.cost_coef[..., 1]
    recip_α = 1 / α
    recip_α_plus1 = 1 / (α + 1)
    recip_A = 1 / A
    e_pH_τT_r = data_param.e_size_draw * (data_param.price_H + t_effective)
    e_p_δ_H = data_param.e_size_draw * (
        data_param.price_H + data_param.time_H * SOCIAL_COST_METHANE_PER_KG
    )
    dwl = ((e_p_δ_H - e_pH_τT_r * recip_α_plus1) * (e_pH_τT_r * recip_A) ** recip_α) - (
        α * recip_α_plus1 * e_p_δ_H * (e_p_δ_H * recip_A) ** recip_α
    )
    return dwl


def dwl_terms_for_objective(self, r):
    """
    Objective function, omitting the parts that don't depend on r.

    Note, despite the name `self`, this isn't actually a class method, but we'll
    use a class method to call it.
    r may have shape (N,) or (1,)
    """
    dwl_terms = (
        self.e_p_δ_H - (self.e * self.recip_α_plus1 * (self.pH + self.τT * r))
    ) * (self.e * (self.pH + self.τT * r) * self.recip_A) ** self.recip_α
    return dwl_terms


def objective_dwl(self, r):
    """
    dwl_terms_for_objective is a separate function so we can take the unweighted
    sum of DWL terms here and the weighted sum later
    """
    return np.sum(dwl_terms_for_objective(self, r))


def dwl_diff1(self, r):
    diff1 = (
        self.recip_A ** (self.recip_α)
        * (self.e) ** (1.0 + self.recip_α)
        * self.τT
        * (self.δH - self.τT * r)
        * (self.pH + self.τT * r) ** (self.recip_α - 1.0)
        * self.recip_α
    )
    return diff1


def dwl_diff2(self, r):
    diff2 = (
        self.recip_A ** (self.recip_α)
        * self.e ** (1.0 + self.recip_α)
        * self.τT ** 2
        * (self.pH + self.τT * r) ** (self.recip_α - 2.0)
        * (-self.α * (self.δH + self.pH) + self.δH - self.τT * r)
        * self.recip_α ** 2
    )
    return diff2


def dwl_uniform_incl_cost(self, r):
    # Just like dwl_target_x_incl_cost, except multiply r by N instead of summing
    return objective_dwl(self, r) + self.audit_cost * self.N * r


def dwl_uniform_incl_cost_diff1(self, r):
    # Just like dwl_target_x_incl_cost_diff1 except sum of terms.
    return np.sum(dwl_diff1(self, r) + self.audit_cost)


def dwl_target_x_incl_cost(self, r):
    return objective_dwl(self, r) + self.audit_cost * cost_target_x(self, r)


def dwl_target_x_incl_cost_diff1(self, r):
    # cost_target_x_diff1 is all 1, so just rely on numpy broadcasting here.
    return dwl_diff1(self, r) + self.audit_cost


def dwl_target_x_incl_cost_diff2(self, r):
    # cost_target_x_diff2 is all 0, so just drop the cost term.
    return dwl_diff2(self, r)


def dwl_target_e_low_detect_incl_cost(self, r):
    return objective_dwl(self, r) + self.audit_cost * cost_target_e_low_detect(self, r)


def dwl_target_e_low_detect_incl_cost_diff1(self, r):
    return dwl_diff1(self, r) + self.audit_cost * prob_leak_times_r_diff1(self, r)


def dwl_target_e_low_detect_incl_cost_diff2(self, r):
    return dwl_diff2(self, r) + self.audit_cost * prob_leak_times_r_diff2(self, r)


def cost_target_e_low_detect(self, r):
    return np.sum(prob_leak_times_r(self, r))


def cost_target_x(self, r):
    return np.sum(r)


def prob_leak_times_r(self, r):
    """
    Value of (1 - q(r)) * r
    (= r * (e * (p * H + τT * r) / A) ** (1 / α))
    """
    return r * (self.e * (self.pH + self.τT * r) * self.recip_A) ** (self.recip_α)


def prob_leak_times_r_diff1(self, r):
    r"""
    Value of d [(1 - q(r)) * r] / dr

    A^{- \frac{1}{\alpha}}
    e^{\frac{1}{\alpha}}
    (p H + r \tau T)^{\frac{1}{\alpha} - 1}
    (\alpha (p H + r \tau T) + r \tau T)
    \frac{1}{\alpha}
    """
    diff1 = (
        self.recip_A ** (self.recip_α)
        * self.e ** self.recip_α
        * (self.pH + self.τT * r) ** (self.recip_α - 1)
        * (self.α * (self.pH + self.τT * r) + self.τT * r)
        * self.recip_α
    )
    return diff1


def prob_leak_times_r_diff2(self, r):
    r"""
    Value of d^2 [(1 - q(r)) * r] / dr^2

    \frac{A^{- \frac{1}{\alpha}}
    e^{\frac{1}{\alpha}}
    \tau
    (p H + r \tau T)^{\frac{1}{\alpha} - 2}
    (2 \alpha p H + \alpha r \tau T + r \tau T)
    \frac{1}{\alpha^2}
    """
    diff2 = (
        self.recip_A ** (self.recip_α)
        * self.e ** self.recip_α
        * self.τT
        * (self.pH + self.τT * r) ** (self.recip_α - 2.0)
        * (self.α * (2 * self.pH + self.τT * r) + self.τT * r)
        * self.recip_α ** 2
    )
    return diff2


def prob_no_leak(self, r):
    r"""
    Value of q(r) = 1 - (e * (p * H + τ * T * r) / A) ** (1 / α)
    """
    prob = 1 - (self.e * (self.pH + self.τT * r) * self.recip_A) ** self.recip_α
    return prob


def prob_no_leak_diff1(self, r):
    r"""
    Value of d q(r) / dr

    q(r) = 1 - (e * (p * H + τ * T * r) / A) ** (1 / α)

    d q(r) / dr=
        -(1 / A) ** (1 / α) * T * τ * e ** (1 / α) *
        (H * p + T * τ * r) ** (1 / α - 1) / (α)

    """
    diff1 = -(
        self.recip_A ** (self.recip_α)
        * self.τT
        * self.e ** (self.recip_α)
        * (self.pH + self.τT * r) ** (self.recip_α - 1)
        * self.recip_α
    )
    return diff1


def prob_no_leak_diff2(self, r):
    r"""
    Value of d^2 q(r) / dr^2
    Where q(r) = 1 - (e * (p * H + τ * T * r) / A) ** (1 / α)
    second diff = = A ** (-1 / α) * (T * τ) ** 2 * e ** (1 / α)
        * (α - 1) * (H * p + T * r * τ) ** (1 / α) /
        (α ** 2 * (H * p + T * r * τ) ** 2)
    """
    diff2 = (
        self.recip_A ** (self.recip_α)
        * self.τT ** 2
        * self.e ** (self.recip_α)
        * (self.α - 1)
        * (self.pH + self.τT * r) ** (self.recip_α - 2)
        * self.recip_α ** 2
    )
    return diff2


def cost_target_x_diff1(self, r):
    r"""
    The first derivative of the target-x budget term wrt r_i
    """
    return np.ones_like(r)


def cost_target_x_diff2(self, r):
    r"""
    The second derivative of the target-x budget term wrt r_i
    """
    return np.zeros_like(r)


def hessian_target_e_low_detect(self, r, lagrange, obj_factor):
    # The callback for calculating the Hessian (second deriv) of the lagrangian
    # lagrange is only one number because only one constraint.
    # fmt: off
    hess_diag = (
        obj_factor * dwl_diff2(self, r)
        + lagrange * prob_leak_times_r_diff2(self, r)
    )
    # fmt: on
    return hess_diag


def hessian_target_x(self, r, lagrange, obj_factor):
    # The callback for calculating the Hessian (second deriv) of the lagrangian
    # `lagrange` doesn't factor in because second deriv of target-x budget is 0
    hess_diag = obj_factor * dwl_diff2(self, r)
    return hess_diag


# Same for incl_cost case:
hessian_target_x_incl_cost = hessian_target_x


def hessian_uniform_incl_cost(self, r, lagrange, obj_factor):
    # Here r has length 1; the hessian is a 1x1 matrix
    assert len(lagrange) == 0
    hess = obj_factor * np.sum(dwl_diff2(self, r))
    return hess


def hessian_target_e_incl_cost(self, r, lagrange, obj_factor):
    assert len(lagrange) == 0
    hess_diag = obj_factor * dwl_target_e_low_detect_incl_cost_diff2(self, r)
    return hess_diag


def no_constraint(self, x):
    """Provide a placeholder constraint function. Always zero."""
    return 0.0


def no_jacobian(self, x):
    """Provide a placeholder constraint function. Always zero."""
    return np.zeros_like(x)


def cost_target_e_high_detect(self, r_big_r_small):
    r_big, r_small = _divide_up_r(self.N, r_big_r_small)
    terms = _probability_weight(
        self.prob_is_large,
        prob_leak_times_r(self, r_big) + prob_no_leak(self, r_big) * r_small,
        r_small,
    )
    return np.sum(terms)


def cost_target_e_high_detect_diff1(self, r_big_r_small):
    r"""
    Budget term is: (1 - q) * r + q * rho
    = r * (e * (p * H + τT * r) / A) ** (1 / α)
      + rho * (1 - (e * (p * H + τT * r) / A) ** (1 / α))

    The z-weighted deriv with respect to r (r_big) is:

    z * A^{- \frac{1}{\alpha}}
    e^{\frac{1}{\alpha}}
    (p H + r \tau T)^{\frac{1}{\alpha} - 1}
    (\alpha (p H + r \tau T) + r \tau T - \rho \tau T)
    \frac{1}{\alpha}

    The z-weighted deriv with respect to \rho (r_small) is:
    1 - A^{- \frac{1}{\alpha}} z \left(H e p + T e r \tau\right)^{\frac{1}{\alpha}}
    """
    r_big, r_small = _divide_up_r(self.N, r_big_r_small)

    diff1_wrt_r_big = self.prob_is_large * (
        prob_leak_times_r_diff1(self, r_big) + prob_no_leak_diff1(self, r_big) * r_small
    )
    diff1_wrt_r_small = self.prob_is_large * prob_no_leak(self, r_big) + (
        1 - self.prob_is_large
    )
    if len(r_big_r_small) == 2:
        diff1_wrt_r_big = np.sum(diff1_wrt_r_big)
        diff1_wrt_r_small = np.sum(diff1_wrt_r_small)

    cost_diff1 = np.concatenate((diff1_wrt_r_big, diff1_wrt_r_small), axis=None)
    assert cost_diff1.shape == r_big_r_small.shape
    return cost_diff1


def hessian_target_e_high_detect(self, r_big_r_small, lagrange, obj_factor):
    r"""
    The hessian is a little tricky. r_big_r_small has length 2 or 2*N for the
    r and rho compondents, so we need to describe the problem to fit ipopt's
    understanding.

    The cross-partial matrix we need to describe is like this:
          ∂r_1 ... ∂r_N ∂rho_1 ... ∂rho_N
    ∂ r_1              |
    ⋮        A        |        B'
    ∂r_N               |
    --------------------------------------
    ∂rho_1             |
    ⋮        B        |        C
    ∂rho_N             |

    Block A: ∂^2 L / ∂^2 r_i^2
    i != j elements are zero

    Blocks B: ∂^2 L / ∂ r_i ∂ \rho_i
    i != j elements are zero


    Block C: ∂^2 L / ∂^2 \rho_i^2
    i != j elements are zero
    """
    r_big, r_small = _divide_up_r(self.N, r_big_r_small)
    # Create the values labeled A, B, and C above
    block_A = obj_factor * self.prob_is_large * dwl_diff2(
        self, r_big
    ) + lagrange * self.prob_is_large * (
        prob_leak_times_r_diff2(self, r_big) + prob_no_leak_diff2(self, r_big) * r_small
    )
    block_B = lagrange * self.prob_is_large * prob_no_leak_diff1(self, r_big)
    block_C = obj_factor * (1 - self.prob_is_large) * dwl_diff2(self, r_small)
    if self.num_vars == 2:
        block_A = np.sum(block_A)
        block_B = np.sum(block_B)
        block_C = np.sum(block_C)
    # Expected output is a 1d array that will be interpreted as the COO matrix
    # values. The order of the output here need to match the order in
    # hessianstructure_target_e_high.
    return np.concatenate([block_A, block_B, block_C], axis=None)


def _divide_up_r(N, r_big_r_small):
    """
    Helper function used in dwl_target_e_high_detect and
    dwl_target_e_high_detect_diff1
    """
    if r_big_r_small.shape == (2,):
        # Maybe avoid making these temp arrays?
        r_big = np.full((N,), r_big_r_small[0])
        r_small = np.full((N,), r_big_r_small[1])
    elif r_big_r_small.shape == (2 * N,):
        r_big = r_big_r_small[:N]
        r_small = r_big_r_small[N:]
    else:
        raise ValueError(f"Bad shape for r_big_r_small: {r_big_r_small.shape}")
    return r_big, r_small


def _probability_weight(prob, x, y):
    return prob * x + (1 - prob) * y


def dwl_target_e_high_detect(self, r_big_r_small):
    r_big, r_small = _divide_up_r(self.N, r_big_r_small)
    dwl = np.sum(
        _probability_weight(
            self.prob_is_large,
            dwl_terms_for_objective(self, r_big),
            dwl_terms_for_objective(self, r_small),
        )
    )
    return dwl


def dwl_target_e_high_detect_diff1(self, r_big_r_small):
    """
    Note: code in this function works with both the uniform and differentiated
    cases by taking the sum if necessary
    """
    r_big, r_small = _divide_up_r(self.N, r_big_r_small)
    diff1_big = dwl_diff1(self, r_big) * self.prob_is_large
    diff1_small = dwl_diff1(self, r_small) * (1 - self.prob_is_large)
    if len(r_big_r_small) == 2:
        # One value for r_big, one value for r_small
        # _divide_up_r already checks that r_big_r_small either has shape (2,)
        # or (2 * N,)
        diff1_big = np.sum(diff1_big)
        diff1_small = np.sum(diff1_small)

    diff1 = np.concatenate((diff1_big, diff1_small), axis=None)
    assert diff1.shape == r_big_r_small.shape
    return diff1


def dwl_target_e_high_detect_incl_cost(self, r):
    dwl = dwl_target_e_high_detect(self, r)
    cost = self.audit_cost * cost_target_e_high_detect(self, r)
    return np.sum(dwl + cost)


def dwl_target_e_high_detect_incl_cost_diff1(self, r):
    dwl_diff1 = dwl_target_e_high_detect_diff1(self, r)
    cost_diff1 = self.audit_cost * cost_target_e_high_detect_diff1(self, r)
    return dwl_diff1 + cost_diff1


def hessian_target_e_high_detect_incl_cost(self, r_big_r_small, lagrange, obj_factor):
    """See comments in hessian_target_e_high_detect
    Same idea, but now we're multiplying by audit_cost instead of lagrange.
    """
    r_big, r_small = _divide_up_r(self.N, r_big_r_small)
    # Block A:
    block_A = obj_factor * self.prob_is_large * dwl_diff2(
        self, r_big
    ) + self.audit_cost * self.prob_is_large * (
        prob_leak_times_r_diff2(self, r_big) + prob_no_leak_diff2(self, r_big) * r_small
    )
    # Blocks B:
    block_B = self.audit_cost * self.prob_is_large * prob_no_leak_diff1(self, r_big)
    # Block C:
    block_C = obj_factor * (1 - self.prob_is_large) * dwl_diff2(self, r_small)
    if self.num_vars == 2:
        block_A = np.sum(block_A)
        block_B = np.sum(block_B)
        block_C = np.sum(block_C)
    # Expected output is a 1d array that will be interpreted as the COO matrix
    # values. The order of the output here need to match the order in
    # hessianstructure_target_e_high.
    return np.concatenate([block_A, block_B, block_C], axis=None)


def hessianstructure_default(self):
    """Define the structure of the hessian (2nd deriv matrix) of the
    lagrangian with respect to r_i, r_j. It's diagonal for almost all cases.
    """
    hs = scipy.sparse.eye(self.num_vars, format="coo")
    return (hs.col, hs.row)


def hessianstructure_target_e_high(self):
    """Define the structure of the hessian (2nd deriv matrix) of the
    lagrangian with respect to i, j. It's no longer diagonal because we're using
    a double-length choice variable r_big_r_small.
    Instead each quadrent of the matrix is diagonal.
    (In the case of uniform r and rho; we end up with a dense 2x2)

    Note that derivatives are symmetric, and ipopt knows this, so we only have
    to deal with the lower triangle.
    """
    assert self.num_vars in {2 * self.N, 2}
    diag = scipy.sparse.identity(self.num_vars // 2, dtype="i", format="coo")
    hs = scipy.sparse.bmat([[diag, None], [diag, diag]], format="coo")
    # The indexes that come out here are in the order written (rows, then cols)
    # In other words, sub-matrices A, B, C
    return (hs.col, hs.row)


class ProblemSetup(object):
    def __init__(self, data_param, audit_info):
        """
        Set up the problem, with a budget constraint.
        Supports target_x and target_e cases.
        """

        # Define and pre-compute some parameters that don't depend on r.
        self.e = data_param.e_size_expect
        self.pH = data_param.price_H
        self.τT = audit_info.τT
        self.A = data_param.cost_coef[..., 0]
        self.α = data_param.cost_coef[..., 1]
        self.recip_α = 1.0 / self.α
        self.recip_α_plus1 = 1.0 / (self.α + 1.0)
        self.recip_A = 1.0 / self.A
        self.δH = data_param.time_H * SOCIAL_COST_METHANE_PER_KG
        self.e_p_δ_H = self.e * (self.pH + self.δH)
        self.audit_cost = audit_info.audit_cost
        self.N = len(self.e)
        self.prob_is_large = data_param.prob_is_large

        if audit_info.audit_frac > 0 and audit_info.audit_cost == 0:
            if audit_info.audit_rule == "target_x":
                self.fn_objective = objective_dwl
                self.fn_gradient = dwl_diff1
                self.fn_constraints = cost_target_x
                self.fn_jacobian = cost_target_x_diff1
                self.fn_hessian = hessian_target_x
                self.fn_hessianstructure = hessianstructure_default
                self.num_vars = self.N
            elif audit_info.audit_rule == "target_e":
                if audit_info.detect_threshold == 0:
                    self.fn_objective = objective_dwl
                    self.fn_gradient = dwl_diff1
                    self.fn_constraints = cost_target_e_low_detect
                    self.fn_jacobian = prob_leak_times_r_diff1
                    self.fn_hessian = hessian_target_e_low_detect
                    self.fn_hessianstructure = hessianstructure_default
                    self.num_vars = self.N
                else:
                    self.fn_objective = dwl_target_e_high_detect
                    self.fn_gradient = dwl_target_e_high_detect_diff1
                    self.fn_constraints = cost_target_e_high_detect
                    self.fn_jacobian = cost_target_e_high_detect_diff1
                    self.fn_hessian = hessian_target_e_high_detect
                    self.fn_hessianstructure = hessianstructure_target_e_high
                    self.num_vars = self.N * 2
            else:
                raise ValueError(f"Bad audit_rule: {audit_info.audit_rule}")
        elif audit_info.audit_frac == 0 and audit_info.audit_cost > 0:
            self.fn_constraints = no_constraint
            self.fn_jacobian = no_jacobian
            if audit_info.audit_rule == "uniform":
                self.fn_objective = dwl_uniform_incl_cost
                self.fn_gradient = dwl_uniform_incl_cost_diff1
                self.fn_hessian = hessian_uniform_incl_cost
                self.fn_hessianstructure = hessianstructure_default
                self.num_vars = 1
            elif audit_info.audit_rule == "target_x":
                self.fn_objective = dwl_target_x_incl_cost
                self.fn_gradient = dwl_target_x_incl_cost_diff1
                self.fn_hessian = hessian_target_x_incl_cost
                self.fn_hessianstructure = hessianstructure_default
                self.num_vars = self.N
            elif audit_info.audit_rule == "target_e":
                if audit_info.detect_threshold == 0:
                    self.fn_objective = dwl_target_e_low_detect_incl_cost
                    self.fn_gradient = dwl_target_e_low_detect_incl_cost_diff1
                    self.fn_hessian = hessian_target_e_incl_cost
                    self.fn_hessianstructure = hessianstructure_default
                    self.num_vars = self.N
                else:
                    self.fn_objective = dwl_target_e_high_detect_incl_cost
                    self.fn_gradient = dwl_target_e_high_detect_incl_cost_diff1
                    self.fn_hessian = hessian_target_e_high_detect_incl_cost
                    self.fn_hessianstructure = hessianstructure_target_e_high
                    self.num_vars = self.N * 2
            else:
                raise ValueError(f"Bad audit_rule: {audit_info.audit_rule}")
        else:
            raise ValueError(
                f"Exactly one of audit_frac and audit_cost should be positive"
            )

    # Define the functions like this so we get proper bound methods (and
    # access to `self`), but allowing the actual function to vary, and
    # without mucking around with MethodType
    def objective(self, x):
        return self.fn_objective(self, x)

    def gradient(self, x):
        return self.fn_gradient(self, x)

    def constraints(self, x):
        return self.fn_constraints(self, x)

    def jacobian(self, x):
        return self.fn_jacobian(self, x)

    def hessian(self, x, lagrange, obj_factor):
        return self.fn_hessian(self, x, lagrange, obj_factor)

    def hessianstructure(self):
        return self.fn_hessianstructure(self)


# Some functions we use to calculate the uniform-budgeted case analytically


def prob_leak_with_policy(t_effective, cost_coef, e_size, price_H):
    """
    Probability of a leak (1 - q) under policy t_effective

    Note that e_size may be either e_size_expect or e_size_draw, depending what
    we're doing with the function.
    """
    A = cost_coef[..., 0]
    α = cost_coef[..., 1]
    prob = (e_size * (price_H + t_effective) / A) ** (1 / α)
    assert is_probability(prob)
    return prob


def dwl_deriv1_wrt_r(r, data_param, τT):
    r"""
    The first derivative with respect to r of the DWL for one well
    Only used for solve_for_λ_uniform

    \frac{
        A^{- \frac{1}{\alpha}}
        e^{\frac{\alpha + 1}{\alpha}}
        \tau T
        (\delta H - \tau T r)
        (p H + \tau T r)^{\frac{1 - \alpha}{\alpha}}
    }{
        \alpha
    }
    """
    # cost_coef is a (1, 2) or (N, 2) array
    A = data_param.cost_coef[..., 0]
    α = data_param.cost_coef[..., 1]
    recip_α = 1 / α
    δ = SOCIAL_COST_METHANE_PER_KG
    H = data_param.time_H

    dwl_diff1 = (
        A ** (-recip_α)
        * (data_param.e_size_expect) ** ((α + 1) * recip_α)
        * τT
        * (δ * H - τT * r)
        * (data_param.price_H + τT * r) ** ((1 - α) * recip_α)
        * recip_α
    )
    return dwl_diff1


def is_budget_binding(audit_info, time_H):
    """
    If we can audit every well with prob (δ * H) / (τ * T), there's no
    benefit to auditing more. (and the converse)
    """
    return audit_info.audit_frac < (SOCIAL_COST_METHANE_PER_KG * time_H / audit_info.τT)


def calc_well_audit_prob_uniform(data_param, audit_info, r_guess):
    assert len(data_param.e_size_expect) > 1
    if not (audit_info.audit_frac > 0 and audit_info.audit_cost == 0):
        raise ValueError("Other cases should go to ipopt")
    if is_budget_binding(audit_info, data_param.time_H):
        r, λ = solve_for_λ_uniform(data_param, audit_info, r_guess)
    else:
        r = SOCIAL_COST_METHANE_PER_KG * data_param.time_H / audit_info.τT
        λ = 0.0
    return r, λ


def solve_for_λ_uniform(data_param, audit_info, r_guess):
    """
    Solve for λ for the uniform-audit policy.
    Returns resulting r and λ

    Uses the analytical result that when the budget is binding, r = audit_frac,
    and λ = mean( d DWL / dr )

    """
    N = len(data_param.e_size_expect)
    r = np.full((N,), audit_info.audit_frac)
    dwl_deriv_vals = dwl_deriv1_wrt_r(r, data_param, audit_info.τT)
    λ = np.mean(dwl_deriv_vals)
    return audit_info.audit_frac, λ


def is_probability(x):
    """Check if `x` is in [0, 1]"""
    return np.all((x >= 0.0) & (x <= 1.0))


def read_constants(json_file="code/constants.json"):
    """
    Read some important constants from `filename`.

    Notes on their values:
    * SOCIAL_COST_METHANE_PER_KG is $2, which is about $51.55 / ton  CO2e. That's low.
    * METHANE_GWP is 29.8, from IPCC AR6
    * TAU_LEVELS are 2 * SOCIAL_COST_METHANE_PER_KG, and $5 per ton CO2e
    * T_LEVELS are in hours, so 1day is 24, 1week is 168, 1month is 730, 3month is 2190, and 1year is 8760.
        * Not all of these are used. If defining additional levels, make sure the cross-product of tau_levels and t_levels is still unique (e.g. can't add 2month because high * 1month == med * 2month)

        If we were defining in python instead of reading json:
        TAU_LEVELS = {
            "high": 2 * SOCIAL_COST_METHANE_PER_KG,
            "med": SOCIAL_COST_METHANE_PER_KG,
            "low": 5 / 1000 * METHANE_GWP,  # low = $5/ton CO2e
        }
        T_LEVELS = {
            "1day": 24.0,  # 1 day
            "1week": 168.0,  # 1 week
            "1month": 730.0,  # 1 month
            "3month": 2190.0,  # 3 months
            "1year": 8760.0,  # 1 year
        }

    * LEAK_SIZE_DEF is 5.0 kg/hr, based on the detection threshold of AVIRIS-NG (*not* the policy detection threshold)
    * AUDIT_COST_TO_COMPARE_FIXED_VS_OPTIM_DWL is $600/audit
    * N_WELL_PADS is the number of well pads in this data. It's checked in some of the code, but included here for convenience.
    * MODEL_NAMES known model names, categorized.
        * note that the "normal" models still get a lognormal treatment
    """
    with open(json_file, "rt") as f:
        const = json.load(f)
    expected_names = {
        "SOCIAL_COST_METHANE_PER_KG",
        "METHANE_GWP",
        "LEAK_SIZE_DEF",
        "AUDIT_COST_TO_COMPARE_FIXED_VS_OPTIM_DWL",
        "N_WELL_PADS",
        "TAU_LEVELS",
        "T_LEVELS",
        "MODEL_NAMES",
        "POLICY_DETECT_THRESHOLD_HIGH",
        "POLICY_DETECT_THRESHOLD_LOW",
    }
    assert expected_names == set(const.keys())
    return const


def abatement_cost_per_pad(prob_leak, data_param):
    """
    Abatement cost for the well.

    Not used in the optimization, since we have analytical expressions for DWL.
    """
    A = data_param.cost_coef[..., 0]
    α = data_param.cost_coef[..., 1]
    α_plus1 = α + 1
    abatement_cost = -(A / α_plus1) * (prob_leak ** α_plus1)

    return abatement_cost


SOCIAL_COST_METHANE_PER_KG = read_constants()["SOCIAL_COST_METHANE_PER_KG"]


if __name__ == "__main__":
    raise RuntimeError("This file isn't meant to be run directly.")
