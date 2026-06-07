"""
Quadratic programming solver for amino acid profile optimization.

Problem formulation
-------------------
Given n foods with weight fractions x_i (x_i >= 0, sum = 1):

    blend_protein   = sum_i x_i * p_i          [g protein per 100g blend]
    blend_aa_j      = sum_i x_i * aa_ij         [g AA j per 100g blend]
    blend_norm_j    = blend_aa_j / blend_protein [g AA j per g protein]

Minimise  sum_j max(0, target_j - blend_norm_j)^2      [deficit-only]
s.t.      sum_i x_i = 1,  x_i >= 0  (optionally x_i <= max_i)

One-sided objective: only deficits are penalised. Excess AAs are
metabolised harmlessly, so overshooting a target does not cost anything.
A symmetric (blend - target)^2 objective would cause the solver to prefer
uniformly mediocre blends over blends that nail most AAs but exceed the
target on a few — the wrong behaviour for dietary optimisation.

The objective is a ratio of linear functions — not pure QP — but SLSQP
handles it without issue (smooth on the feasible simplex for typical
food databases).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .data import AA_COLS, AA_LABELS


@dataclass
class BlendResult:
    food_ids: list[str]
    food_names: list[str]
    weight_fractions: np.ndarray     # sum = 1
    food_grams: np.ndarray           # for given protein_target
    protein_per_100g_blend: float    # g protein per 100g blend
    blend_norm: np.ndarray           # g AA / g protein in blend
    target_norm: np.ndarray          # g AA / g protein target
    gap: np.ndarray                  # blend_norm - target_norm
    coverage: np.ndarray             # min(1, blend/target) per AA — 1.0 = fully met
    rmse: float
    limiting: list[tuple[str, float]]   # (aa_label, coverage%) for deficient AAs
    infeasible_aa: list[str]            # AAs with zero coverage across all foods


def optimize(
    food_ids: list[str],
    food_names: list[str],
    food_proteins: np.ndarray,   # g protein per 100g food, shape (n,)
    food_aa: np.ndarray,         # g AA per 100g food, shape (n, k)
    target_norm: np.ndarray,     # g AA per g protein, shape (k,)
    protein_target: float = 30.0,
    min_fractions: np.ndarray | None = None,
    max_fractions: np.ndarray | None = None,
    n_restarts: int = 8,
) -> BlendResult:
    n, k = food_aa.shape
    p = food_proteins  # (n,)

    # Detect AAs that are zero in ALL selected foods — optimizer can't fix these
    aa_present = food_aa.sum(axis=0) > 0  # (k,) bool
    infeasible = [AA_LABELS[AA_COLS[j]] for j in range(k) if not aa_present[j]]

    lb = min_fractions if min_fractions is not None else np.zeros(n)
    ub = max_fractions if max_fractions is not None else np.ones(n)
    bounds = list(zip(lb, ub))
    constraints = [
        {"type": "eq", "fun": lambda x: np.sum(x) - 1.0, "jac": lambda x: np.ones(n)}
    ]

    def objective(x: np.ndarray) -> float:
        blend_protein = float(x @ p)
        if blend_protein < 1e-10:
            return 1e12
        blend_aa = x @ food_aa          # (k,)
        norm = blend_aa / blend_protein
        # One-sided: only penalize deficits. Excess AAs are metabolised
        # harmlessly — penalising overshoot causes the solver to prefer
        # blends that are uniformly mediocre over blends that nail most
        # AAs but exceed the target on a few.
        deficit = np.maximum(0.0, target_norm - norm)
        return float(deficit @ deficit)

    best = None
    rng = np.random.default_rng(42)
    for _ in range(n_restarts):
        x0 = rng.dirichlet(np.ones(n))
        # Clip to bounds
        x0 = np.clip(x0, lb, ub)
        x0 /= x0.sum()
        res = minimize(
            objective,
            x0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"ftol": 1e-14, "maxiter": 2000},
        )
        if best is None or res.fun < best.fun:
            best = res

    x = np.clip(best.x, 0, 1)
    x /= x.sum()

    blend_protein = float(x @ p)
    blend_aa = x @ food_aa
    blend_norm = blend_aa / max(blend_protein, 1e-10)
    gap = blend_norm - target_norm

    # Coverage: how much of target is met (capped at 1.0 = 100%)
    with np.errstate(divide="ignore", invalid="ignore"):
        coverage = np.where(target_norm > 0, np.minimum(1.0, blend_norm / target_norm), 1.0)

    limiting = [
        (AA_LABELS[AA_COLS[j]], float(coverage[j]) * 100)
        for j in range(k)
        if coverage[j] < 1.0 - 1e-4
    ]
    limiting.sort(key=lambda x: x[1])

    # Food quantities for the protein target
    serving_weight = protein_target / max(blend_protein / 100, 1e-10)
    food_grams = x * serving_weight

    rmse = float(np.sqrt(best.fun / k))   # root mean squared deficit (zero = fully met)

    return BlendResult(
        food_ids=food_ids,
        food_names=food_names,
        weight_fractions=x,
        food_grams=food_grams,
        protein_per_100g_blend=blend_protein,
        blend_norm=blend_norm,
        target_norm=target_norm,
        gap=gap,
        coverage=coverage,
        rmse=rmse,
        limiting=limiting,
        infeasible_aa=infeasible,
    )
