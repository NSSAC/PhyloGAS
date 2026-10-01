"""Scoring of an inferred geographic transition matrix against ABM truth.

Extracted verbatim from BeyondBaseline's ``run_all_scenarios.py`` (the block
that previously ran inline at ~line 1024 inside the per-scenario loop). Pulled
out so it can be called on saved sample sets after the fact, instead of only
during a sampling sweep.

The four metrics, all computed on the off-diagonal entries of the normalised
transition matrices (the diagonal is self-transitions and carries no
geographic-flow signal):

``pearson_r``
    Linear correlation of transition intensities. Sensitive to getting the
    *magnitude* of flow right.
``masked_mae``
    Mean absolute error over edges active in either matrix. Masking avoids the
    huge shared-zero population flattering the score.
``cosine_similarity``
    Angle between the flow vectors: is the *shape* of the flow right,
    independent of overall scale.
``topological_f1``
    Binarised presence/absence of each edge: did we recover the right set of
    geographic connections at all, regardless of intensity.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
from sklearn.metrics import f1_score

from . import mugration


def score_offdiagonals(truth_flat, inferred_flat) -> dict:
    """Score two flattened off-diagonal flow vectors.

    Both arrays must already be aligned to the same alphabet ordering.
    """
    truth_flat = np.asarray(truth_flat, dtype=float)
    inferred_flat = np.asarray(inferred_flat, dtype=float)

    if truth_flat.shape != inferred_flat.shape:
        raise ValueError(
            f"shape mismatch: truth {truth_flat.shape} vs inferred {inferred_flat.shape}. "
            "Both matrices must be aligned to the same alphabet first."
        )

    # pearsonr is undefined when either input is constant (e.g. all zeros).
    if truth_flat.size < 2 or np.all(truth_flat == truth_flat[0]) or np.all(
        inferred_flat == inferred_flat[0]
    ):
        pearson_r = float("nan")
    else:
        pearson_r, _ = pearsonr(truth_flat, inferred_flat)

    active = (truth_flat > 0) | (inferred_flat > 0)
    masked_mae = (
        float(np.mean(np.abs(truth_flat[active] - inferred_flat[active])))
        if np.any(active) else 0.0
    )

    cos_sim = (
        1.0 - cosine(truth_flat, inferred_flat)
        if truth_flat.sum() > 0 and inferred_flat.sum() > 0 else 0.0
    )

    truth_bin = (truth_flat > 0).astype(int)
    inferred_bin = (inferred_flat > 0).astype(int)
    f1 = (
        f1_score(truth_bin, inferred_bin, zero_division=0)
        if truth_bin.sum() > 0 or inferred_bin.sum() > 0 else 1.0
    )

    return {
        "pearson_r": float(pearson_r),
        "masked_mae": float(masked_mae),
        "cosine_similarity": float(cos_sim),
        "topological_f1": float(f1),
        "n_edges_truth": int(truth_bin.sum()),
        "n_edges_inferred": int(inferred_bin.sum()),
    }


def score_matrices(truth_matrix, truth_alphabet, inferred_matrix, inferred_alphabet) -> dict:
    """Align two transition matrices to a common alphabet, then score them."""
    master = list(truth_alphabet)
    truth_norm = mugration.align_and_normalize_matrix(truth_matrix, truth_alphabet, master)
    inferred_norm = mugration.align_and_normalize_matrix(
        inferred_matrix, inferred_alphabet, master
    )
    truth_flat = np.array(mugration.get_off_diagonals(truth_norm, master))
    inferred_flat = np.array(mugration.get_off_diagonals(inferred_norm, master))
    return score_offdiagonals(truth_flat, inferred_flat)
