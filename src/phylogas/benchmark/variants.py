"""Assign benchmark variant labels to transmission components.

Moved from TwinSampler (``label_components.py``) because this produces
*ground truth for benchmarking*, not an ascertainment artefact: the labels are
matched against a real-world importation schedule so that prevalence
estimation has something to be scored against.

The split point is clean. TwinSampler still finds connected components and
emits ``component_id`` -- that is a property of the transmission graph it
builds. PhyloGAS consumes those component ids and decides which variant each
one represents.

Note that this label is deliberately independent of the genomes the painter
assigns. ``variant_benchmark`` is an epidemiological assignment derived from
the importation schedule; a sequence's actual lineage derives from whichever
seed founded its chain. Where both exist their disagreement is itself
measurable.

The two matchers below are ported verbatim from TwinSampler so that results
remain comparable; only the surrounding plumbing is new.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.preprocessing import minmax_scale


# ==========================================================================
# Matchers (ported verbatim from TwinSampler label_components.py)
# ==========================================================================
def mode1_temporal_match(sim_components, real_imports):
    """Assigns variants by finding the closest real importation in time."""
    print("  Applying Mode 1: Temporal & Proportional Matching...")
    if real_imports.empty or sim_components.empty:
        return {}

    # Sort both by tick
    sim_components = sim_components.sort_values('first_tick').reset_index(drop=True)
    real_imports = real_imports.sort_values('tick').reset_index(drop=True)
    
    assignments = {}
    real_imports_used = [False] * len(real_imports)
    
    # For each simulated component, find the best available real import
    for sim_idx, sim_row in sim_components.iterrows():
        best_real_idx = -1
        min_dist = float('inf')
        
        # Find the closest *unused* real import
        for real_idx, real_row in real_imports.iterrows():
            if not real_imports_used[real_idx]:
                dist = abs(sim_row['first_tick'] - real_row['tick'])
                if dist < min_dist:
                    min_dist = dist
                    best_real_idx = real_idx
        
        if best_real_idx != -1:
            assignments[sim_row['component_id']] = real_imports.loc[best_real_idx, 'variant']
            real_imports_used[best_real_idx] = True # Mark as used
            
    return assignments

def mode2_bipartite_match(sim_components, real_imports, time_weight=0.7, max_time_penalty_days=90):
    """
    Assigns variants using optimized bipartite matching on time and size,
    with a heavy penalty for matches outside a reasonable time window.
    (Optimized with NumPy broadcasting for massive component lists).
    """
    print(f"  Applying Mode 2: Bipartite Matching (time_weight={time_weight}, penalty_window={max_time_penalty_days} days)...")
    if real_imports.empty or sim_components.empty:
        return {}

    # Normalize size columns for fair comparison in cost function
    sim_components['norm_size'] = minmax_scale(sim_components['component_size'])
    real_imports['norm_size'] = minmax_scale(real_imports['sample_count'])
    
    # --- VECTORIZED MATRIX CALCULATION ---
    print("  Calculating cost matrix...")
    # 1. Extract values as numpy arrays and reshape for broadcasting
    # sim arrays become shape (N, 1)
    sim_ticks = sim_components['first_tick'].values[:, np.newaxis]
    sim_sizes = sim_components['norm_size'].values[:, np.newaxis]
    
    # real arrays become shape (1, M)
    real_ticks = real_imports['tick'].values[np.newaxis, :]
    real_sizes = real_imports['norm_size'].values[np.newaxis, :]
    
    # 2. Compute differences (Broadcasting automatically creates N x M matrices)
    time_diffs = np.abs(sim_ticks - real_ticks)
    size_costs = np.abs(sim_sizes - real_sizes)
    
    # 3. Calculate Base Cost
    time_normalizer = float(max_time_penalty_days)
    norm_time_cost = time_diffs / time_normalizer
    cost_matrix = (time_weight * norm_time_cost) + ((1 - time_weight) * size_costs)
    
    # 4. Apply massive penalty if out of window
    cost_matrix = np.where(time_diffs > max_time_penalty_days, cost_matrix + 1000, cost_matrix)
            
    # --- SOLVE THE ASSIGNMENT ---
    print(f"  Solving linear sum assignment for {cost_matrix.shape[0]}x{cost_matrix.shape[1]} matrix...")
    sim_indices, real_indices = linear_sum_assignment(cost_matrix)
    
    # --- CREATE THE ASSIGNMENT MAP ---
    assignments = {}
    
    # Extract native arrays for fast lookup
    comp_ids = sim_components['component_id'].values
    variants = real_imports['variant'].values
    
    for i, j in zip(sim_indices, real_indices):
        assignments[comp_ids[i]] = variants[j]
        
    return assignments


# ==========================================================================
# Plumbing
# ==========================================================================
VARIANT_COLUMN = "variant_benchmark"
LEGACY_VARIANT_COLUMN = "variant_label"

_TICK_COLUMNS = ("exposure_tick", "sim_tick", "tick")
_PID_COLUMNS = ("pid", "sim_pid", "alias_pid")


def summarize_components(events: pd.DataFrame) -> pd.DataFrame:
    """Rebuild the (component_id, first_tick, component_size) table.

    TwinSampler already ran connected-component detection, so this is a
    groupby rather than a second graph traversal.
    """
    if "component_id" not in events.columns:
        sys.exit(
            "  Error: no 'component_id' column in the events file.\n"
            "         Run TwinSampler's simulate-linelist first; it emits this."
        )
    tick_col = next((c for c in _TICK_COLUMNS if c in events.columns), None)
    if tick_col is None:
        sys.exit(f"  Error: no tick column; looked for {', '.join(_TICK_COLUMNS)}")
    pid_col = next((c for c in _PID_COLUMNS if c in events.columns), None)
    if pid_col is None:
        sys.exit(f"  Error: no pid column; looked for {', '.join(_PID_COLUMNS)}")

    summary = (
        events.groupby("component_id")
        .agg(first_tick=(tick_col, "min"), component_size=(pid_col, "nunique"))
        .reset_index()
    )
    print(f"  {len(summary):,} components (tick column '{tick_col}', "
          f"size by unique '{pid_col}')")
    return summary


def unroll_schedule(schedule: pd.DataFrame) -> pd.DataFrame:
    """Expand a per-row importation schedule into one row per cluster."""
    required = {"tick", "variant", "clusters", "sample_count"}
    missing = required - set(schedule.columns)
    if missing:
        sys.exit(f"  Error: schedule is missing columns: {', '.join(sorted(missing))}")
    rows = []
    for _, r in schedule.iterrows():
        n = int(r["clusters"])
        avg = r["sample_count"] / n if n > 0 else 0
        for _ in range(n):
            rows.append({"tick": r["tick"], "variant": r["variant"], "sample_count": avg})
    return pd.DataFrame(rows)


def assign_variants(
    events: pd.DataFrame,
    schedule: pd.DataFrame,
    mode: str = "bipartite",
    column: str = VARIANT_COLUMN,
) -> pd.DataFrame:
    """Attach ``column`` to every event, via its component's matched variant."""
    summary = summarize_components(events)
    real_imports = unroll_schedule(schedule)
    if real_imports.empty:
        sys.exit("  Error: the importation schedule expanded to zero clusters.")
    print(f"  {len(real_imports):,} scheduled importation clusters")

    if mode in ("bipartite", "variant_bipartite"):
        mapping = mode2_bipartite_match(summary.copy(), real_imports.copy())
    elif mode in ("temporal", "variant_temporal"):
        mapping = mode1_temporal_match(summary.copy(), real_imports.copy())
    else:
        sys.exit(f"  Error: unknown mode '{mode}' (use bipartite or temporal)")

    out = events.copy()
    out[column] = out["component_id"].map(mapping).fillna("background")
    counts = out[column].value_counts()
    print(f"  assigned {len(mapping):,} components; event counts:")
    for v, n in counts.head(12).items():
        print(f"    {str(v):22s} {n:>10,}")
    return out
