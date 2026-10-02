"""Benchmark metrics that require agent-based-model ground truth.

Moved from BeyondBaseline. The dividing line is ground truth: BeyondBaseline
must stay runnable by a health department on a real line list, so anything
needing the ABM's hidden transmission graph or true infection counts lives
here instead.

What moved:

===================================  ==================================
eval_type                            needs
===================================  ==================================
B_cumulative_infections              true infection counts
C_stride_window_infections           true infection counts
E_stride_variant_prevalence_error    true variant counts + variant_benchmark
F_stride_component_coverage          the transmission graph
M_8_week_rolling_tree_coverage       the transmission graph
*_coverage_size_*                    the transmission graph
===================================  ==================================

What stayed in BeyondBaseline: KL divergence against the line list and the
population, and the equity metrics. Those need only the line list the sampler
was handed, so a surveillance program can compute them on real data.

Functions are ported verbatim so results remain comparable with published
runs; only the plumbing around them is new.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .variants import VARIANT_COLUMN, LEGACY_VARIANT_COLUMN

DATE_FIELD_DEFAULT = "date"


def _kl_dist(p, q, eps: float = 1e-12) -> float:
    """KL(p || q) over the union of their indices."""
    idx = sorted(set(p.index) | set(q.index))
    pv = p.reindex(idx, fill_value=0.0).astype(float).values + eps
    qv = q.reindex(idx, fill_value=0.0).astype(float).values + eps
    pv = pv / pv.sum()
    qv = qv / qv.sum()
    return float(np.sum(pv * np.log(pv / qv)))


# Imported under its historical name by the ported functions below.
kl_dist = _kl_dist


def variant_column(df: pd.DataFrame) -> "str | None":
    """Return whichever variant column the frame carries, preferring the new one."""
    for c in (VARIANT_COLUMN, LEGACY_VARIANT_COLUMN):
        if c in df.columns:
            return c
    return None


# ==========================================================================
# Ported verbatim from BeyondBaseline run_all_scenarios.py
# ==========================================================================
def build_weekly_infections(infections_path, pop_df, start_date, num_weeks_ref, date_col: str = "date"):
    """
    Build weekly infections history aligned to linelist slicing.
    Now requires a real date column in the infections file (default: 'date').
    """
    # Let pandas sniff the delimiter (comma, tab, etc.) and avoid skipping header rows.
    inf = pd.read_csv(infections_path, sep=None, engine="python", dtype={'alias_pid': str, 'alias_contact': str, 'sim_pid': str, 'pid': str, 'contact_pid': str})
    inf.columns = [c.strip() for c in inf.columns]

    inf = normalize_age_group_col(inf, "age_group")
    if date_col not in inf.columns:
        # try case-insensitive match (e.g., 'Date', 'DATE')
        ci_map = {c.lower(): c for c in inf.columns}
        if date_col.lower() in ci_map:
            date_col = ci_map[date_col.lower()]
        else:
            raise ValueError(
                f"Infections file must include a '{date_col}' column "
                f"(case-insensitive). Found columns: {list(inf.columns)}"
            )

    # Parse dates
    inf[date_col] = pd.to_datetime(inf[date_col], errors="coerce")
    if inf[date_col].isna().all():
        raise ValueError(f"Unable to parse any dates in infections column '{date_col}'.")

    # Map pid -> group using population file
    pid_col = "sim_pid" if "sim_pid" in inf.columns else ("pid" if "pid" in inf.columns else None)
    if pid_col is None:
        raise ValueError("Infections file must contain 'sim_pid' or 'pid' to map to demographic groups.")

    pop_pid_col = "sim_pid" if "sim_pid" in pop_df.columns else "pid"
    if pop_pid_col not in pop_df.columns:
        raise ValueError("Population file must have a 'sim_pid' or 'pid' column to map infections to 'group'.")

    pid_group_map = pop_df[[pop_pid_col, "group"]].dropna()
    inf = inf.merge(pid_group_map, left_on=pid_col, right_on=pop_pid_col, how="left")
    inf = inf.dropna(subset=["group"])

    # Weekly counts aligned to linelist weeks
    weekly_inf_hist = []
    cur = start_date
    for _ in range(num_weeks_ref):
        prev_mon = cur - timedelta(days=7)
        prev_sun = cur - timedelta(days=1)
        mask = (inf[date_col] >= prev_mon) & (inf[date_col] <= prev_sun)
        weekly_inf_hist.append(inf.loc[mask, "group"].value_counts())
        cur += timedelta(weeks=1)

    return weekly_inf_hist, inf


def build_weekly_variant_counts(
    infections_path,
    start_date,
    num_weeks_ref,
    date_col: str = "date",
    variant_col: str = "variant_label",
):
    """
    Build weekly *true* variant counts from the infections file.

    Returns
    -------
    weekly_variant_counts : list of pd.Series
        One entry per week. Each Series has index=variant_label, values=counts.
    """
    inf = pd.read_csv(infections_path, sep=None, engine="python", dtype={'alias_pid': str, 'alias_contact': str, 'sim_pid': str, 'pid': str, 'contact_pid': str})
    inf.columns = [c.strip() for c in inf.columns]

    # resolve date column (case-insensitive)
    if date_col not in inf.columns:
        ci_map = {c.lower(): c for c in inf.columns}
        if date_col.lower() in ci_map:
            date_col = ci_map[date_col.lower()]
        else:
            raise ValueError(
                f"Infections file must include a '{date_col}' column "
                f"(case-insensitive). Found columns: {list(inf.columns)}"
            )

    if variant_col not in inf.columns:
        raise ValueError(
            f"Infections file must include a '{variant_col}' column for variant labels. "
            f"Found columns: {list(inf.columns)}"
        )

    inf[date_col] = pd.to_datetime(inf[date_col], errors="coerce")
    if inf[date_col].isna().all():
        raise ValueError(f"Unable to parse any dates in infections column '{date_col}'.")

    weekly_variant_counts = []
    cur = start_date
    for _ in range(num_weeks_ref):
        prev_mon = cur - timedelta(days=7)
        prev_sun = cur - timedelta(days=1)
        mask = (inf[date_col] >= prev_mon) & (inf[date_col] <= prev_sun)
        wk = inf.loc[mask]

        if wk.empty:
            weekly_variant_counts.append(pd.Series(dtype=float))
        else:
            counts = wk[variant_col].value_counts()
            weekly_variant_counts.append(counts.astype(float))

        cur += timedelta(weeks=1)

    return weekly_variant_counts, inf



# ----------------- evaluation helpers -----------------


def build_undirected_adj(df, pid_col="alias_pid", contact_col="alias_contact"):
    """
    Builds an adjacency list for the entire transmission network (undirected).
    Returns: dict {pid: [neighbor_pids]}
    """
    adj = {}
    
    # Ensure strings
    df[pid_col] = df[pid_col].astype(str)
    df[contact_col] = df[contact_col].astype(str)
    
    for _, row in df.iterrows():
        u = row[pid_col]
        v = row[contact_col]
        
        # Initialize
        if u not in adj: adj[u] = []
        
        # Valid edge check (ignore -1 or self-loops)
        if u and v and u != "-1" and v != "-1" and u != "nan" and v != "nan" and u != v:
            if v not in adj: adj[v] = []
            
            # Add undirected edge
            adj[u].append(v)
            adj[v].append(u)
            
    return adj


def calculate_coverage_score(target_population_set, sampled_set, adj_graph):
    """
    Computes Coverage Score = (1 / |Pt|) * Sum(1 / (d(u, S) + 1))
    using Multi-Source BFS.
    """
    if not target_population_set:
        return 0.0
    
    if not sampled_set:
        return 0.0 # d(u, S) is inf, 1/(inf+1) is 0

    # Multi-Source BFS Initialization
    queue = deque()
    distances = {} # Stores d(u, S)
    
    # Initialize with all sampled nodes that exist in the graph
    for s in sampled_set:
        if s in adj_graph: 
            distances[s] = 0
            # FIX: Append tuple (node, distance)
            queue.append((s, 0))
        # Note: If s is not in adj_graph (isolated), it doesn't help reach others, 
        # but it has distance 0 to itself. This is handled implicitly if s in target_population_set.
        # However, for the BFS to run, we only queue valid graph nodes.

    # BFS
    while queue:
        # FIX: Now this unpacks correctly
        current, dist = queue.popleft()
        
        # Explore neighbors
        if current in adj_graph:
            for neighbor in adj_graph[current]:
                if neighbor not in distances:
                    distances[neighbor] = dist + 1
                    # FIX: Append tuple (neighbor, new_distance)
                    queue.append((neighbor, dist + 1))
    
    # Calculate Score
    total_score = 0.0
    
    for u in target_population_set:
        # If u was visited, we have a distance.
        if u in distances:
            d = distances[u]
            total_score += 1.0 / (d + 1.0)
        # If u corresponds to a sampled node that was isolated (not in adj_graph),
        # its distance to S is 0 (since it IS in S).
        elif u in sampled_set:
             total_score += 1.0 # 1 / (0 + 1)
        else:
            # d(u, S) = infinity -> term is 0
            total_score += 0.0
            
    return total_score / len(target_population_set)


def precompute_component_sizes(adj_graph, all_pids):
    """
    Returns a dict {pid: component_size} for every pid in the graph.
    Uses BFS/DFS to find connected components.
    """
    pid_to_size = {}
    visited = set()
    
    # Ensure all pids are in the map, defaulting to size 1 if isolated/missing from graph
    # (Though adj_graph usually contains everyone if built from linelist)
    for pid in all_pids:
        if pid not in pid_to_size:
            pid_to_size[pid] = 1

    for start_node in adj_graph:
        if start_node not in visited:
            # Found a new component, traverse it to count size
            component_nodes = []
            queue = deque([start_node])
            visited.add(start_node)
            component_nodes.append(start_node)
            
            while queue:
                curr = queue.popleft()
                for nbr in adj_graph.get(curr, []):
                    if nbr not in visited:
                        visited.add(nbr)
                        queue.append(nbr)
                        component_nodes.append(nbr)
            
            # Assign size to all members
            size = len(component_nodes)
            for node in component_nodes:
                pid_to_size[node] = size
                
    return pid_to_size


# --------- shared helpers ---------
# Map short codes to long labels; keep existing long labels untouched.
AGE_GROUP_MAP = {
    "p": "Preschool (0-4)",
    "s": "Student (5-17)",
    "a": "Adult (18-49)",
    "o": "Older adult (50-64)",
    "g": "Senior (65+)",
}


def series_auc(ys, xs=None):
    """Trapezoidal AUC over actual week positions; ignores NaNs. Lower is better."""
    y = np.asarray(list(ys), dtype=float)
    if xs is None:
        x = np.arange(1, len(y) + 1, dtype=float)
    else:
        x = np.asarray(list(xs), dtype=float)
        if len(x) != len(y):
            raise ValueError("series_auc requires xs and ys to have the same length.")
    m = np.isfinite(y)
    if m.sum() < 2:
        return float("nan")
    trapz_fn = getattr(np, "trapezoid", np.trapz)
    return float(trapz_fn(y[m], x[m]))

# SCEN_LABELS = {
#     1: "CS-C(LL)", 2: "RS-R(LL)", 3: "RS-C(LL)",
#     4: "CS-C(LL,P)", 5: "RS-R(LL,P)", 6: "RS-C(LL,P)",
#     7: "CS-P", 8: "RS-P",
# }
SCEN_LABELS = {
    1: "1S–1(LL)", 2: "4S–4(LL)", 3: "1S–1(LL,P)",
    4: "4S–4(LL,P)", 5: "1S–P", 6: "4S–P",
}


# ----------------- scenario runner (seeded) -----------------


def normalize_age_group_col(df, col="age_group"):
    """Map age_group codes (p/s/a/o/g) to long labels; leave long labels as-is."""
    if col in df.columns:
        raw = df[col]
        # case-insensitive match on single-letter codes
        mapped = (
            raw.astype(str).str.strip().str.lower()
            .map(AGE_GROUP_MAP)
        )
        # keep original values where no mapping applies (already long labels or NaN)
        df[col] = mapped.where(mapped.notna(), raw)
    return df

# ----------------- load & preprocess -----------------


def _prepare_samples_df(weeks_list):
    if not weeks_list:
        return pd.DataFrame()
    all_samples_df = pd.concat(weeks_list, ignore_index=True)
    if args.date_field in all_samples_df.columns:
        all_samples_df[args.date_field] = pd.to_datetime(all_samples_df[args.date_field], errors="coerce")
    return all_samples_df

# Small helpers for stride-aligned evaluations
def _calendar_week_bounds(week_idx):
    anchor = start_date + timedelta(weeks=week_idx)
    return anchor - timedelta(days=7), anchor - timedelta(days=1)

def _calendar_window_bounds(start_idx, end_idx):
    window_start, _ = _calendar_week_bounds(start_idx)
    _, window_end = _calendar_week_bounds(end_idx)
    return window_start, window_end

def _stride_eval_indices(scfg, n_weeks):
    stride = sampling_stride_weeks(scfg)
    return list(range(stride - 1, n_weeks, stride))

def _sum_hist_window(hist_list, start_idx, end_idx):
    out = pd.Series(dtype=float)
    if not hist_list:
        return out
    upper = min(end_idx, len(hist_list) - 1)
    for j in range(max(0, start_idx), upper + 1):
        out = out.add(hist_list[j], fill_value=0)
    return out

def _filter_df_by_week_window(df, date_col, start_idx, end_idx):
    if df.empty or date_col not in df.columns:
        return df.iloc[0:0].copy()
    window_start, window_end = _calendar_window_bounds(start_idx, end_idx)
    mask = (df[date_col] >= window_start) & (df[date_col] <= window_end)
    return df.loc[mask]



def _cum_kl_vs_stride(hist_list, ref_hist_list, scfg):
    n = min(len(hist_list), len(ref_hist_list))
    xs, ys = [], []
    for end_idx in _stride_eval_indices(scfg, n):
        sample_counts = _sum_hist_window(hist_list, 0, end_idx)
        ref_counts = _sum_hist_window(ref_hist_list, 0, end_idx)
        xs.append(end_idx + 1)
        if sample_counts.sum() == 0 or ref_counts.sum() == 0:
            ys.append(np.nan)
            continue
        ys.append(kl_dist(sample_counts / sample_counts.sum(), ref_counts / ref_counts.sum()))
    return xs, ys

def _window_kl_vs_stride(hist_list, ref_hist_list, scfg, window_weeks=None):
    n = min(len(hist_list), len(ref_hist_list))
    stride = sampling_stride_weeks(scfg)
    win = int(window_weeks or stride)
    xs, ys = [], []
    for end_idx in _stride_eval_indices(scfg, n):
        start_idx = max(0, end_idx - win + 1)
        sample_counts = _sum_hist_window(hist_list, start_idx, end_idx)
        ref_counts = _sum_hist_window(ref_hist_list, start_idx, end_idx)
        xs.append(end_idx + 1)
        if sample_counts.sum() == 0 or ref_counts.sum() == 0:
            ys.append(np.nan)
            continue
        ys.append(kl_dist(sample_counts / sample_counts.sum(), ref_counts / ref_counts.sum()))
    return xs, ys

def _axes_for_algos(n_algo: int, figsize_per_col=(7, 6)):
    """
    Create a 1 x n_algo row of axes, sharing Y.
    Returns: fig, [axes...]
    """
    fig, axes = plt.subplots(1, n_algo, figsize=(figsize_per_col[0]*n_algo, figsize_per_col[1]), sharey=True)
