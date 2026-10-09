"""Tree coverage (Mean Reciprocal Distance) and the truth runner around it.

`calculate_coverage_score` is BeyondBaseline's function, ported verbatim. The
runner computes the same quantity in two halves -- one BFS per sampled set,
then a cheap sum per target set -- so a week's tree-size thresholds and age
groups share one search. These tests pin the split version to the original,
and the runner's windowing to a direct, per-call reimplementation.

    pytest tests/test_tree_coverage.py -v
"""

from __future__ import annotations

import random

import numpy as np
import pandas as pd
import pytest

T = pytest.importorskip("phylogas.benchmark.truth_metrics")
R = pytest.importorskip("phylogas.benchmark.truth_runner")


def _random_graph(rng, n):
    nodes = [str(i) for i in range(n)]
    adj = {u: [] for u in nodes if rng.random() < 0.85}   # some nodes off-graph
    for _ in range(rng.randint(0, 2 * n)):
        a, b = rng.sample(nodes, 2)
        if a in adj and b in adj:
            adj[a].append(b)
            adj[b].append(a)
    return nodes, adj


def test_split_bfs_matches_original():
    rng = random.Random(7)
    for _ in range(300):
        nodes, adj = _random_graph(rng, rng.randint(2, 80))
        sampled = set(rng.sample(nodes, rng.randint(0, len(nodes))))
        dist = T.distances_from(sampled, adj)
        for _ in range(3):
            target = set(rng.sample(nodes, rng.randint(0, len(nodes))))
            assert T.coverage_from_distances(target, sampled, dist) == \
                T.calculate_coverage_score(target, sampled, adj)


def test_treegraph_matches_original():
    """The compiled search on the sparse graph gives the original's values."""
    pytest.importorskip("scipy.sparse.csgraph")
    rng = random.Random(3)
    for _ in range(200):
        n = rng.randint(2, 80)
        nodes = [str(i) for i in range(n)]
        pids = nodes[:]
        contacts = [rng.choice(nodes + ["-1", "nan"]) for _ in nodes]   # forest + roots
        frame = pd.DataFrame({"alias_pid": pids, "alias_contact": contacts})
        adj = T.build_undirected_adj(frame.copy())
        g = T.TreeGraph(frame["alias_pid"], frame["alias_contact"])
        sizes = T.precompute_component_sizes(adj, set(pids))
        assert all(g.size[g.codes([u])[0]] == sizes[u] for u in pids)
        for _ in range(3):
            sampled = set(rng.sample(nodes, rng.randint(0, n))) | ({"ghost"} if rng.random() < .2 else set())
            target = sorted(set(rng.sample(nodes, rng.randint(0, n))) | ({"ghost"} if rng.random() < .2 else set()))
            dist = g.distances(sampled)
            got = T.coverage_from_codes(g.codes(target), np.array(target, dtype=object), sampled, dist)
            assert got == pytest.approx(T.calculate_coverage_score(set(target), sampled, adj), abs=1e-12)


@pytest.mark.parametrize("target,sampled,expected", [
    (set(), {"a"}, 0.0),             # no cases
    ({"a"}, set(), 0.0),             # no samples
    ({"z"}, {"z"}, 1.0),             # sampled node absent from the graph
    ({"z", "a"}, {"z"}, 0.5),        # ...and an unreachable case
    ({"b"}, {"a"}, 0.5),             # one hop: 1/(1+1)
])
def test_edge_cases(target, sampled, expected):
    adj = {"a": ["b"], "b": ["a"]}
    assert T.calculate_coverage_score(target, sampled, adj) == expected
    assert T.coverage_from_distances(target, sampled, T.distances_from(sampled, adj)) == expected


def test_recipe_stride():
    assert R._recipe_stride("4S__surs", 9) == 4
    assert R._recipe_stride("1S-P__lasso_greedy", 9) == 1
    assert R._recipe_stride("custom", 9) == 9


def test_equity_name_matches_beyondbaseline():
    assert R._equity_name("Older adult (50-64)") == "equity_Older_adult_50-64"
    assert R._equity_name("Senior (65+)") == "equity_Senior_65+"


START = "2021-04-05"     # a Monday, like BeyondBaseline's default
MIN_POOL = 5
RACES = {"W": "White", "B": "Black", "L": "Latino"}


def _toy_epidemic(tmp_path):
    """A forest of transmission chains, a line list, two sample sets and a
    population file shaped like BeyondBaseline's input."""
    rng = random.Random(11)
    t0_date = pd.Timestamp("2021-03-29")
    people = []
    for sim_pid in range(3000):
        people.append(dict(sim_pid=str(sim_pid), age_group=rng.choice("psaog"),
                           smh_race=rng.choice(list(RACES)), county=rng.choice(["Albemarle", "Fairfax"]),
                           gender=rng.choice([1, 2])))
    pop = pd.DataFrame(people)
    pop.to_csv(tmp_path / "population.csv", index=False)
    by_pid = pop.set_index("sim_pid")

    rows, pid = [], 0
    for comp in range(160):            # dense enough that every week clears MIN_POOL
        size = rng.choice([1, 2, 5, 20, 150])
        t0 = rng.randint(0, 100)
        ids = []
        for k in range(size):
            parent = rng.choice(ids) if ids else "-1"
            day = t0 + k // 6 + rng.randint(0, 3)
            person = str(rng.randrange(3000))
            d = by_pid.loc[person]
            rows.append(dict(alias_pid=str(pid), alias_contact=parent, sim_tick=day, sim_pid=person,
                             date=(t0_date + pd.Timedelta(days=day)).date().isoformat(),
                             age_group=d.age_group, smh_race=RACES[d.smh_race], county=d.county,
                             sex={1: "male", 2: "female"}[d.gender], component_id=comp,
                             variant_benchmark=("background" if comp % 7 == 0 else f"v{comp % 3}")))
            ids.append(str(pid))
            pid += 1
    inf = pd.DataFrame(rows)
    inf.to_csv(tmp_path / "allevents.csv", index=False)
    ll = inf.sample(frac=0.4, random_state=3)
    ll.to_csv(tmp_path / "linelist.csv", index=False)
    ll.sample(frac=0.2, random_state=1).to_csv(tmp_path / "4S__surs_samples.csv", index=False)
    ll.sample(frac=0.4, random_state=2).to_csv(tmp_path / "1S__dense_samples.csv", index=False)
    return inf


def _score(tmp_path, **kw):
    return R.score_samples([str(tmp_path / "*_samples.csv")], str(tmp_path / "allevents.csv"),
                           linelist=str(tmp_path / "linelist.csv"),
                           start_date=START, min_pool=MIN_POOL, **kw)


def _calendar(dates):
    week0 = pd.Timestamp(START) - pd.Timedelta(days=7)
    return (pd.to_datetime(dates) - week0).dt.days // 7


def test_runner_matches_per_call_reference(tmp_path):
    """Every metric against a direct reimplementation calling the original
    `calculate_coverage_score` per target set. The line list fixes the number
    of weeks; the coverage metrics run over every infection on the true
    transmission graph."""
    _toy_epidemic(tmp_path)
    df = _score(tmp_path)

    L = pd.read_csv(tmp_path / "linelist.csv", dtype={"alias_pid": str})
    L["_wk"] = _calendar(L["date"])
    n = 0
    while (L["_wk"] == n).sum() >= MIN_POOL:
        n += 1
    assert n > 4

    A = pd.read_csv(tmp_path / "allevents.csv", dtype={"alias_pid": str, "alias_contact": str})
    A = T.normalize_age_group_col(A, "age_group")
    A["_wk"] = _calendar(A["date"])
    adj = T.build_undirected_adj(A.copy(), "alias_pid", "alias_contact")
    sizes = T.precompute_component_sizes(adj, set(A["alias_pid"]))

    for recipe in ("4S__surs", "1S__dense"):
        s = pd.read_csv(tmp_path / f"{recipe}_samples.csv", dtype={"alias_pid": str})
        s["_wk"] = _calendar(s["date"])
        s = s[s["_wk"].between(0, n - 1)]
        stride = int(recipe[0])
        ref = {}
        for end in range(stride - 1, n, stride):
            S = set(s.loc[s["_wk"] <= end, "alias_pid"])
            P = A.loc[A["_wk"] <= end]
            for t in R.COVERAGE_SIZE_THRESHOLDS:
                ref.setdefault(f"coverage_size_{t}", []).append(T.calculate_coverage_score(
                    {u for u in set(P["alias_pid"]) if sizes.get(u, 1) > t}, S, adj))
            for ag in sorted(A["age_group"].dropna().unique()):
                ref.setdefault(R._equity_name(ag), []).append(T.calculate_coverage_score(
                    set(P.loc[P["age_group"] == ag, "alias_pid"]), S, adj))
            lo = max(0, end - R.ROLLING_TREE_WEEKS + 1)
            ref.setdefault("8_week_rolling_tree_coverage", []).append(T.calculate_coverage_score(
                set(A.loc[A["_wk"].between(lo, end), "alias_pid"]),
                set(s.loc[s["_wk"].between(lo, end), "alias_pid"]), adj))
            w = max(0, end - stride + 1)
            p_hat = s.loc[s["_wk"].between(w, end), "variant_benchmark"].value_counts().drop("background", errors="ignore")
            p_true = A.loc[A["_wk"].between(w, end), "variant_benchmark"].value_counts().drop("background", errors="ignore")
            p_hat, p_true = p_hat / p_hat.sum(), p_true / p_true.sum()
            idx = sorted(set(p_hat.index) | set(p_true.index))
            ref.setdefault("stride_variant_prevalence_error", []).append(float(
                (p_hat.reindex(idx, fill_value=0) - p_true.reindex(idx, fill_value=0)).abs().sum()))
            den = A.loc[A["_wk"].between(w, end), "component_id"].nunique()
            num = s.loc[s["_wk"].between(w, end), "component_id"].nunique()
            ref.setdefault("stride_component_coverage", []).append(num / den if den else np.nan)

        xs = list(range(stride, n + 1, stride))
        got = df[df["algorithm"] == recipe].set_index("eval_type")
        for name, ys in ref.items():
            assert got.loc[name, "weeks"] == len(xs), name
            assert got.loc[name, "auc"] == pytest.approx(R._series_auc(xs, ys), abs=1e-12, nan_ok=True), name


def test_sampling_start_auto_is_first_allevents_week(tmp_path):
    inf = _toy_epidemic(tmp_path)
    first = pd.to_datetime(inf["date"]).min().normalize()
    L = pd.read_csv(tmp_path / "linelist.csv")
    wk0 = ((pd.to_datetime(L["date"]) >= first) & (pd.to_datetime(L["date"]) < first + pd.Timedelta(days=7))).sum()
    start = R.resolve_sampling_start(str(tmp_path / "allevents.csv"), str(tmp_path / "linelist.csv"),
                                     "auto", min_pool=int(wk0))
    assert start == first + pd.Timedelta(days=7)          # week 0 = first all-events week


def test_sampling_start_auto_skips_thin_lead_in(tmp_path):
    """BeyondBaseline stops at the first week below min_pool, so 'auto' must
    not hand it a thin first week."""
    _toy_epidemic(tmp_path)
    L = pd.read_csv(tmp_path / "linelist.csv")
    first = pd.to_datetime(pd.read_csv(tmp_path / "allevents.csv")["date"]).min().normalize()
    counts = ((pd.to_datetime(L["date"]) - first).dt.days // 7).value_counts()
    pool = int(counts.get(0, 0)) + 1                       # week 0 is one short
    start = R.resolve_sampling_start(str(tmp_path / "allevents.csv"), str(tmp_path / "linelist.csv"),
                                     "auto", min_pool=pool)
    k = int((start - first).days // 7) - 1                 # weeks skipped
    assert k >= 1 and counts.get(k, 0) >= pool and all(counts.get(w, 0) < pool for w in range(k))


def test_sampling_start_explicit_is_kept(tmp_path):
    _toy_epidemic(tmp_path)
    start = R.resolve_sampling_start(str(tmp_path / "allevents.csv"), str(tmp_path / "linelist.csv"),
                                     START, min_pool=MIN_POOL)
    assert start == pd.Timestamp(START)


def test_kl_vs_infections_matches_original(tmp_path):
    """cumulative_infections / stride_window_infections against BeyondBaseline's
    own path-based pipeline: its loader, the verbatim build_weekly_infections
    port, its calendar split, and its kl_dist."""
    rab = pytest.importorskip("scenarios_simulation.run_all_scenarios")
    from scenarios_simulation.sampling_algorithms import kl_dist
    _toy_epidemic(tmp_path)
    df = _score(tmp_path, population=str(tmp_path / "population.csv"))

    start = pd.Timestamp(START)
    feats = rab._normalize_stratifiers(list(R.DEFAULT_STRATIFIERS))
    line_df, pop_df, _, weekly_ll = rab.load_linelist_and_population(
        str(tmp_path / "linelist.csv"), str(tmp_path / "population.csv"), "date", start, MIN_POOL, feats)
    n = len(weekly_ll)
    weekly_inf, _ = T.build_weekly_infections(str(tmp_path / "allevents.csv"), pop_df, start, n)

    for recipe in ("4S__surs", "1S__dense"):
        s = pd.read_csv(tmp_path / f"{recipe}_samples.csv", dtype={"alias_pid": str})
        s = s[["alias_pid"]].merge(line_df[["alias_pid", "group", "date"]], on="alias_pid")
        weeks = rab.split_samples_by_calendar_week(s, "date", start, n)
        hist = [weeks[w]["group"].value_counts() if w in weeks else pd.Series(dtype=float)
                for w in range(n)]
        stride = int(recipe[0])
        for name, win in (("cumulative_infections", None), ("stride_window_infections", stride)):
            xs, ys = [], []
            for end in range(stride - 1, n, stride):
                lo = 0 if win is None else max(0, end - win + 1)
                a = pd.concat(hist[lo:end + 1]).groupby(level=0).sum()
                b = pd.concat(weekly_inf[lo:end + 1]).groupby(level=0).sum()
                xs.append(end + 1)
                ys.append(kl_dist(a / a.sum(), b / b.sum()) if a.sum() and b.sum() else np.nan)
            got = df[(df["algorithm"] == recipe) & (df["eval_type"] == name)]["auc"].item()
            assert got == pytest.approx(R._series_auc(xs, ys), abs=1e-12, nan_ok=True), (recipe, name)


def test_tree_size_uses_true_trees(tmp_path):
    """coverage_size_100 must see the 150-infection trees. Built from the line
    list, those trees break into fragments at every unreported person, and the
    metric is 0 for want of any tree that large."""
    _toy_epidemic(tmp_path)
    df = _score(tmp_path)
    assert (df.loc[df["eval_type"] == "coverage_size_100", "auc"] > 0).all()


def test_coverage_ranks_higher_first(tmp_path):
    _toy_epidemic(tmp_path)
    df = _score(tmp_path)
    for et in ("coverage_size_0", "8_week_rolling_tree_coverage", "stride_component_coverage"):
        sub = df[df["eval_type"] == et]
        best = sub.loc[sub["rank_overall"] == 1, "auc"].max()
        assert best == sub["auc"].max(), et
