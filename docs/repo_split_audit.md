# TwinSampler / BeyondBaseline split: what is safe to delete

Audit performed 2026-09-28. Every shared file was diffed, dated against its
last commit in each repo, and checked for import coupling and symbol loss.

**Headline: the split is clean and safe to finish. Each repo is authoritative
for its own domain, there is zero code coupling between the two halves, and no
functionality is lost by deleting the non-owned half from each side — with one
exception that must be handled first (see "Blocking issue").**

---

## The ownership rule

| Directory | Owner | Rationale |
|---|---|---|
| `scripts/linelist_generation/` | **TwinSampler** | TS commits are ~2.5 months newer |
| `scripts/scenarios_simulation/` | **BeyondBaseline** | BB commits are ~2.5 months newer |
| `utils/`, `Data/`, `LineList/`, `Algorithm_results/`, `demo/` | see below | identical in both |

Divergence runs in **opposite directions**, which is the strongest possible
evidence that each repo is the live copy of its own half:

| File | Diff lines | TwinSampler last commit | BeyondBaseline last commit |
|---|---|---|---|
| `linelist_generation/simulate_linelist.py` | 314 | **2026-05-29** | 2026-03-18 |
| `linelist_generation/label_components.py` | 242 | **2026-05-29** | 2026-03-12 |
| `linelist_generation/generate_linelist.sh` | 8 | **2026-05-29** | 2026-03-17 |
| `scenarios_simulation/run_all_scenarios.py` | 1526 | 2026-03-30 | **2026-05-30** |
| `scenarios_simulation/sampling_algorithms.py` | 58 | 2025-12-15 | **2026-06-11** |
| `scenarios_simulation/scenarios_config.py` | 14 | 2026-03-23 | **2026-05-31** |
| `scenarios_simulation/README.md` | 87 | 2026-03-25 | **2026-04-29** |

Everything else shared (33 files: notebooks, `LICENSE`, `Data/`, `utils/sampling.py`,
`ascertainment_module.py`, `place_resolver.py`, `rucc_utils.py`, ...) is **byte-identical**.

---

## Confirmed: no coupling in either direction

```
scenarios_simulation -> linelist_generation :  no imports
linelist_generation  -> scenarios_simulation:  no imports
```

`BeyondBaseline/pyproject.toml` already excludes `linelist_generation` from the
package (`include = ["scenarios_simulation*"]`), so the code it ships is
already correct — the extra directory is dead weight on disk, not a dependency.

---

## Confirmed: sunsetting TwinSampler's `scenarios_simulation/` loses nothing

- **File level:** BeyondBaseline is a strict superset. Every file in
  `TwinSampler/scripts/scenarios_simulation/` exists in BeyondBaseline.
- **Symbol level:** only two functions exist in TS and not BB —
  `_cum_kl_vs` and `_roll_kl_vs` in `run_all_scenarios.py`.

  These are **not lost**. BeyondBaseline replaced them with stride-aware
  successors `_cum_kl_vs_stride` / `_window_kl_vs_stride`, which support the
  multi-week stride evaluation (`4S-4`) that the older pair could not express.

- BeyondBaseline additionally has 13 functions TwinSampler lacks
  (`_calendar_week_bounds`, `_stride_eval_indices`, `split_samples_by_calendar_week`,
  `target_dist_at_week`, ...) plus whole modules TS never had:
  `mugration_station.py`, `aggregate_results.py`, `lasso_test/`.

**Verdict: `TwinSampler/scripts/scenarios_simulation/` can be deleted outright.**

---

## Blocking issue: do NOT delete TwinSampler's linelist code

Four files exist **only** in TwinSampler:

| File | Status |
|---|---|
| `scripts/linelist_generation/demographics_module.py` | **Required** — imported by TS's `simulate_linelist.py` |
| `scripts/linelist_generation/get_pairwise_seq_id.py` | Standalone utility, unreferenced |
| `Data/get_county_table.py` | Fetches the census county table |
| `Data/county_fips.csv` | Its output; consumed by `place_resolver.py` |

`demographics_module.py` provides `DemographicsLoader`, used at three points in
TwinSampler's `simulate_linelist.py`. BeyondBaseline's copy of
`simulate_linelist.py` is an older snapshot that predates that refactor, so it
does not import it — which is exactly why BB's copy is stale and must be the one
deleted.

Verified: TwinSampler's `simulate_linelist.py --help` runs; it is the live copy.

---

## Recommended actions

### In BeyondBaseline — delete the linelist half

```bash
cd ../BeyondBaseline
git rm -r scripts/linelist_generation
git commit -m "linelist generation is owned by TwinSampler; remove stale copy"
```

Stale by ~2.5 months and missing `demographics_module.py`, so anyone running it
gets an older ascertainment model than TwinSampler produces.

### In TwinSampler — delete the scenarios half

```bash
cd ../TwinSampler
git rm -r scripts/scenarios_simulation
git commit -m "sampling scenarios are owned by BeyondBaseline; remove superseded copy"
```

BeyondBaseline is a strict superset; confirmed above at both file and symbol level.

### Both repos — decide on the duplicated commons

These are byte-identical in both and were simply inherited at fork time:

| Path | Suggestion |
|---|---|
| `utils/sampling.py` | **Orphaned in both** — nothing imports it. Delete from both, or move into whichever package actually needs it. |
| `Data/Virginia_importation_schedule.csv` | Input to TwinSampler. Keep in TS; drop from BB. |
| `Data/ReadMe.md` | Follows the above. |
| `LineList/*.ipynb` (6 notebooks) | Linelist analysis -> TwinSampler. |
| `Algorithm_results/*.ipynb` (5 notebooks) | Sampling-algorithm results -> BeyondBaseline. |
| `demo/demo.ipynb` | Check which half it demonstrates; keep one. |
| `LICENSE`, `scripts/ReadME.md` | Keep in both (correct). |

### Also worth doing in BeyondBaseline

From `UNCOMMITTED_CHANGES.md` in that repo:

- Commit the packaging work. **PhyloGAS's `environment.yml` installs
  BeyondBaseline from GitHub and cannot work without `pyproject.toml`.**
- Add `.gitignore` for `replicate_results*/`, `*.tar.gz`,
  `sim_replicates_*.out.txt`, `__pycache__/`, `.venv/`, `.DS_Store`, `z`,
  `*.egg-info/`. **Do not commit `lasso_greedy50.tar.gz` or
  `replicate_results.tar.gz`** — they are permanent repo bloat.
- Decide where `development_dialog.json` (7 MB) lives. Valuable provenance, but
  it is a transcript; consider `docs/` outside git, or Git LFS.

---

## Suggested order

1. Commit BeyondBaseline's packaging work (unblocks PhyloGAS install).
2. Add `.gitignore` to BeyondBaseline; unstage the artifacts.
3. `git rm -r scripts/linelist_generation` in BeyondBaseline.
4. `git rm -r scripts/scenarios_simulation` in TwinSampler.
5. Resolve the duplicated commons table above.
6. Delete the stale genetic painter from `synthetic_biosurveillance`
   (see `salvage_audit.md`).

Steps 3 and 4 are independent and safe in either order.
