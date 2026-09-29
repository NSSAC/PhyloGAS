# Review guide

What changed across the two working sessions, ordered by how carefully you
should read it. Nothing has been committed.

---

## Already committed by you — not under review here

Commits `d6b1d75` ("speed up and fixing bug where the importation sequence is
being overwritten...") through `b8d014b` ("location based adjustment for bash
scripts") are on `origin/main` and contain **all** of the following:

- the re-importation / single-pass fix in `genetic_painter.py`
- the `rate_limited.py` rewrite
- the `simple.py` / `poor.py` dtype fixes
- all painter performance work (`S1` storage, eviction, CSV dtypes, date LUT,
  compression flags)
- the SLURM script resource, `set -u` ordering, and path-anchoring changes

Verified by diffing every moved file against `HEAD`: they are byte-identical.
Nothing from those sessions is sitting uncommitted.

---

## Tier 1 — the only content change to existing code

### `src/phylogas/painter/genetic_painter.py` — 15 lines, 2 hunks

This is the entire delta to previously-committed code. Everything else in this
round is a pure rename or a new file.

```diff
-from mutational_models import registry as model_registry
+try:
+    from .mutational_models import registry as model_registry
+except ImportError:
+    from mutational_models import registry as model_registry
```
Needed because the module is now inside a package. The fallback keeps the
`run.03.vadelta.*` scripts working when they invoke the file directly —
verified by running `run.03.vadelta.a` in place after the move.

```diff
     output_file_prefix = args.output_prefix
+    _out_parent = os.path.dirname(os.path.abspath(output_file_prefix))
+    if _out_parent:
+        os.makedirs(_out_parent, exist_ok=True)
```
Bug found while testing the Snakemake pipeline: the painter died after loading
the network and seed FASTA because the output directory did not exist. Minutes
of wasted work on a real run.

Inspect with:
```bash
git diff HEAD -M --find-renames=50% -- src/phylogas/painter/genetic_painter.py
```

---

## Tier 2 — context for behaviour you already committed

These are already on `origin/main`; listed so the rationale is recorded
somewhere durable, not because they need re-review.

### Compression defaults — `--compression_level` defaults to **1**, not 6

This is the single biggest runtime lever and it changes your output file sizes.

| preset | throughput | ratio |
|---|---|---|
| 6 (old lzma default) | ~243 rec/s | ~1500x |
| 1 (new default) | ~9,000 rec/s | ~320x |

Level 6 cost ~6 of your 8.68 hours. Level 1 files are a few times larger.
**If archive size matters more than wall time, set `compression_level: 6` and
`compression_threads: 8`** — 13.2s vs 70.4s per 315 MB, keeping most of the ratio.

Your call; I picked speed as the default.

### 5. Memory and IO changes in the painter

- Genome storage `<U1` -> `S1`: 117 KB -> 29 KB per sequence, and
  `tobytes()` replaces `''.join(tolist())` (486 us -> 2.7 us)
- Superseded genomes evicted on reinfection (proved safe: 0 references to an
  evicted genome across 29,009 edges)
- EpiHiper CSV read: 4 of 5 columns, `int32` tick, categorical `exit_state`
  (~7.5 GB -> ~1.9 GB), vectorised `isin` instead of a per-row lambda (26x)
- tick->date lookup table instead of 8M `Timestamp` objects (62s -> 0.25s)

Net: **8.68 h -> ~20 min**, measured peak RSS ~17 GB (was requested at 1 TB).

### 6. `src/phylogas/painter/test/*.sh` — SLURM scripts

- Resources: `bii-largemem`/1TB/24h/30 cores -> `standard`/64G/4h/8 cores
- `set -eo pipefail` then `set -u` **after** conda activation (the
  `GDAL_DATA: unbound variable` failure you hit)
- Scripts anchor to their own location via `BASH_SOURCE`; the wrapper uses
  `SLURM_SUBMIT_DIR` because sbatch copies the batch script to a spool dir
- `code=` path corrected — I had it wrong twice; now verified from three cwds

**Re-check the resource numbers against a real `seff` before trusting them at
full scale.** My 17 GB projection is extrapolated from a 1.79M-record run.

---

## Tier 3 — new scaffolding (skim for design agreement)

| File | What to check |
|---|---|
| `pyproject.toml` | Dependency pins; the `phylogas` entry point |
| `src/phylogas/cli.py` | Verb names (`train`/`paint`/`prep-seeds`/`subset-fasta`/`run`/`validate-config`/`fetch-data`) |
| `src/phylogas/config.py` | `{placeholder}` expansion semantics |
| `config.template.yaml` | **Key names** — these become your public API |
| `Snakefile` | Rule granularity, `resources:` defaults |
| `command_map.md` | The old->new map you asked for |
| `README.md` | New "Implementation Status" table marks what is *not* built |

**Design decision worth your veto:** I used argparse, not Typer/Click as
vision.txt suggests. Neither was installable in my test environment, and
argparse is stdlib so `phylogas --help` works before any dependency resolves.
Same UX. Easy to swap later.

---

## Tier 4 — documents to read, not code

| File | Purpose |
|---|---|
| `../BeyondBaseline/UNCOMMITTED_CHANGES.md` | **What your lost-context changes were for.** Answer: packaging into 6 CLI commands. Complete and working. |
| `docs/repo_split_audit.md` | TwinSampler/BeyondBaseline split — what is safe to delete |
| `docs/salvage_audit.md` | What was rescued from `synthetic_biosurveillance` |

---

## File moves (history preserved via `git mv`)

```
src/genetic_painter/      -> src/phylogas/painter/
src/seq_prep/             -> src/phylogas/seqprep/
src/filter_variant_time.py-> src/phylogas/filter_variant_time.py
(new)                        src/phylogas/popprep/   <- rescued from synthetic_biosurveillance
```

`git log --follow` still works. Legacy scripts verified still runnable after
the move.

---

## What I verified, and what I did not

**Verified by execution:**
- Full `train` -> `paint` Snakemake run: 29,010 sequences, matching metadata,
  correct no-op on re-run
- Parent-child divergence after the fix: mean 0.166, max 3, **0 pairs >8 diffs**
- Clean-venv install; all 7 subcommands respond
- All three mutation models; metadata augmentation; xz/bgzf/uncompressed paths
- BeyondBaseline's 6 entry points install and run
- Legacy `run.03.vadelta.*` scripts still work post-restructure

**NOT verified — needs your attention:**
- Anything at full 5.35M-record scale (largest local test: 3.3M)
- The `simulate_linelist` and `sample_scenarios` Snakemake rules — written from
  your `vision.txt` command lines, never executed (needs TwinSampler and
  BeyondBaseline installed)
- `nextstrain_build` — needs a local ncov checkout
- Resource sizing on the real cluster
- `phylogas fetch-data` is a **stub**; `phylogas benchmark` **does not exist**

---

## Suggested commit sequence

Painter correctness and performance are **already committed** (`d6b1d75`..`b8d014b`).
What remains is the restructure:

1. **Move into a package** — the 21 `git mv` renames + the 15-line
   `genetic_painter.py` change. Stage with `git add -A src/` so git records
   them as renames; verify with
   `git diff HEAD -M --stat` (every moved file should show `0`).
2. **Packaging** — `pyproject.toml`, `.gitignore`, the `__init__.py` files.
3. **CLI + config** — `src/phylogas/cli.py`, `src/phylogas/config.py`,
   `config.template.yaml`.
4. **Pipeline** — `Snakefile`.
5. **Rescued code** — `src/phylogas/popprep/` (from `synthetic_biosurveillance`).
6. **Docs** — `command_map.md`, `docs/`, README update.

1 and 2 belong together; 3-6 are independent.

### Verifying the renames before you commit

```bash
git add -A src/
git diff HEAD -M --find-renames=50% --stat -- src/
```

Every relocated file should report `0` changed lines. The only non-zero entry
should be `genetic_painter.py` at `15 +-`. If you see large insert/delete pairs
instead of `{old => new}` arrows, the rename detection failed and the move
should be redone.
