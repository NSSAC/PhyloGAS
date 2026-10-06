# Where chain identity and per-infection metadata should live

Status: **open**. Nothing here is implemented. Written 2026-10-06 so the
decision is made once rather than re-derived.

## The question

`mu_lineage` needs to know which transmission chain each painted sequence
belongs to. Today that comes from joining TwinSampler's all-events table on
`alias_pid`. It works, but it raised a design question worth settling: should
the painter consume TwinSampler's metadata table -- giving one source carrying
the chain id, the seed id, demographics and geography together -- or should
each package compute what it needs?

The underlying concern is duplication. `simulate_linelist.py` already carries
a deliberate second copy of the `strain` id construction, with
`test_strain_id_contract.py` existing only to stop the two drifting:

> The canonical definition is PhyloGAS's `phylogas/ids.py`; this copy is
> deliberate so TwinSampler runs without PhyloGAS installed, and
> `test_strain_id_contract.py` pins both to the same fixture table.

Adding a second implementation of *chain identity* would repeat that pattern.
Cross-checking two implementations against each other is worth doing once as
an audit; it is not worth maintaining as permanent machinery.

## Two constraints that narrow the options

**1. The painter must traverse every infection, whatever it writes.**
From `_OutputSet`'s docstring in `genetic_painter.py`:

> The painter always walks the whole transmission tree -- a child's genome is
> derived from its parent's, so nothing can be skipped during computation.
> These objects decide only what reaches disk.

So TwinSampler's table can filter or decorate the painter's *output*. It can
never define the painter's computational universe. Any proposal phrased as
"the painter consumes the all-events table instead of the graph" has to mean
"consumes it for output decisions", because the graph walk is not optional.

**2. TwinSampler must run without PhyloGAS.** It imports nothing from
`phylogas` and does not declare it as a dependency -- by design, so a health
department can run it on a real line list. Its `component_id` therefore
cannot be removed.

Constraint 2 decides the ownership question. If the **painter** owns chain
identity, TwinSampler still computes components for its standalone use, so
there are necessarily two implementations. If **TwinSampler** owns it, there
can be exactly one. The option that actually eliminates the duplication is
the one where TwinSampler owns chain identity and the painter consumes it.

That is the opposite of what the shape of the pipeline suggests -- the painter
is the component that literally traces the chains -- but standalone operation
is a hard requirement and duplication-elimination is the goal.

## There is already a precedent

`--linelist-filter` has the painter reading a TwinSampler-shaped file and
restricting an output set to the identifiers in it. `_build_output_sets`
parses `label=path`, `_read_filter_keys` pulls ids or pids out of whichever
column it finds, and one pass emits the full tree plus any number of subsets.

So "the painter consumes TwinSampler output" is not a new architecture. The
extension being considered is from *filtering* on that file to *decorating*
from it -- carrying selected columns through onto the records the painter
writes. That is a smaller change than the framing suggested.

## What the join costs today

Measured on the 36-week Virginia Delta run:

| | |
|---|---|
| painted records | 5,350,395 |
| allevents infections (distinct `alias_pid`, 1:1 with `strain`) | 5,295,971 |
| painted infections with no allevents row | 54,424 (1.02%) |
| allevents rows before de-duplication | 8,902,620 (1.68 per infection) |

Three costs:

- **Memory.** The `alias_pid -> component_id` dict is roughly a gigabyte, on
  every benchmark run.
- **Coverage.** The 1.02% with no allevents row drop out of `mu_lineage`,
  counted as `skipped_no_chain`. They are infections whose clinical state
  fell outside the window or did not classify as a severity.
- **A dependency on another package's graph construction.** This has already
  gone wrong once: before `--target_variant`, the tick filter severed the
  graph into 491,535 components against 3,641 importations. The clock
  measurement should not be the thing that discovers such a bug.

None of these are blocking. The join works and `mu_lineage` is in
`clock_estimates.csv` today.

## Open questions

**1. What does the painter write, and keyed how?**
The universe is settled by constraint 1 -- every `E2` exposure is traversed.
But the output could stay as it is (one record per traversed infection, its
own metadata), or carry columns joined from TwinSampler's table, or be
TwinSampler's table with the painter's columns merged in. The third makes the
painted set's row count depend on TwinSampler's filters, which changes what
`all_infections` means.

**2. Does the painter stay runnable with only a graph?**
Useful for testing, and it is how the painter runs today. If yes, it needs a
no-TwinSampler path -- which is fine as an explicitly degraded test mode, but
becomes a second implementation if it has to produce correct chain ids.

**3. If chain identity is joined in, what happens to the 1% without a row?**
Written with a null chain, or dropped from the output? Dropping changes the
painted set; writing nulls keeps it and leaves the gap visible downstream,
which is the current behaviour of the benchmark join.

## If the painter were to own it instead

Recorded because it was considered and is cheap, even though it does not
eliminate the duplication.

The painter maintains `current_sequences[infection_id]`, assigning each child
from its parent: `current_sequences[infection_id] = mutate(current_sequences[contact_infection_id])`.
A parallel dict propagated identically -- `founder_of[child] = founder_of[parent]`,
each index case its own founder -- gives every record its founding infection
with no join, full coverage, and no cross-package dependency. Roughly ten
lines plus a metadata column.

Two details if this is ever done:

- Propagate the founder's **`alias_pid`**, not `real_strain`. There are 3,322
  seeds for 3,306 importations, so a seed sequence can be reused; grouping by
  seed would merge chains imported on different dates into one fixed-effects
  group and bias `mu_lineage`. The founder's infection id is unique per chain
  by construction. Carrying both gives the seed-trend rows exact import dates
  and seed divergences.
- It costs a repaint. Sequences would not change -- the dict never touches the
  RNG -- so it is metadata-only, but Snakemake cannot know that.

## Recommendation

Leave it as it is for now; the measurement is unblocked. When it is taken up,
the version worth building is TwinSampler owning chain identity with the
painter decorating from its table, extending `--linelist-filter` rather than
reordering the pipeline. Answer question 2 first: whether the graph-only path
has to produce correct chain ids or may be an explicitly degraded test mode.
That single answer determines whether one implementation is achievable at all.

## Background

- `docs/clock_modes.md` sections 1a and 1b: what `mu_lineage` and the
  seed-trend quantities are, and why the pooled rate is the one comparable to
  TreeTime.
- `src/phylogas/benchmark/lineage.py`: the join as currently implemented.
