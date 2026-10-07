# Variant overlays

An overlay is an importation schedule whose `variant` labels are laid over a
simulation, so the prevalence benchmark can test whether a sampling strategy
recovers a chosen prevalence profile. `phylogas assign-variants` matches each
simulated transmission chain's first tick to the schedule's ticks and copies
that row's `variant` into `variant_benchmark`.

Format: `tick,date,variant,clusters,sample_count` (one row per day and
variant; `clusters` importations that day).

Select one with `benchmark.variant_schedule` in the config. Without it, the
default is the schedule `prep-seeds` writes from the real importations, with
each cluster's sublineage as its label.

| file | contents |
|---|---|
| `Virginia_importation_schedule.csv` | The original overlay, from TwinSampler's `Data/`: 2,568 importation clusters of XBB/JN.1-era variants (JN.1, XBB.1.5, XBB.1.16, EG.5, ...), 2022-11-29 to 2024-10-18, on its own ticks 0-689. Manufactured labels: unrelated to the Delta wave the simulation models, by design. |
