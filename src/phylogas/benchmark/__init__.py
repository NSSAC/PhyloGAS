"""Benchmarking: compare inferred phylodynamics against the simulation's truth.

These metrics all require ground truth from the agent-based model, which is
why they live in PhyloGAS rather than BeyondBaseline. BeyondBaseline is the
sampling-strategy engine and must stay runnable by a health department on a
real linelist with no ABM at all.

Not every metric here is sequence-facing. Mugration inference, for example,
works purely on the transmission graph. The dividing line is *ground truth*,
not *sequences*.
"""
