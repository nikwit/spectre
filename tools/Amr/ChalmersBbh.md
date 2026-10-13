<!--
Distributed under the MIT License.
See LICENSE.txt for details.
-->

# First BBH AMR pilot

Use `chalmers-amr` from `nikwit/spectre` on the cluster. Build
`EvolveGhBinaryBlackHole` with the cluster environment. The first pilot uses
`SpecSpectral`: p refinement and coarsening on a fixed h mesh. There is no
Chalmers refinement request, junk-radiation gate, neighbor target smoothing,
or adaptive filter in this experiment.

Start from an already rendered, working **rectangular** BBH inspiral input and
its real initial data. An equal-mass, nonspinning early inspiral is the first
choice. The generator preserves the gauge, moving maps, horizon controllers,
constraint damping, local time stepping, boundary conditions, and steep
Hypercube filter. All blocks must use tensor-product Legendre bases.

```sh
build/bin/python-spectre tools/Amr/ChalmersBbh.py /path/Inspiral.yaml \
  --output /path/bbh-amr-inputs --checkpoint-hours 1.5
build/bin/EvolveGhBinaryBlackHole \
  --input-file /path/bbh-amr-inputs/p-smoke.yaml --check-options
```

The checkpoint wallclock limit must leave room before the allocation ends.
Use `--input-run-dir` when relative initial-data paths refer to a different
directory than the base YAML. Validate each generated input with the cluster
executable, use the normal cluster launcher, and run each case in a separate
output directory. The generator submits no jobs.

Run in this order:

1. `p-smoke`: 20 time units. Check initialization, controls, excision speeds,
   timestep behavior, and actual point additions and removals.
2. `restart-control` and `restart-checkpoint`: matched 20-unit runs, except
   the latter writes a checkpoint at 10 units after that time's AMR phases.
   Restart separately with the same binary and node/core layout using the
   usual `+restart` launch. Compare fields, meshes, subsequent decisions,
   cooldown timestamps, and target history. Keep output files separate.
3. `fixed`, `fixed-fine`, `spec`, and `spec-tight`: 200 units by default.
   `fixed-fine` adds two points per direction. The adaptive pair differs
   only in spectral and constraint targets, all multiplied by
   `--tight-factor` (default 0.25) for `spec-tight`.

All adaptive cases permit p coarsening and preserve the initial h mesh.
AMR starts at the first scheduled evaluation, without waiting for junk to
leave the wave zone. The default interval is 0.5, with a two-call cooldown
and a minimum elapsed time of 1 after each realized p change. These times
are in evolution units (M for total mass one); select
`--time-between-changes` explicitly for other mass scales. SpEC also uses a
shorter elapsed-time restriction in its frequent-trigger mode, so 1 is a
pilot setting rather than a universal SpEC constant.

The spatial spectral target is a background plus two Gaussian contributions:

```
w_A = max_over_element exp(-|x_grid - center_A|^2 / width_A^2)
log10(T) = (log10(T_background) + w_A log10(T_A) + w_B log10(T_B))
           / (1 + w_A + w_B)
```

Centers come from the domain's object X coordinates and center-of-mass Y/Z
offset. Each width is the absolute value of the other object's X coordinate,
matching the simple SpEC BBH input prescription. This implementation uses
SpEC's grid-coordinate profile convention, without its optional mapping to
inertial distances. The profile is sampled at the first AMR evaluation and
retained across p changes and restarts. With `--spectral-target 1e-4` and
`--near-hole-factor 0.01`, the Gaussian target is 1e-6; a unit-weight isolated
Gaussian at its center blends to 1e-5, not 1e-6. Set the factor to 1 for a
uniform target.

Other defaults: 7–20 points per direction, h policy limits 1–4, and projected
C2/C3 RMS targets 0.002/0.004. These are provisional BBH settings,
not calibrated accuracy requirements. Native nodal RMS targets cannot be
copied numerically
from SpEC's topological norm; see GeneralizedHarmonicExperiments.md. The
stored target shift remains bounded by [-1, 0]. Regional minimum-mode arrays
can be set in the generated criterion if needed; initialize at or above each
floor. Requests beyond mesh limits stop for inspection. Keep this pilot in
early inspiral, before the base input changes slab size near merger.

`AmrConstraints` records RMS/max gauge, two-index and three-index constraints;
`AmrFields` saves double-precision fields every 5 units. AMR statistics and
terminal decisions record work and point-count changes. Existing horizon,
excision, timestep, and waveform outputs are retained. Add
`--mesh-observations` for a full mesh history at every AMR call.
Compare regional constraints near each hole and in the wave zone, horizon
masses/spins, waveform
differences after junk, and cost against both fixed references. A 200-unit
pilot tests integration; it does not establish full-inspiral waveform accuracy.

The local parsing fixture in `build-agents/chalmers-bbh-preparation` has no
real BBH initial data and must not be evolved. No BBH evolution has yet been
performed locally. Older checkpoints predate the expanded AMR history layout;
start new runs with the rebuilt executable.
