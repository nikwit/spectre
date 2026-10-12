<!--
Distributed under the MIT License.
See LICENSE.txt for details.
-->

# First BBH AMR pilot

Fetch `chalmers-amr` from `nikwit/spectre` on the cluster. The commits include
the AMR implementation, tests and experiment tools. The BBH executable now
registers `ChalmersWithSpecCoarsening`, initializes its history, and projects
that history on mesh changes. Local build directories and simulation outputs
are not versioned; rebuild with the cluster environment.

Start from an already rendered, working **rectangular** BBH inspiral input and
its real initial data. An equal-mass, nonspinning early inspiral is the first
choice. Keep the usual damped-harmonic gauge, moving maps, horizon controllers,
constraint damping and steep filter. The generator preserves these settings
and local time stepping. Cylindrical and spherical-harmonic grids are not
supported by this criterion yet.

On the cluster, build `EvolveGhBinaryBlackHole` in the normal cluster build.
Generate inputs with that build's Python wrapper (replace the example paths):

```sh
build/bin/python-spectre tools/Amr/ChalmersBbh.py /path/Inspiral.yaml \
  --output /path/bbh-amr-inputs --checkpoint-hours 1.5
build/bin/EvolveGhBinaryBlackHole \
  --input-file /path/bbh-amr-inputs/p-smoke.yaml --check-options
```

The checkpoint time must leave room before the allocation ends (1.5 hours is
only an example for a two-hour allocation). Use `--input-run-dir` if relative
input-data paths refer to a different directory than the base YAML. Validate
all generated inputs with the cluster executable and use the cluster's usual
launcher. Run each case in a separate output directory. There are no submission
scripts or automatic jobs.

Run in this order:

1. `p-smoke`: 20 time units, p refinement and coarsening, fixed h mesh. Check
   initialization, control systems, timestep behavior, and actual p changes.
2. `hp-smoke`: 20 units with splitting and joining enabled. Check successful
   adaptation, excision characteristic speeds, horizon finding, and constraints.
3. `restart-control` and `restart-checkpoint`: identical 20-unit runs, except
   the latter writes a checkpoint at 10 units, after that time's AMR phases,
   and continues. Restart that checkpoint separately with the **same binary
   and node/core layout**, using the usual `+restart` launch. Compare fields,
   mesh, subsequent AMR decisions and target shifts at matching times. Keep
   output files separate from the uninterrupted continuation.
4. Only after these pass, run `fixed`, `fixed-fine`, `pileup`, and
   `stored-target` to 200 units. `fixed-fine` adds two points per direction;
   the adaptive pair differs only in the target-shift limit (0 versus -1).
   Zero uses the previous instantaneous constraint feedback. Both adaptive
   cases retain pile-up removal and hp coarsening.

Times are relative to `Evolution.InitialTime` and correspond to M for total
mass one. Defaults: AMR every 0.5, 7–20 points per direction, h levels 1–4,
Chalmers target 1e-4, relative spectral target 1e-4, and projected C2/C3 targets
0.002/0.004. These are provisional pilot settings, **not BBH-calibrated accuracy
requirements**. Requests beyond the mesh limits stop the pilot for inspection.
Keep it in early inspiral, before the base input changes slab size near merger.

`AmrConstraints` records RMS/max gauge, two-index and three-index constraints;
`AmrFields` saves double-precision fields every 5 units. AMR statistics and
terminal decisions record cost, coarsening and target changes. The base horizon,
excision and timestep diagnostics are retained. Add `--mesh-observations` if
full mesh history every AMR call is needed; it increases output substantially.
Compare constraint peaks, horizon masses/spins, fields and work against both
fixed references. Also inspect regional constraints near each excision;
global RMS values alone can conceal local deterioration. A 200-unit pilot is
an integration check, not a waveform-accuracy or full-inspiral validation.

Local preparation on 2026-10-12: the BBH target built successfully and all eight
inputs passed native `--check-options`, using a pipeline-rendered parsing
fixture in `build-agents/chalmers-bbh-preparation`. That fixture deliberately
has no real BBH initial data and must not be used for an evolution. No BBH
simulation or checkpoint/restart comparison has been run yet.
