<!--
Distributed under the MIT License.
See LICENSE.txt for details.
-->

# Generalized-harmonic AMR experiments

Continuation of the [initial estimator and wave experiments](README.md).
For the first cluster BBH test, see [ChalmersBbh.md](ChalmersBbh.md).

## Filtered parity Chalmers AMR in stationary GH

`ChalmersSingleHole.py` now accepts `--criterion parity` to select
`ChalmersParityIgnoringHighestMode` in `EvolveGhSingleBlackHole`. Use
`--filter hypercube --half-power 256` with this variant. The default remains
`--criterion legacy`; `--criterion filtered` selects the fit that discards the
highest mode without separating parity. Both filtered estimators require the
Hypercube filter in this driver. `--h-coarsening` enables sibling joins, and
`--radial-level 1` supplies an initially over-refined radial mesh for exercising
those joins. Polynomial coarsening still requires a positive
`--coarsening-factor`.

The October 2026 GH return test uses stationary mass-one Schwarzschild in
Cartesian ingoing Kerr-Schild coordinates, analytic-Christoffel gauge, inner
excision at radius 1.8 with outgoing-characteristic checks, and exact analytic
outer ghost data at radius 6. The two radial shells meet at radius 3. The
initial data, damping, gauge, and boundaries are the same as the earlier GH
control described above. The filter is now much steeper (HalfPower 256), and
the criterion omits its highest mode and fits even and odd sequences
separately. The monitored fields are the spacetime metric, Pi, and Phi;
there is no additional GH constraint-based coarsening veto.

The over-refined pilot exposed a separate bookkeeping bug in GH sibling
joins: `control_system::FutureMeasurements` stores a NaN in its old-measurement
slot before the first measurement, indefinitely for inactive control systems.
Identical copied schedules therefore compared unequal. Its equality operator
now recognizes matching initial sentinels, while still comparing all actual
measurement times and counters. Regression tests cover initial, updated,
inactive, empty, serialized, and unequal schedules. This changes neither the
AMR estimator nor the control-system update schedule. Active moving-map
control systems have not been validated by this stationary test.

`ChalmersSingleHoleParity.yaml` is a short native input regression with
initial radial refinement, filtering, parity AMR, and coarsening enabled.
`--mesh-observations` records lapse, coordinates, element IDs, and extents in
`Meshes.vol` at every 0.1M AMR check; the full 50 GH components remain in
`Fields.vol` at the requested observation cadence. `--defer-analysis` permits
serial timed evolutions followed by a separate analysis stage.

### Stationary GH results
Eight comparison evolutions completed on ten PEs each, including fixed and
adaptive runs to 50M. The common local modal target was 1e-5, sigma threshold
0.5, coarsening factor 0.1 where enabled, and CFL 0.5 (0.25 for the timestep
control). The adequately resolved start used 14 points per direction on each
of the 12 base elements. Its measured initial full-state error was 2.1003e-7.
The two radial-coarsening controls started with 24 elements, and the coarse
start used 10 points per direction on 12 elements.

Errors below are independently integrated coordinate-volume L2 norms with
full symmetric-tensor multiplicities. They are sampled maxima: full fields
were saved every 1M in the long and 10M controls, and every 0.5M in the 5M
controls. Meshes were saved every 0.1M. The local AMR target is not a bound on
these global norms, on accumulated evolution error, or on GH constraints.

<table>
<tr>
<th>
Run
</th>
<th>
End time / M
</th>
<th>
Final elements / points
</th>
<th>
Max full-state error
</th>
<th>
Max three-index constraint
</th>
<th>
Max horizon mass error
</th>
</tr>
<tr>
<td>
Fixed grid
</td>
<td>
50
</td>
<td>
12 / 32,928
</td>
<td>
4.3896e-07
</td>
<td>
5.7089e-06
</td>
<td>
2.0563e-09
</td>
</tr>
<tr>
<td>
Parity, coarsening enabled
</td>
<td>
50
</td>
<td>
12 / 29,400
</td>
<td>
5.2759e-07
</td>
<td>
6.3773e-06
</td>
<td>
4.0036e-09
</td>
</tr>
<tr>
<td>
Filtered fit without parity separation
</td>
<td>
10
</td>
<td>
32 / 46,280
</td>
<td>
1.9748e-06
</td>
<td>
1.8490e-05
</td>
<td>
3.0156e-08
</td>
</tr>
<tr>
<td>
Radial refinement: p-coarsening only
</td>
<td>
5
</td>
<td>
24 / 45,864
</td>
<td>
4.4500e-07
</td>
<td>
5.7279e-06
</td>
<td>
4.2605e-09
</td>
</tr>
<tr>
<td>
Radial refinement: h+p coarsening
</td>
<td>
5
</td>
<td>
12 / 28,224
</td>
<td>
1.6056e-06
</td>
<td>
1.2945e-05
</td>
<td>
1.3919e-08
</td>
</tr>
<tr>
<td>
Coarse start (10 points/direction)
</td>
<td>
5
</td>
<td>
12 / 19,008
</td>
<td>
1.0119e-04
</td>
<td>
3.3978e-04
</td>
<td>
5.6301e-07
</td>
</tr>
</table>

The 50M parity run made 18 directional decreases in polynomial order, all
radial, at 0.1M and 0.2M. It saved 10.7% of the final grid points relative to
the fixed grid, with a 20.2% increase in peak full-state error and an 11.7%
increase in peak three-index constraint. Both runs settled to bounded errors
and constraints. The long run started at the minimum h-level and did not
perform any joins or splits. Its 501 mesh snapshots show no reversals.

Against the filtered fit without parity separation, on their common 0--10M
interval, the parity fit used 35.8% fewer time-averaged grid points and 36.5%
fewer final points, while reducing peak full-state error by a factor of 3.74
and peak three-index constraint by a factor of 2.90. The comparison without
parity separation made 16 splitting decisions (20 additional elements) and
178 p-decreases on surviving element IDs. It also had no immediate reversals;
its disadvantage here was unnecessary subdivision for the achieved accuracy.

The radially over-refined case made 48 p-decreases and 12 sibling joins, all
joins at 0.3M, reducing 65,856 points to 28,224 (57.1%) and 24 elements to 12.
It had no refinement reversals through 5M. With joins disabled, p-coarsening
alone retained 24 elements and 45,864 points after 102 p-decreases. Allowing
joins therefore saved a further 38.5% of the final points, but raised peak
full-state error by a factor of 3.61 and peak three-index constraint by a
factor of 2.26. This is a real accuracy cost, particularly for derivatives;
it is not evidence of constraint-preserving coarsening. At this stage the
joined trajectory was tested only to 5M, and its constraint had not fully
flattened. The follow-up below extends that trajectory and isolates the
join, filter, and resolution effects before a dynamical black-hole test.

From the coarse start, the criterion made 60 p-increases and no splits or
reversals, reaching 19,008 points by 0.2M. This demonstrates refinement, but
not recovery to the accuracy of the well-resolved runs: global full-state
error started at 7.8705e-5, peaked at 1.0119e-4, and was still 3.8224e-5 at
5M. Both initial/evolution error and the looser final mesh permitted by
hysteresis can contribute; this experiment does not separate them. A modal
truncation estimate should not be treated as a repair mechanism for an
under-resolved start or as a global accumulated-error estimator.

The refinement-only 10M control retained the fixed mesh. All 101 mesh
snapshots matched the fixed-grid prefix, and the full nodal GH state was
identical at 10M (integrated difference exactly zero). The half-CFL control
also matched all 101 parity-AMR meshes. Its full-state difference at 10M was
5.2768e-11, only 0.0100% of the reference error. This is a combined timestep
and substep-filter-frequency check, not a formal temporal convergence order
measurement. The joined versus unjoined radial controls differed by
1.5852e-6 at 5M when integrated on their common partition.

All independent diagnostics passed the 20/28-point quadrature comparison;
the largest change in any reported diagnostic was below 2.6e-15. Mapped
coordinates,
domain volume, initial nodal metric, initial discrete Phi constraint, complete
observation times, and mesh/statistics agreement were checked. Native gauge
and four-index constraint norms are retained in each `Analysis.json`.

Nine targeted CTests passed: parsing and executing the original GH AMR input,
its coarsening input, and the new parity/join input; the FutureMeasurements
and control-system Trigger unit tests; and the Chalmers Python tests. Source
formatting and `git diff --check` also passed. The new input regression ran
on one PE; the comparison evolutions ran on ten PEs.

The complete local records are in `build-agents/chalmers-gh-parity/`:
`Summary.json`, `CommonWindow10.json`, `Comparisons.json`, `Validation.json`,
`BuildMetadata.json`,
the per-run inputs and HDF5 data, `LongComparison.png`, `Adaptation.png`,
`run_study.py`, and `analyze_study.py`. Source snapshots and a patch relative
to the start of this study are retained there. The initial empty-field
output-option rejection and pre-fix join assertion are preserved as setup
failures. A small auxiliary regression was deliberately stopped after it
briefly overlapped the start of the fixed-grid run; see `TimingNotes.json`.
Wall times are single-run desktop measurements, and that fixed-grid timing
is not an uncontaminated benchmark. Grid counts are the primary cost measure.

A representative run can be reproduced after building
`EvolveGhSingleBlackHole` and refreshing the Python package:

```sh
build-agents/bin/python-spectre tools/Amr/ChalmersSingleHole.py run \
  --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/gh-parity-reproduction \
  --strategy hp05 --criterion parity --filter hypercube --half-power 256 \
  --points 14 --tolerance 1e-5 --coarsening-factor 0.1 --h-coarsening \
  --cfl 0.5 --end-time 50 --observation-interval 1 --processes 10 \
  --mesh-observations --defer-analysis --timeout 14400
build-agents/bin/python-spectre tools/Amr/ChalmersSingleHole.py analyze \
  build-agents/gh-parity-reproduction
```

Add `--radial-level 1 --end-time 5 --observation-interval 0.5` for the join
control. Omit `--h-coarsening` to isolate p-coarsening, or set
`--coarsening-factor 0` and omit `--h-coarsening` for refinement only. Select
`--criterion filtered` for the comparison without parity separation.

### Separating GH join transients from resolution error

The follow-up records are in `build-agents/chalmers-gh-joins/`. The native
executable and AMR criterion are unchanged from the preceding study.
`ChalmersSingleHole.py` adds three diagnostic/control options:

- `--grid-points Nx Ny Nz` overrides the scalar `--points` setting and permits
  an independently initialized control on the final anisotropic parent mesh.
- `--early-observation-until T` saves full fields at the existing 0.1M AMR
  boundaries through T, preserving the evolution time-step boundaries.
- `--constraint-diagnostics` also records native nodal maxima of the gauge,
  three-index, and four-index constraint norms.

The radially refined start contains 24 elements and 65,856 points. At
coarsening factor 0.1, all 12 radial sibling pairs join at 0.3M after 48
directional p-decreases. The resulting 12-element mesh has 14 x 14 x 12
points per element. A third, fixed-mesh control starts independently from
analytic data on precisely that parent mesh, with the usual discrete
initial Pi/Phi constraint setup. Joined, unjoined, and fixed-parent controls
all run to 50M, with ten PEs per evolution and no simultaneous evolutions.

<table>
<tr>
<th>
Control
</th>
<th>
Final elements / points
</th>
<th>
Peak full-state error
</th>
<th>
Peak three-index constraint
</th>
<th>
Grid-point-steps
</th>
</tr>
<tr>
<td>
h+p coarsening, factor 0.1
</td>
<td>
12 / 28,224
</td>
<td>
1.60561e-6
</td>
<td>
1.34124e-5
</td>
<td>
231,465,024
</td>
</tr>
<tr>
<td>
p-coarsening only, factor 0.1
</td>
<td>
24 / 45,864
</td>
<td>
4.47092e-7
</td>
<td>
5.76446e-6
</td>
<td>
508,944,576
</td>
</tr>
<tr>
<td>
Fixed parent mesh
</td>
<td>
12 / 28,224
</td>
<td>
1.60487e-6
</td>
<td>
1.34124e-5
</td>
<td>
225,792,000
</td>
</tr>
</table>

The errors are sampled coordinate-volume L2 norms, with full tensor
multiplicities and independent 20/28-point quadrature checks. Native
`L2Norm` diagnostics are unweighted nodal RMS values and are retained
separately. Full fields are sampled every 0.1M through 1M and every 1M
thereafter; meshes are sampled every 0.1M. Grid-point-steps sum the grid
size at each GTS step, excluding the final endpoint. They are a cost proxy,
not measured CPU work or a controlled wall-time benchmark.

The joined and fixed-parent runs settle to the same discrete steady state.
Their full-state difference decreases from 3.92079e-7 at 1M to 3.94176e-8
at 5M and 1.80960e-12 at 50M. The last difference is 1.13e-6 of the
fixed-parent error. The joined run's state error and three-index constraint
each vary by less than 2e-7 relative over the sampled 40--50M interval.
There are no further mesh changes after the initial joins and no refinement
reversals. These controls attribute the persistent accuracy cost to the
coarse filtered mesh; the additional join-history difference decays.
The factor-0.1 join saves 54.5% of grid-point-steps relative to its p-only
control, but its larger solution and constraint errors remain a real cost.

An exact-state replay separates the immediate operators. A control ending
at 0.3M stops before the final-time AMR action, and native child-to-parent
projection matrices reproduce the longer run's post-join state exactly.

| Stage | Full-state error | Three-index constraint |
| --- | ---: | ---: |
| Before joining | 3.16899e-7 | 1.56325e-6 |
| Native projected parent state | 3.52199e-7 | 1.67930e-6 |
| One native filter on that state, offline | 8.69702e-7 | 3.52498e-6 |
| Evolved to 0.4M | 1.02110e-6 | 3.20315e-6 |

The raw join increases the three-index constraint by 7.4%. Applying one
actual alpha-36, HalfPower-256 filter offline increases it a further factor
of 2.10, predominantly in the outer shell. This probe is not an evolved RK
substep. Together with the fixed-parent control, it identifies the
interaction between filtering and the reduced parent resolution as an
important source of the observed accuracy loss.

### Calibrating the GH coarsening budget

Paired 5M screening runs vary the coarsening factor at fixed refinement
target 1e-5, with and without h joins. Factors 0.03 and 0.01 reject every
join but still reduce polynomial order. For each of these factors, all 15
sampled full states are bitwise identical between the h-enabled and p-only
runs. Their final grids contain 49,392 and 50,568 points, respectively.

Replaying the existing criterion on the common pre-join state suggests an
intermediate factor of 0.05. The worst component's estimated parent error
is 4.3594e-7 in the inner shell, dominated by an angular tail, and 7.7762e-7
in the outer shell, dominated by a radial tail. The 5e-7 join budget should
therefore accept only the six inner-shell parents. The native run makes
exactly those six joins at 0.3M and retains radial refinement outside r=3.

<table>
<tr>
<th>
Factor-0.05 control, 0--5M
</th>
<th>
Final elements / points
</th>
<th>
Mean points
</th>
<th>
Peak state error
</th>
<th>
Peak three-index constraint
</th>
<th>
Grid-point-steps
</th>
</tr>
<tr>
<td>
h+p coarsening
</td>
<td>
18 / 38,808
</td>
<td>
40,172.2
</td>
<td>
4.42280e-7
</td>
<td>
5.83787e-6
</td>
<td>
36,230,208
</td>
</tr>
<tr>
<td>
p-coarsening only
</td>
<td>
24 / 49,392
</td>
<td>
50,144.6
</td>
<td>
4.34730e-7
</td>
<td>
5.66155e-6
</td>
<td>
69,769,728
</td>
</tr>
</table>

Enabling the selective joins saves 19.9% of time-averaged points and 48.1%
of grid-point-steps, with increases of 1.74% in peak state error and 3.11%
in peak three-index constraint. The h+p case makes 66 directional
p-decreases as well as the six joins; the p-only case makes 84
p-decreases. Neither has splits, p-increases, or reversals. These are
matched short-window comparisons, not a general optimality claim for 0.05.
The sigma threshold remains 0.5; the parameter varied here is the
coarsening budget factor, not the h-versus-p smoothness threshold.

The selected factor-0.05 h+p run also completes 50M. Its final state error
is 4.44194e-7 and its final three-index constraint is 5.86709e-6; their
sampled peaks are 4.44354e-7 and 5.87605e-6. Over 40--50M their relative
ranges are 2.45e-6 and 9.69e-7, respectively. It retains 18 elements and
38,808 points after 0.4M, with no subsequent changes or reversals, and uses
315,647,808 grid-point-steps. All 15 full-state observations and all 51
meshes through 5M exactly match the short run. The matched factor-0.05
p-only control is only 5M long; the 50M p-only control above uses factor
0.1, so a comparison of those long runs changes both settings.

At factor 0.05, tightening the refinement target to 3e-6 or 1e-6 rejects
all joins and retains more radial points, but produces no p-increases.
Their 5M peak state errors are 4.33746e-7 and 4.33678e-7, with peak
three-index constraints 5.65272e-6 and 5.65236e-6. Both retain the initial
angular resolution. This is an accuracy plateau under those two target
changes, not evidence of target-proportional convergence.

A further 5M test at target 1e-7 breaks that plateau. It makes 36 angular
p-increases and 42 radial p-decreases, with no joins, splits, or reversals.
The final grid has 64,062 points, 2.7% fewer than the initial grid despite
the angular refinement. Final state error is 1.87328e-7, sampled peak
state error is 2.91313e-7 (at 0.1M), and sampled peak three-index constraint
is 1.53153e-6. This verifies a response to a sufficiently tight target and
simultaneous refinement and coarsening in different directions. It is a
short screening run, not a convergence-order measurement or a global
error guarantee.

The operator replay also tests hypothetical parent radial extents 13 and
14, retaining angular extents 14. After one filter, their state errors are
3.52199e-7 and 3.18205e-7, and their three-index constraints are 1.67930e-6
and 1.56765e-6. One extra radial point therefore largely avoids the
immediate filter-induced jump seen with extent 12. These are offline
projections only: the native join still takes the maximum child extent
and does not increase the parent's polynomial order. A useful next
algorithmic experiment is to consider a join with a simultaneous p-increase
and evaluate the resulting filtered parent against the same budget.
These stationary controls do not currently motivate a new GH-specific
constraint veto, and they do not validate dynamical horizons or binary
black holes.

All 15 evolutions in this follow-up completed successfully. The final
checks verify finite diagnostics, complete field/constraint/horizon
observations, mesh statistics, unchanged executable and driver hashes,
and quadrature agreement (maximum diagnostic change below 2.6e-15).
Native and independently computed final nodal C3 maxima agree to better
than 3.1e-14 in the selected long run. The GH parity-input parsing CTest,
driver formatting/import checks, and whitespace checks pass. One launch
failed during Charm++ startup before any evolution; its input and log
are preserved separately, and an unchanged retry succeeded.

`build-agents/chalmers-gh-joins/Results.md` contains the numerical tables
and four inspected plots. `SummaryAll.json`, `JoinReplay.json`,
`JoinBudget.json`, `CalibrationValidation.json`, `FinalValidation.json`,
the run manifests, inputs, native logs, HDF5 data, and source snapshots
retain the evidence. Only the diagnostic Python driver and this
documentation changed in this follow-up; the native criterion is unchanged.

To reproduce the selective-join control with the additional diagnostics:

```sh
build-agents/bin/python-spectre tools/Amr/ChalmersSingleHole.py run \
  --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/gh-selective-join-reproduction \
  --strategy hp05 --criterion parity --filter hypercube --half-power 256 \
  --points 14 --radial-level 1 --tolerance 1e-5 \
  --coarsening-factor 0.05 --h-coarsening --cfl 0.5 \
  --end-time 50 --observation-interval 1 --early-observation-until 1 \
  --constraint-diagnostics --mesh-observations --processes 10 \
  --defer-analysis --timeout 14400
build-agents/bin/python-spectre tools/Amr/ChalmersSingleHole.py analyze \
  build-agents/gh-selective-join-reproduction
```

## Dynamical Schwarzschild gauge relaxation

`ChalmersGaugeRelaxation.py` starts a mass-one, nonspinning hole with
Kerr-Schild spatial geometry and damped harmonic gauge from the initial time.
The GH initializer sets Pi and Phi consistently with the gauge and discrete
three-index constraint. The ensuing change is a coordinate relaxation;
differences from the stationary Kerr-Schild spacetime components are therefore
not evolved errors. The driver removes that analytic error observer.

The static cubed-sphere domain extends from radius 1.2 to 20, with logarithmic
shells at 2, 4, and 8. Excision uses `DemandOutgoingCharSpeeds`. The outer
condition is `ConstraintPreservingBjorhus`, type
`ConstraintPreservingPhysical`, with no incoming wave profile, as used in the
BBH input. The damped harmonic amplitudes are `[1,0,1]`, exponents `[2,2,2]`,
and spatial decay width is 17.0152695482514. The Kerr-Schild template supplies
positive gamma2 damping. Filtering uses `Hypercube` with half-power 256 on
every substep. The parity-separated directional-energy criterion ignores the
highest mode, has sigma threshold 0.5, and permits p and optional h coarsening
with budget factor 0.05. Global Dormand-Prince 5 time stepping, CFL 0.5,
and AMR every 0.1M match the earlier stationary controls.

`ChalmersGaugeRelaxationAnalysis.py` measures the three-index constraint on
independent Gauss grids, compares evolved numerical solutions on their common
element partition, checks excision outflow, and records horizon mass, spin,
and coordinate radius. Its tensor norms use coordinate-volume integration
and symmetric-component multiplicities; native `L2Norm` observations instead
use nodal RMS. `ChalmersMeshHistory.py` counts realized splits, joins, and
p changes from mesh snapshots at each AMR check, including balancing effects.
These counts do not represent unfulfilled requests from the criterion.

The small input
`tests/InputFiles/GeneralizedHarmonic/ChalmersGaugeRelaxation.yaml` exercises
the gauge, physical outer boundary, filter, and AMR combination through 0.2M.
For a longer matched comparison, use separate output directories:

```sh
build-agents/bin/python-spectre tools/Amr/ChalmersGaugeRelaxation.py \
  --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/gh-relaxation-hp \
  --strategy hp05 --points 12 --radial-level 1 \
  --tolerance 1e-4 --coarsening-factor 0.05 --h-coarsening \
  --end-time 5 --processes 10
build-agents/bin/python-spectre tools/Amr/ChalmersGaugeRelaxationAnalysis.py \
  analyze build-agents/gh-relaxation-hp
```

For p-only, replace `--strategy hp05 --h-coarsening` with `--strategy p`.
For a smaller initial mesh, use `--grid-points 12 12 8 --radial-level 0`.
Fixed controls use `--strategy fixed --points 12` or `--points 14` with radial
level zero. Reference convergence and time-step convergence must be checked
before treating a fixed run as an accuracy reference. No moving maps or active
horizon controls are exercised by this setup.

For this exactly spherical setup, a fixed `CartoonSphere1D` calculation can
check the numerical reference much more cheaply than increasing every 3D
angular order. `ChalmersGaugeRelaxationReference.py` uses
`EvolveGhNoBlackHole3D`: it shares the GH initializer and evolution actions with
the single-hole executable, and supports the Cartoon domain without requiring
horizon-control anchors. It retains the same physical initial data, gauge,
outer boundary, damping, and filter. It removes horizon interpolation events
and has no AMR. The cubed-sphere experiments continue to use the single-hole
executable and observe the horizon.

`ChalmersSphericalReference.py` reconstructs all Cartesian GH components from
14 parity-even radial coefficient functions. It verifies spherical symmetry
on the axis and the logarithmic maps, then compares with the 3D solution on a
common radial partition. Independent radial differentiation supplies a
three-index-constraint norm. A spherical areal-radius mass diagnostic is also
recorded. The reference must be checked for both spatial and time-step
convergence; a small horizon-mass error alone is insufficient.

For example, the following fixed reference uses 16 radial points per element
and refinement levels `[2,3,4,5]` in the four shells. Compare with `[1,2,3,4]`
and with halved CFL before interpreting its differences as accuracy estimates:

```sh
cmake --build build-agents --target EvolveGhNoBlackHole3D -j 2
build-agents/bin/python-spectre tools/Amr/ChalmersGaugeRelaxationReference.py \
  --executable build-agents/bin/EvolveGhNoBlackHole3D \
  --output build-agents/gh-relaxation-reference \
  --points 16 --radial-levels 2 3 4 5 --cfl 0.5 \
  --end-time 5 --processes 8
build-agents/bin/python-spectre tools/Amr/ChalmersSphericalReference.py \
  analyze build-agents/gh-relaxation-reference
build-agents/bin/python-spectre tools/Amr/ChalmersSphericalReference.py \
  compare-3d build-agents/gh-relaxation-hp \
  build-agents/gh-relaxation-reference \
  --output build-agents/gh-relaxation-hp-reference.json
```

The initial seven-run matrix completes 5M in every cubed-sphere case. The
smaller-start hp run increases radial p in the innermost shell at 0.3M, splits
six outer-shell elements at 0.6M, increases radial p in the 2--4M shell at
0.7M, and lowers both angular orders in the outermost children at 2.6M.
It makes six splits, twelve directional p increases, and twelve p decreases.
The initially finer hp run at target 1e-4 instead makes six radial joins in
4 < r < 8 and 138 p decreases. At target 3e-5 it rejects those joins,
increases angular p in six elements, and splits six outer-shell elements
at 4.6M. No run has an immediate refinement/coarsening reversal. These are
realized mesh changes; they do not establish a complete late-time h
split/join cycle.

The finest fixed spherical reference has 60 radial elements with 16 points
each and levels `[2,3,4,5]`. Refining from `[1,2,3,4]` changes its full state
by at most 2.87267e-4 in coordinate-volume L2. The final cubed-sphere error
estimates change by at most 0.0291% when those two references are exchanged.
The measured final state differences and work are:

| Run | Final state L2 difference | Grid-point-steps |
| --- | ---: | ---: |
| Fixed 12^3 | 1.29310e-2 | 71,788,032 |
| Fixed 14^3 | 9.92614e-3 | 156,605,568 |
| Initially finer hp, target 1e-4 | 6.17195e-3 | 112,276,296 |
| Matched p-only, target 1e-4 | 6.17210e-3 | 123,092,712 |
| Initially finer hp, target 3e-5 | 5.89304e-3 | 129,636,450 |
| Smaller-start hp, target 1e-4 | 1.39386e-2 | 32,686,272 |
| Matched smaller-start p-only | 2.12105e-2 | 28,133,568 |

Allowing the six selective joins saves 8.79% of grid-point-steps with almost
unchanged reference error; the direct final hp-minus-p state difference is
4.26674e-5. From the smaller starting mesh, hp reduces final state error by
34.3% while using 16.2% more grid-point-steps. These are two different accuracy
and cost tradeoffs. Local wall times are retained in the study, but are not
controlled benchmarks and do not scale consistently with this work proxy.
Tightening the target from 1e-4 to 3e-5 improves final state error by only
4.52%, rather than in proportion to the target change.

Constraints require a separate assessment. The smaller-start hp run's final
integrated three-index constraint is 9.40615e-3, versus 2.99055e-3 for p-only.
The excess is dominated by the radial constraint in the outer shell and
already exists before angular coarsening. Its sampled peak is nevertheless
slightly lower than p-only's, so final and peak norms should both be retained.
The fixed 14^3 run also has a larger final integrated constraint than fixed
12^3 despite smaller state error and much smaller near-hole nodal peaks.
A single constraint sample, nodal peak, or horizon mass is therefore
insufficient to rank these meshes.

The sampled excision light-speed margin stays above 0.34093. Maximum horizon
mass error is 2.69e-9 for fixed 14^3, about 6.52e-8 for the initially finer
adaptive runs, and 5.74e-7 for the smaller-start controls. The horizon's
coordinate radius grows from 2 to about 2.062M. This is a successful short
gauge-transient and adaptation test, with an unresolved constraint-versus-
state-error tradeoff. Before a BBH pilot, isolate that outer-shell radial
constraint error and test a longer relaxation interval.

`build-agents/chalmers-gh-relaxation/Results.md` contains the complete tables,
reference checks, inspected figures, and validation record. The run directories
retain all inputs, output, and failed reference startup attempts. The initial
single-hole Cartoon launch lacked a required horizon-control grid anchor;
the supported no-horizon executable resolves that issue. A 24-point reference
launch exceeded this build's 21-point Legendre limit; the completed references
instead refine h at 16 points per element.

Halving CFL on the finest graded spherical reference changes the full state
by at most 8.41e-10; the uniform level-3 check gives 1.70e-9. Time-step error
is negligible at the scale of the measured spatial errors. The input parsing
and execution CTests, Python formatting/import checks, analytic diagnostic
checks, finite-data/completeness checks, and quadrature checks all pass.

## Directional constraint feedback in the GH hp decision

`gh::ChalmersWithConstraintFeedback` wraps the filtered, parity-separated
Chalmers criterion in the single-hole executable. It evaluates the current
three-index reduction constraint along each normalized column of the physical
Jacobian:

\[
 C_d = \sqrt{\frac{1}{N}\sum_{q=1}^{N}\sum_{a,b}
 \left[\frac{J^i{}_d(\partial_i g_{ab}-\Phi_{iab})}
 {\sqrt{\sum_i(J^i{}_d)^2}}\right]^2_q}.
\]

This is a nodal RMS of Euclidean tensor components. It includes the two copies
of off-diagonal spacetime components. It is not a physical-volume integral or
an invariant norm. The existing `Constraints` criterion also projects the
derivative index, but sums only independent tensor components; its numerical
targets therefore should not be copied without adjusting the normalization.

For spectral target `T` and constraint target `C_target`, the feedback uses
`T_d = T * max(MinimumTargetFactor, min(1, C_target/C_d))`. Zero constraint
gives factor one. Both absolute and relative spectral targets are scaled.
The smallest directional factor is used for the common fit-negligibility
budget and the prospective-parent field budget. Directional refinement and
p-coarsening compare against their own scaled targets.

The controller has three independently reproducible settings:

- `baseline`: factor one and no constraint veto. The extra constraint
  calculation and logging should leave the original mesh decisions unchanged.
- `veto`: factor one; reject p reduction in a direction whose current constraint
  exceeds the target. Reject a join if any direction is excessive, since parent
  projection/filtering can affect all derivatives. The usual sibling agreement
  requires both children to accept a join.
- `tighten`: apply the same veto and the bounded directional target adjustment.
  This first experiment is instantaneous and has no history-dependent state,
  temporal predictor, or delayed relaxation of the target.

All refinement flags are preserved by the veto. This is a check on the current
constraints, not a calculation of prospective constraint injection. The
existing field projection/filter-loss budgets still apply. A low current
constraint alone does not certify the constraint of the coarsened state.

Example using the smaller starting mesh from the gauge-relaxation study:

```sh
build-agents/bin/python-spectre tools/Amr/ChalmersGaugeRelaxation.py \
  --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/gh-feedback-tighten \
  --strategy hp05 --grid-points 12 12 8 --radial-level 0 \
  --tolerance 1e-4 --coarsening-factor 0.05 --h-coarsening \
  --constraint-feedback tighten --constraint-target 3e-5 \
  --minimum-target-factor 0.1 --end-time 5 --processes 10
build-agents/bin/python-spectre tools/Amr/ChalmersGaugeRelaxationAnalysis.py \
  analyze build-agents/gh-feedback-tighten
build-agents/bin/python-spectre tools/Amr/ChalmersGhFeedbackAnalysis.py \
  build-agents/gh-feedback-tighten
```

The feedback analysis writes per-element directional constraints, native veto
counts, converged outer-shell radial/tangential constraint integrals, and
compressed spectral diagnostic records. The latter contain the filtered,
parity-aware Chalmers estimate, the weighted spectral-tail estimate, and the
existing SpEC-style multi-fit slope and pile-up count. The ordinary multi-fit
slope is not parity-aware: these records are diagnostic only and do not change
the controller. The highest filtered mode is excluded from those estimates.
Full power spectra are retained so parity zeros and small-amplitude components
can be distinguished from a meaningful plateau.

Native feedback logs are evaluated before AMR. Volume observations at the
same time follow AMR; a p change preserves the element ID but changes its
nodes and hence its nodal RMS. Compare native and reconstructed constraints
only on unchanged elements. Use native logs for actual decision-time factors
and the volume spectra to assess the resulting representations.

The three matched runs complete 5M on ten processes. They use the same
smaller-start gauge-relaxation problem, spectral target 1e-4, constraint target
3e-5, and coarsening factor 0.05. Tightening is limited to a factor of ten.
The final norms below are coordinate-volume L2 norms; the state difference
uses the previously converged, fixed spherical reference. Work is measured
in grid-point-steps, not wall time.

<table>
<tr>
<th>
Feedback
</th>
<th>
Final three-index constraint
</th>
<th>
Final state difference
</th>
<th>
Grid-point-steps
</th>
</tr>
<tr>
<td>
Baseline
</td>
<td>
9.40615e-3
</td>
<td>
1.39386e-2
</td>
<td>
32,686,272
</td>
</tr>
<tr>
<td>
Coarsening veto
</td>
<td>
9.40615e-3
</td>
<td>
1.39386e-2
</td>
<td>
32,686,272
</td>
</tr>
<tr>
<td>
Veto and target tightening
</td>
<td>
8.10147e-3
</td>
<td>
1.42440e-2
</td>
<td>
33,479,256
</td>
</tr>
</table>

No constraint veto fires in these runs. Baseline and veto have identical mesh
histories, and the new baseline reproduces the original run's fields and
coordinates exactly. Thus this evolution does not demonstrate a benefit from
the veto; the unit tests exercise its rejection paths.

Tightening adds twelve angular p increases in the innermost shell at 2.5M
and six radial p increases in the outer children at 3.9M. All three runs retain
the six outer-shell splits and twelve angular p reductions. None has an
immediate refinement/coarsening reversal. The minimum realized target factor
is 0.180054, so the factor-of-ten bound is not reached.

At 5M, tightening lowers the integrated constraint by 13.87% at 2.43% more
work, but increases the state error by 2.19%. The early sampled peak constraint
is unchanged at 9.43816e-3. The final outer-shell radial constraint falls from
9.39972e-3 to 8.09518e-3; the tangential part changes from 1.68963e-4 to
1.70793e-4. Exchanging the finest spherical reference for the next-coarser
reference changes the tightened run's final state-error estimate by only
0.00176%, much less than the measured deterioration. All runs retain an
excision light-speed margin above 0.34093 and a maximum horizon mass error
of 5.74e-7.

A separate 0.5M pilot starts with 48 elements and 12^3 points per element,
with tightening and the veto enabled. It completes six h joins and 126
directional p reductions, with no immediate reversals. This establishes that
both forms of coarsening remain executable when the constraints permit them;
it does not establish their long-term accuracy.

The feedback mechanism works, but this short experiment does not show an
overall accuracy improvement. Directional constraints can remain above their
target because tightening the field target is an indirect response, not an
enforced constraint bound. A longer evolution after the extra radial
refinement is needed before judging whether its benefit persists. A subsequent
comparison could test persistent target adjustment or prospective constraint
injection, while retaining separate state-error and coarsening checks.

`build-agents/chalmers-gh-feedback/` retains inputs, native decisions, passive
spectra, the exact run executable, source snapshots, reference comparisons,
and `FeedbackComparison.png`/`.pdf`. `Summary.json`, `Comparisons.json`, and
`Validation.json` hold the aggregate measurements. All eleven targeted CTests
pass, including the new input parse/execute tests. Reconstructed directional
constraints agree with the six-significant-digit native logs to 4.08e-6
relative on unchanged elements. Independent quadrature checks and reference
sensitivity checks support the reported differences.

## Extending the GH feedback comparison to 20M

The extension in `build-agents/chalmers-gh-feedback-long/` uses the same saved
single-hole executable and the same physics, field target, constraint target,
filter, coarsening factor, and mesh limits as the 5M comparison. The inputs
change only the completion time and the associated AMR, observation, and
step-to-time sequences. The original runs had no restart checkpoints, so the
extension starts from the same initial data. Full field observations continue
at 1M intervals, with the original early observations retained; mesh snapshots
and criterion evaluations continue every 0.1M.

The three-dimensional runs use three PEs for baseline and tightening and two
for the replacement veto run. Two one-PE references fill the remaining slots;
after the coarser reference finishes, its PE runs the half-CFL check. A veto
launch on three PEs failed inside Charm++ startup before application
initialization. Its crash report and output are retained, and its replacement
uses identical input options. Analysis replaces completed simulations with at
most two workers per main run, keeping the total allocation within ten.

Both fixed spherical references use 16 points per element, with radial levels
`[1,2,3,4]` and `[2,3,4,5]`. Their full-state difference at 20M is 1.80477e-5;
its maximum over the observed interval is 2.87267e-4. Halving CFL from 0.5 to
0.25 on the coarser reference changes the state by at most 2.95505e-9 through
20M. This is a temporal check of the reference, not a separate time-step
convergence test of each adaptive run. The adaptive-run error estimates are
compared against both references at 10M and 20M.

All three main runs complete 20M and give the following results.
Constraint and state norms are coordinate-volume L2 norms, with the state
error measured against the finer spherical reference.

<table>
<tr>
<th>
Strategy
</th>
<th>
Final three-index constraint
</th>
<th>
Final state error
</th>
<th>
Grid-point-steps
</th>
<th>
Final elements / points
</th>
</tr>
<tr>
<td>
Baseline
</td>
<td>
4.29977e-3
</td>
<td>
6.10163e-2
</td>
<td>
177,171,072
</td>
<td>
48 / 55,680
</td>
</tr>
<tr>
<td>
Coarsening veto
</td>
<td>
4.29977e-3
</td>
<td>
6.10163e-2
</td>
<td>
177,171,072
</td>
<td>
48 / 55,680
</td>
</tr>
<tr>
<td>
Veto and target tightening
</td>
<td>
4.64448e-3
</td>
<td>
6.53114e-2
</td>
<td>
221,367,312
</td>
<td>
42 / 60,162
</td>
</tr>
</table>

No coarsening veto fires throughout the extension. Baseline and veto have
identical fields, coordinates, and element extents at all 26 full-field
observations, including 20M. Their mesh-change counts and work are identical.

At 20M, tightening increases the final constraint by 8.02%, the state error
by 7.04%, and work by 24.95%. Its state error was 3.88% smaller at 10M, but
2.38% larger at 15M: the apparent improvement was temporary. Exchanging the
two spherical references changes the error estimates at 10M and 20M by less
than 0.00044%, well below these differences. The sampled peak integrated
constraint remains 9.43816e-3 in both runs, attained during the early transient.

The baseline makes 24 splits and 18 directional p increases; tightening makes
18 splits and 90 p increases. Both split at 10.2M, including refinements
required by balance. Baseline additionally splits the six outermost children
at 12M. Tightening instead increases the angular extents near the hole, reaching
17 points, below the configured maximum of 20. It takes 4,535 global steps,
compared with baseline's 3,900, as well as using more points per step on
average.
These figures account for the increased work despite its smaller element count.

The final inner-shell integrated constraint decreases from 3.20636e-3 to
2.95872e-3, while the outer-shell radial part increases from 2.84635e-3 to
3.56145e-3 (25.1%). This identifies an important difference in where the two
strategies spend resolution; it does not by itself establish the cause of
the evolved state-error difference. Maximum horizon mass errors are very
similar, 1.64543e-6 and 1.64226e-6.

Both runs retain all twelve angular p reductions at 2.6M and add six radial
p reductions at 18.9M. Neither has an immediate reversal detected by the
mesh-history checks. There are no h joins in these initially coarse runs;
the earlier finer-start pilot remains the h-coarsening test.

The target factor reaches its lower bound of 0.1 near the hole at 6.6M and
remains bounded there thereafter. The tightened run's largest final directional
nodal constraint is about 117 times `ConstraintTarget`. The indirect target
adjustment therefore does not enforce that constraint target. This longer
test supports keeping the instantaneous capped feedback experimental. The
next useful diagnosis is its allocation
of p refinement near the hole and h refinement in the outer shell, with both
state error and constraint error retained as separate measures.

All six successful evolutions reach 20M, including the three reference runs.
The three-dimensional fields and both original-CFL references reproduce the
earlier outputs through 4M exactly. Independent constraints, quadrature,
finite-data and observation-completeness checks pass. The reconstructed
directional constraints agree with the rounded native logs within 4.81e-6
relative. The sampled excision light-speed margin remains positive throughout.
No C++ algorithm changes or rebuilds are made for this extension.

The study directory retains the inputs and launch/analysis scripts, native
decisions and spectra, `Summary.json`, `Comparisons.json`, `LateTime.json`,
`BaselineVetoEquivalence.json`, and `FinalValidation.json`. The inspected
figures are `FeedbackComparison20M.png` and `LateTimeTradeoffs.png`, each also
available as PDF. They show the full evolution and the changes relative to
baseline after 5M, respectively.

## Auditing the saved 20M spectra for pile-up

`ChalmersPileUpAudit.py` audits the two distinct trajectories in the 20M
feedback study. The veto trajectory is already verified identical to baseline
at all 26 field observations. This audit changes no evolution code or AMR
decisions and launches no simulations:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  build-agents/bin/python-spectre tools/Amr/ChalmersPileUpAudit.py \
  build-agents/chalmers-gh-feedback-long/baseline-20 \
  build-agents/chalmers-gh-feedback-long/tighten-20 \
  --output build-agents/chalmers-gh-pileup-audit
```

The audit reconstructs the transverse Legendre-L2 energy spectra used by
Chalmers from the saved fields. It removes the highest filtered mode and
examines each parity separately, using the criterion's negligibility and
roundoff conventions. The ordinary SpEC diagnostic on the saved power
monitors is retained for comparison; its raw count must not be interpreted
as a count of meaningful plateau modes in a spectrum with parity zeros.

An exploratory two-window screen requires at least five significant samples
in one parity. It fits a prefix ending at the first of the final three samples,
trying prefix starts 0, 1 and 2 within that parity. Among positive fits spanning
at least 0.25 decades it selects the smallest RMS residual. The default screen
requires a prefix drop of at least 0.5 decades, a prefix fit RMS at most 0.2
decades, a terminal slope magnitude at most 20% of the prefix slope, and a
terminal range at most 0.5 decades. This is an audit heuristic, not a port of
the SpEC controller or a calibrated refinement rule. Looser and stricter
screens are also reported. Short spectra remain inconclusive.

| Diagnostic | Baseline | Target tightening |
| --- | ---: | ---: |
| Field/direction spectra examined | 137,700 | 129,600 |
| Ordinary SpEC pile-up count greater than 3 | 388 | 744 |
| Parity-aware plateau candidates, default screen | 0 | 1,704 |
| Candidates, stricter / looser screens | 0 / 0 | 600 / 3,048 |
| Element/direction observations containing default candidates | 0 | 120 |
| Candidate directions passing existing p-coarsening checks | 0 | 0 |

Every raw SpEC flag has a negligible parity. All 388 baseline flags and 482
of the tightening flags also have total resolved directional variation below
the common refinement budget. The raw counts alone are consequently poor
evidence for excessive resolution. Neither the raw flags nor the plateau
candidates coincide with an above-target Chalmers estimate in the saved states.

The tightening candidates occur in the two angular directions of the six
innermost elements, at the observations from 11M through 20M, with extents 15
or 17. Some persist across five consecutive 1M observations at unchanged
extent. For example, at 13M, block 0's direction-0 `SpacetimeMetric_xx` has
an even-parity prefix slope of 0.736 decades per degree and a terminal slope
of 0.112. Its Chalmers error estimate is 0.768 times the 1e-5 target. This is
resolved decay followed by a substantially flatter tail, not an alternating
parity-zero artefact. Its origin is not established by the spectrum alone.

The native projection/filter-aware p-coarsening routine rejects every one
of the 120 affected element/direction observations after checking all 50
components. The largest component's candidate error is 68.6--108.4 times the
directional coarsening budget (target times 0.05/3). Thus a plateau is not a
license to discard those modes while preserving the current accuracy budget.

Decision-time coverage is limited. AMR occurs every 0.1M, whereas fields are
saved at 26 times, after AMR. The audit reproduces 4,860 native directional
decisions on elements whose meshes are unchanged; all are `DoNothing`.
All saved states are below the refinement threshold, so the hypothetical
p-veto changes zero snapshot decisions. Crucially, none of the native 18
baseline or 90 tightening p-increase requests has an exact saved pre-AMR
spectrum. The preceding 15M observation contains plateaus in the directions
of the 24 tightening requests at 15.8M and 15.9M, but the field driving those
requests is unknown. The bracketing interval can be 1M. These observations
cannot establish which actual p increases would have been prevented, or any
resulting savings or accuracy change.

Nine analytic classification checks and three directional Parseval checks
pass. Reconstructed power monitors agree with the saved native values within
9e-16 relative, and reconstructed Chalmers errors and slopes within 1.4e-15
in the reported target-scaled/absolute measures. These checks validate the
reconstruction, not the empirical plateau thresholds. Outputs include
`Summary.json`, per-run `Audit.json`, `SpectraAudit.json.gz`, `Directions.json`
and `FollowUp.json`, plus the inspected `PileUpAudit.png` and PDF.

The useful next diagnostic is to record the pre-AMR spectra and driving
components whenever p refinement is requested, especially at 15.8M and
15.9M, and test the proposed veto offline against those requests. The current
audit supports retaining the existing coarsening safeguards and keeping
pile-up detection passive until that decision-time evidence is available.

## Measuring accuracy and constraints after p-coarsening

`ChalmersCoarseningExperiment.py` directly tests the saved plateau states rather
than accepting or rejecting changes using the Chalmers spectral budget. It
examines all six innermost elements at every integer time from 11M to 20M:
120 candidates remove one point in one angular direction, and 60 remove one
point in both angular directions. The angular extents are 15 before 16M and
17 thereafter. The radial resolution and outer elements stay fixed.

Each candidate uses the native full-element p-projection matrix and then the
native exponential filter with alpha 36 and half-power 256 in all directions.
There is no reset of Pi or Phi to satisfy constraints. An additional control
reapplies the filter at the original resolution to separate filter roundoff
from the effect of changing resolution.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  build-agents/bin/python-spectre tools/Amr/ChalmersCoarseningExperiment.py \
  build-agents/chalmers-gh-feedback-long/tighten-20 \
  --references \
    build-agents/chalmers-gh-feedback-long/cartoon-graded4-20 \
    build-agents/chalmers-gh-feedback-long/cartoon-graded5-20 \
  --output build-agents/chalmers-gh-coarsening-measured \
  --times 11 12 13 14 15 16 17 18 19 20 --workers 2
```

State error is measured against both spherical numerical references on common
radial partitions. The three-index constraint is recomputed from the candidate
metric derivatives and Phi. The gauge constraint uses the candidate Pi/Phi
and a freshly computed damped-harmonic gauge source. These are coordinate-volume
L2 norms with Euclidean component contractions and symmetric multiplicities,
not invariant norms. Both element-local and combined inner-shell changes are
recorded. This experiment checks the gauge and three-index constraints, not
the complete GH constraint hierarchy or subsequent evolution.

<table>
<tr>
<th>
Candidate
</th>
<th>
Local state-error change
</th>
<th>
Local three-index constraint change
</th>
<th>
Local gauge-constraint change
</th>
</tr>
<tr>
<td>
One angular direction, 120 tests
</td>
<td>
-0.06571% to -0.000291%
</td>
<td>
+0.1204% to +3.0658%
</td>
<td>
-3.6491% to +5.3338%
</td>
</tr>
<tr>
<td>
Both angular directions, 60 tests
</td>
<td>
-0.13128% to -0.000580%
</td>
<td>
+0.2406% to +6.0453%
</td>
<td>
-6.6573% to +10.5730%
</td>
</tr>
</table>

Both references agree that the state error decreases for every candidate.
However, every candidate increases the three-index constraint. The gauge
constraint sometimes improves and sometimes worsens. Thus none passes the
proposed test of preserving accuracy and both measured constraints, even
though the actual changes are much more informative than the previous
69--108-times spectral-budget rejection.

Combining identical reductions across all six inner elements gives these
whole-domain changes. The other elements' fields and norm contributions are
unchanged; these are instantaneous measurements, not evolved results.

<table>
<tr>
<th>
Time
</th>
<th>
Reduction
</th>
<th>
Total grid points removed
</th>
<th>
Global three-index constraint change
</th>
<th>
Global gauge-constraint change
</th>
</tr>
<tr>
<td>
15M
</td>
<td>
One angular direction
</td>
<td>
900 (1.574%)
</td>
<td>
+0.3378%
</td>
<td>
-0.2410%
</td>
</tr>
<tr>
<td>
15M
</td>
<td>
Both angular directions
</td>
<td>
1,740 (3.043%)
</td>
<td>
+0.6728%
</td>
<td>
-0.4333%
</td>
</tr>
<tr>
<td>
20M
</td>
<td>
One angular direction
</td>
<td>
1,020 (1.695%)
</td>
<td>
+0.6720%
</td>
<td>
-0.0731%
</td>
</tr>
<tr>
<td>
20M
</td>
<td>
Both angular directions
</td>
<td>
1,980 (3.291%)
</td>
<td>
+1.3411%
</td>
<td>
-0.0812%
</td>
</tr>
</table>

At 20M the full-domain state error changes by only -0.00000517% for one
direction and -0.00001032% for both. The principal measured tradeoff is fewer
points with a modest, consistently positive three-index constraint increment.
No runtime saving was measured in this snapshot experiment. The initial
decision not to evolve these candidates used a strict no-increase condition.
That condition is **not SpEC's coarsening rule** and is superseded by the
tolerance-based experiment below. A small constraint increase alone is not
grounds for rejecting coarsening.

The two references change the reported inner-shell state-error change by at
most 1.42e-5 percentage points. Angular/radial error quadratures (20,40) and
(28,48) agree to 2.82e-15 in component norms. Constraint quadratures at orders
20 and 28 agree to 5.08e-10 absolute or better. Refiltering alone changes norms
by less than 3.4e-9 percent. Native projection/filtering agrees with independent
modal truncation/filtering on a random polynomial to 5e-14; the gauge-constraint
formula agrees with the separate unit-test reference to 5.6e-16. The recomputed
gauge source matches saved native GaugeH to 1.34e-14, and the original global
constraints/nodal gauge maxima reproduce the previous analysis/native output.

The study directory contains ten `Snapshot*.json` files, `GlobalControls.json`,
`Summary.json`, `Validation.json`, the aggregation and plotting scripts, and
the inspected `MeasuredCoarsening.png` and PDF. Independent projection and
differentiation generally do not commute, but a new coupled metric/Phi
projection is not a prerequisite for testing SpEC's established policy.

## SpEC-style tolerance and pile-up coarsening

`gh::ChalmersWithSpecCoarsening` is a separate, opt-in 3D experiment. It wraps
the existing GH refinement criterion and replaces its p-coarsening checks with
the two paths in SpEC's `ChangeSpectralGrid`. It also applies SpEC's permission
checks to proposed p increases, so a discarded mode is not immediately restored
solely because the two estimators use different normalizations. Chalmers still
requests refinement and chooses h versus p. This is not a port of SpEC's full
controller: neighbor smoothing, adaptive filters, and topology-change rules
are not implemented here. The optional stored local target adjustment is
described below; setting `TruncationErrorAdjustmentLimit: 0` retains the
original comparison policy in this section.

The native diagnostic uses two RMS power spectra, one over independent metric
components and one over the combined independent Pi/Phi components. Each
component first uses `PowerMonitors::power_monitors`. It therefore averages
over transverse modal indices, rather than using the Legendre-L2 directional
energy of the existing Chalmers estimator. The highest, effectively killed
mode is excluded. Relative truncation error uses a Gaussian-weighted geometric
mean near the last retained modes divided by the larger of the first two
coefficients. Repeating that estimate with one fewer unfiltered mode gives
the p-down estimate, clamped to be at least the current estimate.

The pile-up count uses SpEC's lowest-third fit and the additional fits starting
at each of the first three samples. Fits are weighted by their slope errors;
local windows of up to five samples contribute `exp(-32*(local/main)^2)`.
The implementation retains the original mode indices across parity gaps.
Samples below `100*epsilon*max_power` are excluded to avoid treating parity
zeros and roundoff as physical spectra. A flat or growing whole spectrum, a
short spectrum, or a tail entirely at roundoff cannot authorize the pile-up
override. These finite-data safeguards are explicit adaptations, not claims
of bitwise agreement with legacy SpEC behavior on degenerate spectra.

Both two-index and three-index GH constraints are recomputed from the current
fields, including fresh gauge derivatives. The derivative index is projected
along normalized element-map Jacobian columns. The nodal RMS then sums
independent tensor components once, using SpECTRE's existing
`Constraints_detail::logical_constraints`. This differs from the earlier
feedback criterion's symmetric multiplicities. It also uses the RMS directly,
without the extra square root in the legacy SpEC topological aggregation.
Consequently these tolerances must be specified in the documented convention;
SpEC input-file numbers cannot simply be copied.

The decisions are:

- Ordinary p-down requires `10^0.1 * E_next <= spectral_target` and
  `1.25 * C_predicted_next <= constraint_target` for both constraints.
  Current and next spectral errors must also be below the target. The
  derivative spectrum predicts C_ia and the metric spectrum predicts C_iab.
- `C_predicted = C_current^2/C_previous`, with SpEC's `2*C_current` fallback
  when the old value is tiny. Multiply by `max(1,E_next/E_current)` to estimate
  the effect of removing a mode. No monotonic-constraint condition is imposed.
- With `ForceDropPileUpModes`, a minimum pile-up count above 3 across the
  nondegenerate monitored tensor groups can override both ordinary budgets.
  A drop still requires `pile_up > 4.1` or `E_next < spectral_target`, avoiding
  immediate reversal near the threshold. Minimum resolution still applies.
- Proposed p increases require an excessive spectral error or predicted
  constraint and pile-up at most 3, or 1.5 when the forced policy is enabled.
  This wrapper does not independently request refinement when Chalmers has
  made no refinement request.
- Resolution changes have a two-AMR-call cooldown, measured from the realized
  mesh change. The numerical coarsening helper supports a minimum extent of
  five, but this wrapper requires at least seven points: the Chalmers parity
  fit needs three samples per parity plus the filtered mode. The first
  evaluation does not coarsen. History is per element, serialized on
  migration, preserved across
  p changes, and reset on h changes. Constraint forecasts never cross a mesh
  change. Observing a criterion does not advance its history.
- Chalmers h-splits retain priority. Joins retain the existing sibling and
  prospective-parent checks and additionally require both current constraint
  norms to meet their budgets in every direction. An h-join does not acquire
  additional p-coarsening on another axis, since the prospective-parent check
  assumed the existing polynomial extents.

`ChalmersSpecCoarsening.py` runs the comparison and writes native pre-AMR
diagnostics to `SpecDecisions.csv` and `SpecAnalysis.json`. Its snapshot audit
uses the same C++ diagnostic through Python bindings, including the prescribed
damped-harmonic gauge and gamma2 profile:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  build-agents/bin/python-spectre tools/Amr/ChalmersSpecCoarsening.py \
  --audit-source build-agents/chalmers-gh-feedback-long/tighten-20 \
  --output build-agents/chalmers-gh-spec-coarsening/calibration
```

Across the 10M, 15M, and 20M snapshots (114 element observations), both tensor
groups exceed the pile-up threshold in 24 directional observations. The maximum
minimum count is 5.243. Maximum projected RMS constraints are 0.001430 for C_ia
and 0.003042 for C_iab. The initial evolution comparison uses rounded-up
budgets `[0.002, 0.004]` and a relative spectral target of `1.e-4`. These are
calibrated experiment parameters, not universal production tolerances.

The two runs start from the same Kerr-Schild data and retain the established
refinement/constraint-feedback policy until 10M. From 10M onward they enable
ordinary coarsening or ordinary plus forced pile-up coarsening. Each uses five
Charm++ PEs, for ten total. The original tightening trajectory and the two
converged spherical references provide controls. The field projection and
steep alpha-36, half-power-256 filter are unchanged.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  build-agents/bin/python-spectre tools/Amr/ChalmersSpecCoarsening.py \
  --executable \
    build-agents/chalmers-gh-spec-coarsening/EvolveGhSingleBlackHole-run \
  --output build-agents/chalmers-gh-spec-coarsening/tolerance-20 \
  --policy tolerance --start-time 10 --end-time 20 \
  --constraint-targets 0.002 0.004 --processes 5
```

Use `--policy pileup` and a distinct output directory for the forced variant.
The initial ordinary attempt exited in Charm++ startup before initialization;
the unchanged retry is stored in `tolerance-20-retry`. The manifest preserves
the executable hash, parameter choices, controls, and source snapshot.

The focused C++ tests check exponential and plateau spectra, parity zeros,
filtered modes, manufactured two-/three-index constraints, tolerance-based
coarsening with rising constraints, forecast rejection, forced overrides,
both-tensor agreement, cooldown, minimum modes, and serialized history.
The four GH tests and eight affected AMR action/event/projector tests pass.
The 0.3M two-PE execution smoke test also completes successfully; it exercises
the diagnostic and history path but makes no new coarsening decisions.

### Completed 20M comparison

Both five-PE runs complete successfully. Their 21,294 native directional
diagnostic/decision records are identical, as are the fields at all 26 saved
times. The maximum minimum pile-up count is 1.000805, below either pile-up
threshold, so the forced override is never exercised. Ordinary tolerance-based
coarsening and its p-addition permission checks suffice in this evolution.

<table>
<tr>
<th>
Quantity
</th>
<th>
Original Chalmers + feedback
</th>
<th>
Ordinary SpEC-style policy
</th>
<th>
With forced pile-up enabled
</th>
</tr>
<tr>
<td>
Final elements
</td>
<td>
42
</td>
<td>
42
</td>
<td>
42
</td>
</tr>
<tr>
<td>
Final grid points
</td>
<td>
60,162
</td>
<td>
48,870
</td>
<td>
48,870
</td>
</tr>
<tr>
<td>
Grid-point steps, entire run
</td>
<td>
221,367,312
</td>
<td>
194,312,748
</td>
<td>
194,312,748
</td>
</tr>
<tr>
<td>
Grid-point steps, 10M onward
</td>
<td>
144,817,446
</td>
<td>
117,751,590
</td>
<td>
117,751,590
</td>
</tr>
<tr>
<td>
Final state L2 difference from fine reference
</td>
<td>
0.0653114135
</td>
<td>
0.0653114078
</td>
<td>
0.0653114078
</td>
</tr>
<tr>
<td>
Final three-index constraint L2
</td>
<td>
0.00464448
</td>
<td>
0.00367892
</td>
<td>
0.00367892
</td>
</tr>
<tr>
<td>
Final gauge constraint L2
</td>
<td>
0.000344017
</td>
<td>
0.000331773
</td>
<td>
0.000331773
</td>
</tr>
<tr>
<td>
Realized p-increases / p-decreases
</td>
<td>
90 / 18
</td>
<td>
42 / 48
</td>
<td>
42 / 48
</td>
</tr>
<tr>
<td>
Realized h-splits / h-joins
</td>
<td>
18 / 0
</td>
<td>
18 / 0
</td>
<td>
18 / 0
</td>
</tr>
<tr>
<td>
Immediate p/h reversals
</td>
<td>
0 / 0
</td>
<td>
0 / 0
</td>
<td>
0 / 0
</td>
</tr>
</table>

The new policy saves 18.77% of the final points and 18.69% of grid-point steps
after activation (12.22% over the whole run). This is a work proxy, not a
wall-time speedup: the old and new runs use different PE counts and concurrent
loads. The final three-index and gauge constraint norms decrease by 20.79% and
3.56%, respectively. Peak three-index constraints are unchanged because that
peak precedes activation. Excision remains outflow, and the maximum horizon
mass error remains about 1.64e-6.

The final solution-error change is only -8.85e-6 percent. Switching between the
two converged spherical references changes the measured error by 2.56e-7 in
absolute units, much more than the 5.78e-9 difference between these policies
and the control. Thus the supported conclusion is unchanged solution accuracy,
not a measured accuracy improvement. Independent constraint integrals at
Gauss orders 20 and 28 agree to better than 7e-16 at 15M and 20M for the new
runs. All 15 saved observations before 10M reproduce the original fields and
coordinates bit-for-bit.

The first new action removes one point in each angular direction from all six
inner elements at 10M (extents 13 to 12); later radial coarsening reduces other
elements.
The controller also suppresses repeated p-addition requests from the tighter
componentwise estimator. These results test that combined policy, rather than
isolating point removal from its corresponding p-addition permission rule.
They support using ordinary tolerance-based coarsening in the next GH test;
they do not establish an additional benefit from forced pile-up removal.

The build directory `build-agents/chalmers-gh-spec-coarsening` contains
`Summary.json`, `CrossPolicyComparison.json`, the per-run validation and
reference comparisons, and `SpecCoarseningComparison.png` / `.pdf`.
`analyze_study.py` and `plot_results.py` preserve the analysis and plotting
commands. The post-launch h-join and minimum-seven-point guards are checked
against both complete decision histories: neither run requested a join,
used fewer than seven points, or requested a p decrease from seven points.
Consequently those guards do not alter the benchmark trajectories.
`SourceSnapshot`, `FinalSourceSnapshot`, `FinalSourceHashes.json`, and the
archived benchmark GH library document the source versions used.
The final guarded source builds successfully, passes all four focused GH tests
and the input parse, and completes a fresh 0.3M two-PE smoke run with 144 native
directional evaluations. Together with the eight affected AMR framework tests,
these checks are recorded in `Verification.json`.

## Tilted Kerr gauge relaxation

`ChalmersTiltedKerr.py` extends the matched GH experiment to a mass-one Kerr
hole with dimensionless spin `0.5 * (1,2,3) / sqrt(14)`. The hole has no boost
and the cubed-sphere maps remain static. The initial spatial metric is
Kerr-Schild; GH initializes Pi and Phi consistently with damped harmonic gauge
and the discrete three-index constraint. Damping, filtering, excision, the
physical Bjorhus outer boundary, and the 0.1M adaptation schedule are inherited
from the Schwarzschild test.

`--policy original` records native SpEC diagnostics while preserving the
existing Chalmers/constraint-feedback decisions. `--policy ordinary` enables
the tolerance-based coarsening and corresponding p-addition permission rules
at 10M. Both retain the existing relative spectral target `1e-4` and projected
constraint targets `[0.002, 0.004]`; these are transferred experimental
tolerances, not newly calibrated Kerr production settings. Forced pile-up
removal is disabled. `--policy fixed` disables AMR for full 3D reference runs.
`--radial-levels` sets the refinement level independently in the four existing
radial shells, without changing their block maps. `--angular-points` and
`--radial-points` likewise accept four extents, ordered from the inner shell
outward. This permits fixed reference meshes that resolve the tilted hole
strongly near excision without using the same angular degree in the far field.

The analysis now checks the initial metric against the native Kerr-Schild
solution, tracks spin magnitude and direction as well as mass, and uses
`ChalmersGaugeRelaxationAnalysis.compare` to compare full 3D fields on the
common element partition. A stationary Kerr-Schild field comparison is not an
evolved accuracy estimate once the gauge has relaxed. The spherical Cartoon
reference is not applicable to this experiment.

The study is stored in `build-agents/chalmers-gh-tilted-kerr`. The initial
5M pilot uses the existing checked executable; the longer matched evolutions
use an optimized `release-debug` build with assertions and NaN initialization
enabled. A short overlap check against the debug pilot validates that build.
Runs never overwrite existing directories.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  build-agents/bin/python-spectre tools/Amr/ChalmersTiltedKerr.py \
  --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/chalmers-gh-tilted-kerr/pilot-5 \
  --policy original --end-time 5 --processes 10 --horizon-l 6
```

The 5M pilot completes with positive excision characteristic speeds (minimum
0.1998), 70 p-increases, no h changes, and no immediate mesh reversals. It grows
from 27,648 to 37,236 points. The independently integrated three-index
constraint
peaks at 0.01022 and ends at 0.003503. These are coordinate-volume L2 norms;
the native projected nodal RMS tolerances use different normalization.

The pilot also exposes a horizon-observer resolution effect: its measured spin
jumps when the horizon finder reduces its angular resolution from L=7 to L=6.
Matched 0.3M runs with horizon floors L=6 and L=10 have bitwise identical
evolved
fields at all four saved times, while the maximum spin-magnitude error falls
from 4.309e-4 to 2.811e-6 and the mass error from 6.663e-5 to 8.072e-7.
The driver therefore defaults to `--horizon-l 10`. Horizon feedback controls
remain disabled. `OptimizedValidation.json` also records bitwise agreement
between the checked optimized and pilot executables at 0, 0.1, and 0.2M;
the 0.3M terminal checkpoint is excluded from that overlap check because a
terminal run does not perform the next AMR phase.

The planned fixed 20M references use radial levels `[0,0,1,2]`. The medium
reference has angular extents `[18,16,14,12]` and radial extent 10; the fine
reference has angular extents `[20,18,16,14]` and radial extent 12. They have
92,880 and 145,440 points respectively. A third evolution repeats the fine
mesh with CFL 0.25 instead of 0.5. Comparisons at 5, 10, 15, and 20M integrate
all GH evolved fields over their common three-dimensional element partition.
Reference convergence must be checked from these evolved fields, not inferred
from grid size alone.

A uniform-angular 307,200-point reference was stopped at 0.499M after a runtime
assessment; its incomplete output is retained in `reference-fine-20` and is
not a completed accuracy reference. A separate standard-release build passed
the four focused tests and reproduced all saved fields in both an adaptive
0.3M test and a fixed-grid 0.3M test exactly. Its modest speed advantage did
not justify changing the checked production build. The recorded manifests,
`ReleaseValidation.json`, and `BuildSetupNotes.txt` give the build details.

The original-policy 20M attempt stops at 11.8M: blocks B4 and B5 request a
p-increase in direction 1 while already at 21 points, reaching the supported
degree-20 limit. The observed minimum tensor pile-up in those directions is
4.568. This is a refinement-policy error, not a completed evolution or evidence
of an excision failure. `ControlFailure.json` preserves the triggering
diagnostics. A p-to-h fallback at the resolution cap is not implemented.

The run also prompted an input-bound correction in both coarsening drivers:
`NumPolynomialModes` is converted to an extent by adding one, so its lower
bound must be 6 to agree with the wrapper's seven-point floor. The failed
control never requested a decrease from eight or fewer points, and is
unaffected by this correction. The early `ordinary-20` attempt is retained as
an interrupted setup trial; the corrected run is `ordinary-min7-20`.

Failed evolutions can be inspected explicitly with
`ChalmersGaugeRelaxationAnalysis.py analyze RUN --allow-incomplete`. This checks
that saved observations form a complete prefix, retains the failed return
code, records the separate final field/mesh/time-step observation times, and
writes `PartialAnalysis.json` and `PartialHistory.csv`. Partial field
comparisons require both `--allow-incomplete` and an explicit `--times` list.
The default analysis continues to reject failed runs.

The corrected `ordinary-min7-20` attempt passes the control's 11.8M failure:
the SpEC permission check vetoes the two p-additions because both tensor
spectra have pile-up. It reduces the mesh from 45,454 to 43,890 points before
angular h-refinement starts at 13.1M. At 13.2M it encounters a macOS worker
stack overflow in `db::copy_items` while sending parent data to children.
This is a parallel AMR infrastructure failure, not evidence of a numerical
instability. The saved field prefix ends at 13M.

A forced angular-split reproducer fails on four PEs with the old executable
and succeeds on one PE. Disassembly shows a 185,104-byte frame for the
old DataBox copy routine; the worker stack is 536,576 bytes, including its
callers. `DataBox::copy_items` now deserializes each copied item directly into
the result tuple rather than keeping all constructor-argument temporaries
alive together. The checked optimized build passes all 34 selected DataBox
and AMR unit tests. The four-PE split reproducer then succeeds, with bitwise
identical fields and coordinates to the old one-PE run at all four saved
times through 0.3M. An adaptive smoke test also gives bitwise identical
fields to the earlier executable. `StackFixValidation.json` and
`StackTestsFinal.log` record these checks. The full rerun is isolated in
`ordinary-stackfix-20`; the failed output remains preserved.

The stack-fixed rerun completes the full 20M evolution on four PEs. It ends
with 38 elements and 99,218 grid points. The native criterion records 176
p-increase, 26 p-decrease, and 14 split requests; no joins are requested.
The final apparent-horizon mass is 1.00000132351 and dimensionless spin
0.499998220095. These observables remain close to their expected values,
but they do not replace checks of the evolved fields and constraints.

Independent integration uses Gauss orders 28 and 36. The earlier 20/28 pair
is insufficient for some degree-20 fields in this experiment. Spot checks
with orders 20, 28, 36, and 44 confirm convergence of the higher pair;
`QuadratureChecks.json` records the tested elements and times. Reported norms
must be regenerated at these orders rather than mixing older analyses.

The fine fixed mesh completes 20M, but its three-index constraint grows
from 6.21e-4 at 10M to 0.312 at 20M, mostly in the innermost shell. Completion
alone therefore does not validate it as an accuracy reference. The medium
mesh and fine half-CFL evolution test spatial and time-step sensitivity.
Filtering is applied every substep, so changing CFL also changes how often
the filter is applied per unit physical time. The time-step comparison
measures the complete filtered evolution algorithm, not only the ODE solver.

### Completed tilted-hole study

The adaptive SpEC-policy run, fine fixed grid, and fine half-CFL grid all
complete 20M. The medium fixed grid stops near 19.78M when the excision
outflow check detects an ingoing characteristic speed of about -1.58e-5 in
inner blocks B4/B5. Its last recorded step is 19.7809608M and its last full
field snapshot is 19M. `MediumFailure.json` records the failure. Its saved
prefix is useful for comparisons through 19M, but is not a completed 20M
reference. The original policy's separate 11.8M degree-limit failure also
prevents a full-length original-versus-coarsening comparison.

At the last common integer time before the original run stops:

| Quantity at 11M | Original policy | Ordinary SpEC policy |
|---|---:|---:|
| Grid points | 45,916 | 43,890 |
| Grid-point steps over 10–11M | 18,181,600 | 17,692,480 |
| Three-index constraint L2 | 0.00446624 | 0.00460408 |
| Gauge constraint L2 | 0.000397897 | 0.000440733 |
| Full-state difference from fine grid | 0.0157897575 | 0.0157898002 |
| Inner-shell state difference from fine grid | 0.000342699 | 0.000344664 |

The policy saves 4.41% of the points and 2.69% of the grid-point-step work
proxy over this short common interval. The three-index and gauge constraints
increase by 3.09% and 10.77%, respectively; this is a measured tradeoff, not
unchanged constraints. The full-state difference changes by 4.28e-8, while the
difference between fixed resolutions is 0.002835. Restricting to `1.2<r<2M`,
the policy changes the measured difference by 1.97e-6 (0.57%), smaller than
the fixed-resolution difference of 1.367e-4. These comparisons do not resolve
a change in solution accuracy at that level. All 15 saved observations before
10M reproduce the original fields and coordinates bit-for-bit.

The final adaptive mesh has 176 realized p-increases, 26 p-decreases, 14
h-splits, and no joins or immediate p/h reversals. The full log contains
16,224 unique directional records covering every element and all three
directions at all 199 evaluations. No changes occur before activation and no
p-decrease violates the resolution floor. The stack-fixed rerun reproduces
all 19 saved field/coordinate snapshots of the crashed attempt, through 13M,
bit-for-bit. `Validation.json` and `StackfixLongOverlap.json` retain these
checks.

At 20M the adaptive three-index constraint is 0.0695840, after leveling off
between 18M and 20M. The fine fixed grid reaches 0.311634; halving CFL changes
that value to 0.311642. The corresponding adaptive and fine gauge constraints
are 0.00108477 and 0.0170376. The adaptive excision boundary remains outflow,
with minimum sampled outgoing light speed 0.199779. Maximum horizon mass and
spin-magnitude errors are 1.61e-6 and 2.81e-6. Good horizon observables alone
do not establish accuracy of all evolved fields.

Radial integration localizes 99.997% of the adaptive inner-shell squared
three-index norm, and 99.918% of the fine-grid value, to `1.2<r<1.4M`.
The horizon radius is about 1.93–1.99M at 20M, so this is deep inside the
horizon, near excision. The adaptive constraint norm outside `r=2M` is
0.00327821; the fine-grid value there is 0.000117392. The smaller global
adaptive norm therefore does not imply smaller error in every region.
`RadialConstraintProfile.json` records converged integrals in five radial
bands, including an independent check of each band's geometric volume.

The late inner-shell numerical reference is not well established: at 19M,
the medium/fine state difference in that shell is 0.02105, larger than the
adaptive/fine difference of 0.007337. Outside `r=2M`, those differences are
0.001423 and 0.01363. Halving the fine-grid CFL changes the full state by
5.69e-7 at 20M and its exterior part by only 1.32e-9. Thus time-step error is
small compared with the measured spatial differences; it does not explain
the inner-shell growth.

All three-index and field-difference integrals pass the 28/36-point
quadrature check. The nonlinear gauge-constraint integral adaptively raises
its order through 44 points for late fine-grid snapshots, retaining the
relative convergence tolerance of 2e-5. Recomputed damped-harmonic gauge
sources agree with saved native values to about 3e-14. The initial failed
28/36 gauge checks and successful higher-order retries are retained in the
study logs.

This is a completed single-hole experiment, not a BBH accuracy validation.
A useful next comparison is to enable SpEC's forced pile-up removal on the
same tilted-hole problem: the smaller of the two monitored tensor pile-up
counts reaches 6.21 on the ordinary trajectory and 8.00 in the fixed fine
snapshots.
Unlike the Schwarzschild test, this case can exercise that override.
A better resolved inner-shell reference and the original criterion's
p-limit-to-h fallback remain separate follow-up work.

`Summary.json`, `TiltedKerrResults.png` / `.pdf`, and
`TiltedKerrSpatialChecks.png` / `.pdf` contain the collected results.
`aggregate_results.py`, `plot_results.py`, `plot_spatial_checks.py`, and
`radial_constraints.py` preserve the analysis. The executable/source
manifests, the 34-test DataBox/AMR verification, and all failed or interrupted
runs are retained in the study directory. Wall times are not used to claim a
speedup because PE counts, concurrent workloads, and deliberate pauses differ.

### Tilted Kerr with forced pile-up removal

The matched forced-pile-up run completes 20M on ten PEs. The driver now accepts
`--policy pileup`. The only difference between its production input and the
ordinary control is `ForceDropPileUpModes: true`; the executable and linked
libraries are unchanged. Activation remains at 10M. All 15 saved snapshots
before activation reproduce the ordinary fields and coordinates bit-for-bit,
despite using ten rather than four PEs. The field comparison at 10M is also
identical.

Run from the repository root:

```sh
build-agents/bin/python-spectre tools/Amr/ChalmersTiltedKerr.py \
  --executable \
    build-agents/chalmers-gh-tilted-kerr/EvolveGhSingleBlackHole-stackfix \
  --output build-agents/chalmers-gh-tilted-kerr/pileup-20 \
  --policy pileup --start-time 10 --end-time 20 --processes 10
```

The saved output already occupies this directory; choose a new directory to
repeat the experiment. Results below compare `pileup-20` with
`ordinary-stackfix-20` in `build-agents/chalmers-gh-tilted-kerr`.

| Quantity | Ordinary | Forced pile-up |
|---|---:|---:|
| Final grid points | 99,218 | 52,716 |
| Final elements | 38 | 30 |
| Grid-point steps, 10–20M | 676,744,488 | 155,224,468 |
| Global steps, 10–20M | 8,554 | 3,408 |
| Realized p-decreases, full run | 26 | 54 |
| Realized p-increases, full run | 176 | 160 |
| h-splits / h-joins | 14 / 0 | 6 / 0 |
| Three-index constraint L2 at 20M | 0.0695840 | 0.00749437 |
| Gauge constraint L2 at 20M | 0.00108477 | 0.000307604 |
| Three-index L2 outside r=2M at 20M | 0.003278206 | 0.003278205 |
| Maximum horizon mass error | 1.61e-6 | 1.70e-6 |
| Maximum spin-magnitude error | 2.81e-6 | 2.81e-6 |

Final points fall by 46.87%, and grid-point-step work after activation falls
by 77.06%. The first h-split moves from 13.1M to 17.2M, and the forced run
also takes fewer global steps. Work means the sum of grid points for steps
starting in `[10M,20M)`; it excludes diagnostics and AMR overhead. Differing
PE counts make elapsed-time speedup an inappropriate comparison.

The final three-index and gauge norms decrease by 89.23% and 71.64%.
The improvement is concentrated near excision: the three-index norm in
`1.2<r<1.4M` falls from 0.0695057 to 0.00673273. The exterior constraint
is essentially unchanged. Direct integration over exterior blocks gives a
forced-versus-ordinary full-state L2 difference of 3.435e-7 at 20M; this
check avoids subtracting two nearly equal whole-domain and inner-shell
squared norms. The inner-shell state difference is 0.00124270. Both runs
retain outflow at excision, with minimum sampled outgoing light speed
0.199779.

The late inner-shell reference remains unresolved. At 20M the forced/fine
inner-shell state difference is 0.00826358, versus 0.00834838 for ordinary/fine.
Outside r=2M the corresponding differences from the fixed fine grid are both
about 0.0138346. These results demonstrate lower interior constraints and
lower work without a resolved exterior accuracy penalty; they do not
establish convergence of the interior solution or validate BBH evolution.

The native decision audit reproduces all 14,796 directional records over
199 evaluations, with no resolution-floor or cooldown violation. Of the
50 p-down requests after activation, 28 bypass ordinary error budgets.
Two exceed the spectral budget, 18 the predicted three-index budget, and
16 the predicted two-index budget; these categories overlap. The largest
ratios to the ordinary budgets are 1.050 for the spectral estimate and
1.882 for the constraint estimate. Replaying ordinary decisions on the
same forced trajectory identifies those 28 additional removals and 52
additional p-addition vetoes. The flag changes both forced removal and
p-addition hysteresis, so the whole-run improvement is a result of that
combined policy.

There is increased cycling despite compliance with the cooldown. After
activation, the forced run has 24 realized remove-then-add pairs on the
same element/direction, versus six for ordinary. Twelve forced-run pairs
are separated by 0.2M, the earliest permitted return; the ordinary minimum
is 1.2M. Twenty of the 24 pairs start with a forced budget override, and
all additions are requested by the criterion rather than neighbor balancing.
Neither run has a reversal on consecutive 0.1M checks. Thus zero immediate
reversals does not mean zero short cycles. A useful next controlled test is
to reduce these short cycles while preserving the demonstrated savings.

All field and three-index integrals pass the 28/36-point quadrature check;
the forced-run gauge integrals pass the same pair at relative tolerance
2e-5. Recomputed gauge sources agree with saved values within 2e-14.
`PileupSummary.json`, `TiltedKerrPileup.png` / `.pdf`,
`PileupExteriorVerification.json`, and `PileupCycleAudit.json` retain the
results. `analyze_pileup.py`, `summarize_pileup.py`, `plot_pileup.py`,
`verify_pileup_exterior.py`, and `audit_pileup_cycles.py` preserve the analysis;
`PileupManifest.json` records the executable, library hashes, and input check.

### How SpEC limits repeated removal and addition

The follow-up audit reads the local SpEC `ChangeSpectralGrid.cpp` and the
shipped BBH inspiral/ringdown inputs. `SpecControllerAudit.json` records the
exact revision and source hashes. SpEC combines several controls:

- A minimum number of calls between extent changes. The BBH inputs use two,
  which the experimental wrapper already has.
- An additional physical-time restriction in the inspiral input. Its setup
  chooses `min(1, 2*min(M_A,M_B))` using initial horizon masses, and multiplies
  this by 0.01 in frequent-trigger mode. Non-quasicircular orbits select that
  mode automatically. The shipped ringdown input does not set this option;
  the class default is zero. This is not a universal one-M wait in SpEC.
- Pile-up hysteresis: forced removal above three modes and permission to add
  only below or at 1.5 modes. Ordinary removal also has the existing spectral
  factor `10^0.1` and constraint factor `1.25` margins.
- A stored, constraint-driven shift of the logarithmic spectral target.
  For large constraints it tightens by increments based on the convergence
  slope (typically half a spectral mode, with an alternative 0.8-mode
  placement relative to the measured error). Recovery requires decreasing
  constraints and both current and predicted constraints below half their
  tolerance; it then relaxes by only 0.1 mode. Inspiral bounds the log10 shift
  to `[-1,0]`. Its quasi-local step also limits target differences between
  related neighboring directions to a factor of two.
- Consistent use of the shifted spectral target for refinement and the
  prospective pile-up-removal check. In particular, when the pile-up count
  is at most 4.1, the estimated error after removing a mode must meet the
  shifted target. The ordinary permission budgets remain a separate gate.
- Adaptive filtering in the inspiral input, with separate activation and
  deactivation thresholds. Ringdown disables adaptive filtering. The current
  SpECTRE experiment retains its fixed steep filter.

Relevant SpEC functions are `SetAddDropLimits`,
`AdjustTruncErrorShift_Local`, `AdjustTruncErrorShift_QuasiLocal`,
`ComputeGridChangesWithOrWithoutBuffer`, and `FilterChangeRequests`.
The time restriction is assigned by `InputFiles/Bbh/DoMultipleRuns.input`
and enforced separately from the call counter. These mechanisms reduce
cycling; they do not constitute a guarantee that SpEC cannot cycle.

All twelve rapid returns in the two-call tilted-Kerr run are driven by the
three-index constraint. On re-addition the spectral estimate is only
0.412–0.773 times the unshifted target, while the projected three-index
constraint exceeds its target. The derivative-spectrum pile-up has fallen
below 1.5, although the metric spectrum can still have substantial pile-up.
The wrapper therefore permits Chalmers to restore a mode that the forced
drop just removed. For example, B4 direction 0 removes at 15M with
`E_next/T=0.625`, then adds at 15.2M with `E/T=0.634` and `C3/C3_target=1.156`.

`replay_spec_target_shift.py` applies SpEC's local target-update equations
to the saved native diagnostics, using this experiment's constraint
normalization and budgets. It initializes the shift to zero at 10M, retains
previous constraints across p changes as SpEC does, and resets on new h
element IDs. It omits neighbor smoothing and adaptive filtering. This is a
diagnostic replay on unchanged fields, not a SpEC evolution or a new
SpECTRE trajectory.

The replay reaches shift -1 for all twelve short-cycle removals; their next
errors are 3.81–7.33 times the shifted target. The consistent shifted-target
guard would reject all twelve of those proposals. It also rejects 22 of the
30 observed pile-up drops in total, so retaining the previous savings cannot
be assumed. `SpecTargetShiftReplay.json` stores every update and the selected
cycles. A live test of this controller must apply the same target to both
refinement and coarsening, rather than adding only a veto on coarsening.

The simpler timing experiment uses `--calls-between-changes 10`, with the
same forced policy and 0.1M AMR cadence. This emulates a 1M interval between
p changes for this fixed-cadence test. It is not a general physical-time
implementation: h-splits still retain priority and h changes reset the
per-element history. The default remains two calls. The matched run is
`pileup-cooldown10-20`; `CooldownManifest.json` verifies that
`CallsBetweenChanges: 2 -> 10` is its only input change from `pileup-20`.

The ten-call run completes 20M on ten PEs with the same frozen executable
and libraries. Increasing the wait reduces the number of remove/add pairs,
but does not resolve the competing refinement and coarsening decisions:

| Quantity | Two calls (0.2M) | Ten calls (1M) |
|---|---:|---:|
| Remove/add pairs after 10M | 24 | 14 |
| Pairs at the earliest permitted return | 12 | 8 |
| Final grid points | 52,716 | 52,816 |
| Grid-point steps, 10–20M | 155,224,468 | 161,137,378 |
| Global steps, 10–20M | 3,408 | 3,540 |
| Realized p-decreases, full run | 54 | 50 |
| Realized p-increases, full run | 160 | 154 |
| h-splits / h-joins | 6 / 0 | 6 / 0 |
| Three-index constraint L2 at 20M | 0.00749437 | 0.0202548 |
| Gauge constraint L2 at 20M | 0.000307604 | 0.00201117 |
| Three-index L2 outside r=2M at 20M | 0.0032782053 | 0.0032782051 |
| Maximum horizon mass error | 1.6965e-6 | 1.7084e-6 |

Work increases by 3.81%, while the final three-index and gauge norms grow
by factors of 2.70 and 6.54. Most of the additional three-index violation
is in `1.2<r<1.4M`, where the norm is 0.0199845. The direct exterior state
difference between the two runs at 20M is only 3.614e-7, the maximum spin
error remains 2.8114e-6, and excision retains outflow. The final inner-shell
difference from the fine reference changes from 0.00826358 to 0.00828505;
the unresolved late interior convergence issue still applies.

The larger wait also permits stronger pile-up to develop. At 12.8M and
19.5M, four direction-wise drops use the `pile_up > 4.1` branch, which
bypasses the next-spectral-error check. The maximum bypassed ordinary
spectral and predicted-constraint budget ratios are 18.37 and 39.49. These
observations do not isolate the cause of the later constraint increase,
but reinforce that a longer cooldown alone is not a satisfactory remedy.
Keep the default at two calls. The next controller experiment should use
SpEC's stored local target shift consistently for refinement and removal,
including its slow relaxation when constraints are small and decreasing.
This preserves the possibility of coarsening after recovery. Test that
local mechanism before adding neighbor smoothing or adaptive filtering.

All 14,796 native directional decisions are reproduced by the corrected
Python replay, with no cooldown violation; the minimum realized p-change
gap after activation is 1M. All fifteen pre-10M volume snapshots match the
ordinary control bitwise, and the 10M field difference from the two-call
run is zero. The field, constraint, radial-profile, and direct exterior
integrals pass the independent 28/36-point quadrature comparisons. The
first analysis attempt exposed a Python variable-shadowing error in the
new configurable-cooldown replay; correcting it required no evolution
changes or rerun. The initial failure and successful recovery logs are
retained separately.

`CooldownSummary.json`, `CooldownCycleAudit.json`,
`CooldownExteriorVerification.json`, and `TiltedKerrCooldown.png` / `.pdf`
store the completed results in `build-agents/chalmers-gh-tilted-kerr`.
The `analyze_cooldown.py`, `audit_cooldown_cycles.py`,
`verify_cooldown_exterior.py`, `summarize_cooldown.py`, and
`plot_cooldown.py` scripts preserve the analysis. To repeat the evolution,
choose a fresh output directory:

```sh
build-agents/bin/python-spectre tools/Amr/ChalmersTiltedKerr.py \
  --executable \
    build-agents/chalmers-gh-tilted-kerr/EvolveGhSingleBlackHole-stackfix \
  --output build-agents/chalmers-gh-tilted-kerr/pileup-cooldown10-repeat \
  --policy pileup --calls-between-changes 10 --end-time 20 --processes 10
```

### Stored local target adjustment

`ChalmersWithSpecCoarsening` now accepts
`TruncationErrorAdjustmentLimit`. Set it to `-1` for SpEC's inspiral bound,
or `0` to retain the previous instantaneous-feedback/coarsening experiment.
The tilted-Kerr driver exposes this as `--target-shift-limit -1`; its default
is zero so existing comparison commands retain their meaning. Older YAML
inputs used with the new executable need the new key explicitly set to zero.

Each element stores one logarithmic target shift per logical direction.
At each active AMR evaluation, the local controller uses the minimum usable
convergence slope and maximum relative tensor error. Its two constraint
forecasts use the previous evaluation even across a p change, following
SpEC's target-update rule. Excessive constraints tighten the target by the
half-mode / 0.8-mode placement rules described above. A constraint permits
relaxation only when it is decreasing and both its current and predicted
values are below half its tolerance. Relaxation is one tenth of a mode per
evaluation. The most restrictive constraint update wins, and the shift is
clamped to `[TruncationErrorAdjustmentLimit, 0]`. A missing or nondecaying
spectral fit leaves the previous shift unchanged.

The same factor `10^shift` scales Chalmers' directional componentwise field
budgets and the wrapper's relative tensor target. These remain different
error estimators and normalizations; the stored factor is shared. In this
mode the stored feedback replaces the nested criterion's instantaneous
`ConstraintTarget`, `MinimumTargetFactor`, and `EnableCoarseningVeto`
behavior. Its field settings still select h versus p and evaluate candidate
sibling joins. P-addition requires error above the shifted tensor target and
the existing permission based on the base target or excessive constraints;
a large constraint alone cannot bypass the shifted spectral check. In this
mode its constraint forecast retains history across p changes too. Ordinary
removal and the prospective pile-up check use the same shifted target. As in
SpEC,
the separate ordinary permission margins use the base target, and pile-up
above 4.1 retains its exceptional removal path. The wrapper continues to
require a Chalmers refinement request before adding points.

The shift updates during p cooldowns. It is serialized with the per-element
AMR history, survives p changes, and resets for h splits/joins. Read-only
observations use temporary history, and repeated evaluation at the same
time and mesh does not advance the controller. Before `StartTime`, or with
`EnableCoarsening: false`, the shift remains zero and the earlier refinement
policy applies. The prospective-removal constraint forecast still avoids
comparing different meshes; this differs deliberately from the persistent
constraint history used to adjust the target.

Verbose diagnostics now append `target_shift` and `effective_target` to each
directional record. `ChalmersSpecCoarsening.py` preserves both in
`SpecDecisions.csv` and remains able to read older logs. The two-call cooldown,
pile-up thresholds, positive constraint damping, and fixed steep filter are
unchanged.

Validation for the implementation is saved under
`build-agents/chalmers-gh-tilted-kerr/target-shift-implementation`. The checked
optimized `release-debug` build passes all seven focused GH tests and the
`ChalmersSpecCoarsening.yaml` parse and execute tests. Coverage includes the
local update branches, the known drop/re-add conflict, the base permission
gate, slow recovery allowing p-down again, independent directional factors,
serialization, p retention, h reset, observation idempotence, delayed
activation, and an actual smooth-field refinement request from the wrapper.

A ten-PE tilted-Kerr smoke run completes 1M with the new controller active
from the start (`--target-shift-limit -1 --start-time 0 --end-time 1`). It
records six tightening and eight relaxation updates over nine evaluations.
`audit_native.py` independently reproduces all 648 directional target updates
and decisions. The mesh ends at 31,620 points in 24 elements, with 30
p-increases, no p-decreases or h changes, and positive excision outflow
(minimum sampled light speed 0.199779). Coarsening after recovery is checked
by the unit regression; this short evolution does not exercise it. The
existing coordinate-volume constraint integration passes its 28/36-point
quadrature check. This validates the implementation, not a reduction in
late-time cycles or cost; the matched 20M study is reported below.

The previous frozen `EvolveGhSingleBlackHole-stackfix` and its libraries are
unchanged. New trials use `build-release-debug/bin/EvolveGhSingleBlackHole`.
For the matched trial, leave `StartTime` at 10M and all other settings at
their defaults (choose a fresh output directory when repeating):

```sh
build-agents/bin/python-spectre tools/Amr/ChalmersTiltedKerr.py \
  --executable build-release-debug/bin/EvolveGhSingleBlackHole \
  --output build-agents/chalmers-gh-tilted-kerr/target-shift-20 \
  --policy pileup --target-shift-limit -1 --end-time 20 --processes 10
```

### Stored target: matched 20M tilted-Kerr evolution

The ten-PE `target-shift-20` run completes to 20M with the stored local target
enabled at 10M. The two-call cooldown, steep filter, damping, and physical
setup match `pileup-20`. The only YAML difference is the new
`TruncationErrorAdjustmentLimit: -1` option. This selects the stored-feedback
policy described above, including its replacement of instantaneous feedback
and its shared refinement/removal factor; the experiment does not isolate
each of those internal changes separately.

| Measurement | Previous forced pile-up | Stored target |
| --- | ---: | ---: |
| Final points / elements | 52,716 / 30 | 53,354 / 30 |
| Grid-point steps, 10--20M | 155,224,468 | 161,483,704 |
| Global steps, 10--20M | 3,408 | 3,570 |
| Realized directional p-additions / removals, whole run | 160 / 54 | 152 / 42 |
| Realized h-splits / joins | 6 / 0 | 6 / 0 |
| Remove-then-add pairs after activation | 24 | 14 |
| Pairs returning at the earliest allowed 0.2M | 12 | 2 |
| Final three-index constraint L2 | 0.00749437 | 0.00769070 |
| Final gauge constraint L2 | 0.000307604 | 0.000537807 |
| Final three-index constraint L2 outside r=2M | 0.00327821 | 0.00255188 |
| Maximum horizon-mass error | 1.69653e-6 | 1.69720e-6 |

The stored target reduces rapid cycling by 83.3% and all remove/add pairs
by 41.7%. Coarsening remains active: 38 of the 42 directional removals occur
at or after activation. The first h-split moves from 17.2M to 18M. All 14
re-additions are requested by the criterion, rather than induced by neighbor
balancing, and all realized p changes respect the 0.2M minimum interval.

This is a tradeoff, not an unqualified improvement. Work rises by 4.03%,
final points by 1.21%, and the final three-index and gauge norms by 2.62%
and 74.84%, respectively. The near-excision three-index norm in
`1.2<r<1.4M` rises from 0.00673273 to 0.00724799, while the exterior norm
decreases. Compared with ordinary coarsening, the stored-target run retains
a 76.14% reduction in work and an 88.95% reduction in the final three-index
norm. Work is the sum of grid points over steps beginning in `[10M,20M)`;
it excludes AMR/diagnostic overhead and is not a measured wall-clock speedup.

Both remaining rapid cycles are B4/B5 direction 0: removal at 15.5M and
addition at 15.7M. Their shift is already -1 at both decisions. The removal
uses SpEC's exceptional path for pile-up above 4.1: the minimum tensor
pile-up is about 4.109, and the estimated post-removal error is 7.157 times
the shifted target. The subsequent addition has error 7.045 times that
target. Thus the residual rapid reversal is permitted by the large-pile-up
exception, not caused by target relaxation or a cooldown violation. Six
of the 14 pairs, including these two, start with that exceptional path.

The native audit independently reproduces all 14,664 directional target
updates and decisions over 199 evaluations. It finds 58 tightening updates
and no relaxation updates. At the final evaluation, all 18 directions of
the six innermost elements have shift -1; other directions have zero shift.
This evolution therefore demonstrates continued coarsening with the new
policy, but does not exercise relaxation after constraint recovery. The
earlier 1M smoke run and unit regression cover relaxation separately.

All fifteen pre-10M volume snapshots match the ordinary control bitwise;
the 5M and 10M field differences from the previous forced run are also zero.
All independent field, constraint, radial-profile, and direct exterior
integrals pass the 28/36-point quadrature checks. Recomputed damped-harmonic
gauge sources agree with saved values to 2.14e-14 or better. The minimum
sampled outgoing excision light speed remains positive at 0.199779, and
the maximum dimensionless-spin error remains 2.81135e-6. The executable and
linked library hashes match the launch manifest after the run.

The exterior solution changes measurably: the directly integrated stored-
versus-previous-forced full-state L2 difference outside r=2M is 0.00194357
at 20M. Its exterior difference from the fixed fine reference is 0.0147375,
versus 0.0138346 for the previous forced run. The corresponding inner-shell
differences are 0.00826859 and 0.00826358. These are reference differences,
not established error estimates: the late interior reference convergence
problem remains unresolved. Reduced cycling alone does not establish an
accuracy improvement or readiness for BBH evolution. Retain this as an
experimental option; investigate the strong-pile-up exception and the gauge
constraint increase before choosing it as the default policy.

`TargetShiftSummary.json`, `TargetShiftCycleAudit.json`,
`TargetShiftBehavior.json`, `TargetShiftExteriorVerification.json`, and
`TiltedKerrTargetShift.png` / `.pdf` retain the results in
`build-agents/chalmers-gh-tilted-kerr`. `TargetShiftManifest.json` records the
launch configuration and binary/library hashes;
`TargetShiftEarlyRepeatability.json`
and `target-shift-20/TargetShiftAudit.json` retain the independent checks.
The `analyze_target_shift.py`, `audit_target_shift.py`,
`audit_target_shift_cycles.py`, `describe_target_shift.py`,
`verify_target_shift_exterior.py`, `summarize_target_shift.py`, and
`plot_target_shift.py` scripts preserve the analysis. No evolution code was
changed or rebuilt for this matched run.
