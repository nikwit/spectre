<!--
Distributed under the MIT License.
See LICENSE.txt for details.
-->

# Chalmers hp-refinement experiments

This directory contains the offline approximation study and the subsequent
scalar-wave and initial single-hole evolution experiments. These use SpECTRE's
Legendre Gauss-Lobatto nodes and nodal-to-modal transforms. The offline results
measure approximation potential; the evolution studies measure evolved errors
and wall time. Later GH constraint-feedback and SpEC-style controller studies
continue in
[GeneralizedHarmonicExperiments.md](GeneralizedHarmonicExperiments.md).
The cluster BBH starting point is [ChalmersBbh.md](ChalmersBbh.md).

The current cube runner defaults to the directional-energy variant described
in the directional-energy follow-up section. Earlier sections retain the
original stripe-based
experiments and their failures. The 1D criterion and input options are
unchanged.

## Reproduce

With a configured build that has Python bindings enabled:

```sh
cmake --build build-agents --target PyLinearOperators -j 2
ctest --test-dir build-agents --output-on-failure \
  -R '^Unit.LinearOperators.Python.Chalmers$'
MPLCONFIGDIR="$PWD/build-agents/matplotlib-cache" \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersStudy \
  --output build-agents/chalmers-study
```

Replace `build-agents` with the configured build directory. CMake installs the
Python source into the build package; rerun configuration after editing it if
`PY_DEV_MODE` is off. The study writes `comparisons.csv`, `summary.json`, and
three PNG figures. Results stay in the ignored build directory. The CSV retains
individual failures and undetermined decisions, rather than dropping them from
the experiment.

On the development machine, the existing build's extra Python packages were
found in `build-agents/lib/python3.13/site-packages`. Configuration needed that
directory in `PYTHONPATH`; automatic Python dependency bootstrapping was
disabled
to use those installed packages. A read-only ccache with a writable temporary
directory was used to rebuild bindings inside the sandbox.

## Estimator and explicit implementation choices

[Chalmers et al. (2019), section 2.4.6](https://doi.org/10.1016/j.jcpx.2019.100023)
fit individual Legendre coefficient sequences to `A exp(-sigma k)`. For highest
retained degree `p`, the squared indicator is

```text
eta^2 = 2 c_p^2 / (2p+1) + sum_{k=p+1}^infinity 2 A^2 exp(-2 sigma k)/(2k+1).
```

The directional result takes the largest indicator and smallest decay rate over
transverse modal sequences. Above tolerance, `sigma > 1` chooses p enrichment;
otherwise it chooses h refinement. This prototype extends that 2D construction
to 3D tensor products. The indicator is not a global error bound.

The paper does not fully specify the fitting details. This experiment makes the
following choices explicit; results should not be read as a reproduction of the
authors' implementation:

- Fit the natural logarithm of every significant coefficient amplitude, with
  equal weights. No tail-window optimization is performed.
- Use standard Legendre normalization `P_k(1)=1`, verified against the bindings.
- Ignore coefficients below `100 eps` times the largest coefficient of the
  tensor. Skip symmetry zeros while preserving their actual mode indices.
- Require three significant coefficients for a fit. Return `undetermined` when
  an active sequence has insufficient information; do not silently substitute h.
- Treat two absent final coefficients as a numerically terminated tail. This is
  a convention for sampled polynomials, not evidence that the analytic function
  is resolved.
- Return an infinite extrapolated error for nondecaying fits. Retain fit
  residuals as diagnostics, without automatically excluding a noise plateau.
- Sum the extrapolated tail with a hypergeometric expression. Unit tests check
  it against an independent long sum.

These estimates use reference-coordinate norms. Mapping Jacobians are included
in the *measured* errors for parent and child approximations, not in the local
spectral indicator. Initial tolerance is an absolute `1e-6`; the threshold of
one
is not tuned to the results.

## Experiment design

The 1D sweep contains 376 parent approximations on `[-1,1]`:

- 3, 4, 5, 6, 8, 10, 12, and 16 LGL points;
- Gaussians with widths 0.03 through 1.6 and five centers, including the origin
  and positions between nodes;
- sinusoids with five frequencies and two phases;
- constants and quadratics as resolved controls.

For each parent, sample the analytic function independently on an element with
one extra point and on two half-sized children with the original point count.
This tests approximation potential; it does not simulate interpolation of an
already erroneous evolved solution. Compute physical L2 errors on independent
256-point Gauss grids, then repeat with 512 points.

The 3D sweep contains nine cube approximations, with three directional candidate
comparisons each. It uses an anisotropic Gaussian with widths `(0.12, 0.5,
1.3)`,
the same Gaussian rotated by 45 degrees in the xy plane, and an oblique
sinusoid.
The extents are `(4,6,8)`, `(8,8,8)`, and `(12,10,8)`. Reference quadrature uses
64 and 128 points per dimension on each candidate element. A directional change
can have little effect on the global error when another direction dominates.
Therefore a directional hold with global error above tolerance is not itself a
false negative in 3D.

Quadrature verification requires that doubling the order changes each candidate
error by no more than `1e-7 * max(error, tolerance)`. The program exits with an
error if this check fails, and stores `quadrature_verified` in the summary.

For an explicit cost proxy, compare

```text
log(error_before / error_after) / number_of_added_grid_points.
```

Errors in this score are floored at `1e-13`. In direction d, p adds
`number_of_parent_points / extent[d]` points, while a binary h split adds an
entire parent's point count. This score measures approximation gain per added
point, not wall time, CFL restrictions, communication, or final cost at a fixed
accuracy. In particular, h can give a much smaller error even when this score
prefers p. The postprocessing also reports thresholds 0.5 and 1.5 for
sensitivity;
they are not calibrated or validated alternatives.

## Initial observations

The 1D reference run at tolerance `1e-6` gives:

| Quantity | Count |
| --- | ---: |
| All comparisons | 376 |
| h / p / hold / undetermined | 266 / 23 / 52 / 35 |
| Measured parent error above tolerance | 328 |
| Explicit h or p decisions among those | 283 |
| Selected candidate reduces error | 262 / 283 |
| Selected candidate reduces error by at least a factor of two | 193 / 283 |
| Selected candidate maximizes the gain-per-added-point proxy | 142 / 283 |
| Holds despite measured error above tolerance | 10 |

The outcomes demonstrate useful refinement choices but do not establish cost
optimality. For example, a Gaussian of width 0.2 centered at 0.23 on eight
points
has `sigma=-0.072` and selects h. Its L2 error changes from 0.213 to 0.0194
with h,
versus 0.190 with p. Conversely, width 0.8 at the same center and point count
has
`sigma=0.672` and still selects h: errors are `6.43e-4` initially, `3.16e-4`
after
p, and `3.78e-6` after h. The latter h candidate is more accurate, but p has the
slightly better gain per added point.

All 27 directional comparisons on the nine 3D cubes completed. On the
anisotropic
Gaussian with eight points per direction, the directional decay rates are
approximately `(-0.091, 0.290, 1.040)`, giving `(h, h, p)`. Splitting the narrow
x direction reduces the global error from 0.247 to 0.0647, whereas enriching x
by one point gives 0.269. Refining y or z barely changes the global error
because
x remains under-resolved. These are directional diagnostics, not a demonstration
of an optimal 3D adaptive mesh.

Doubling reference quadrature changes the L2 errors by at most `2.91e-14` in 1D
and `2.59e-13` in 3D; all comparisons pass the convergence check. The Chalmers
test file's 12 tests and the existing `Unit.Spectral.Python.Spectral` CTest
entry
pass. The initial run was made on branch `chalmers-amr`, based on develop commit
`f820072985bad5a2f54c9ac00e974e1892802a43`.

The ten false holds consist of seven undersampled Gaussians and three aliased
sinusoids. A width-0.05 Gaussian centered at zero on eight points has estimated
error about `1.05e-8` but actual L2 error about 0.250. The samples miss the
peak.
A sinusoid sampled as a constant has an analogous blind spot. Independent
quadrature and a dedicated unit test verify this failure; changing only the h/p
threshold cannot repair it.

The synthetic plateau test is also significant: adding a high-mode plateau to
an otherwise rapidly decaying spectrum can change the decision from p to h.
The rule alone does not distinguish numerical contamination from insufficient
spatial resolution. A high fit residual is observable, but no protective policy
has been chosen in this experiment.

Fit-window sensitivity, initial-resolution checks, and unreliable fits remain
open questions. These offline results do not establish recovery from an
under-resolved evolution or suitability for binary black holes.

## Native scalar-wave implementation

The C++ estimator is in `src/NumericalAlgorithms/LinearOperators/Chalmers.*`.
Its Python binding, `chalmers_estimate`, is compared against the independent
Python prototype in one, two, and three dimensions. The C++ infinite sum uses
a stable logarithmic expression at small sigma and a positive series with a
remainder bound otherwise. It does not add a SciPy or GSL dependency to C++.

`amr::Criteria::Chalmers` monitors selected real tensor components and chooses
p when `sigma > SmoothnessThreshold`, otherwise h, above
`AbsoluteTarget + RelativeTarget * max(abs(component))`. A split wins when
different monitored components disagree. `PRefineOnly: true` uses the same
estimator and adaptation phase as a control. `CoarseningFactor: 0.0` retains
the original refinement-only behavior. The later coarsening study below
introduces optional p reduction.
Undetermined active fits stop with a diagnostic; refinement caps also stop the
experiment rather than silently clamping or changing the h/p decision.

There is now a joint criterion type, `hp`. It runs once per
`EvaluateAmrCriteria` / `AdjustDomain` cycle, including when pure p criteria are
configured to run only in the `RefineMesh` event. Joint criteria are skipped by
that event. This avoids reevaluating the h/p choice after the first mesh change.
Existing projection, 2:1 balance, and refinement-limit machinery is reused.
The initial registration was limited to scalar-wave executables; the later
single-hole study adds the GH single-hole executable.

The first evolution pilot exposed a needed extra convention: an exponentially
small far tail can have a nondecaying fitted spectrum and hence an infinite
extrapolated indicator. Before fitting, the native criterion ignores a
component when `sqrt(2^Dim) * max(abs(component)) <= AbsoluteTarget`. For the
Legendre Gauss/LGL interpolant this bounds its reference-domain L2 norm.
It does **not** bound an unsampled analytic peak, and does not resolve the
sampling blind spot identified above. This amplitude guard is an implementation
choice, not a claim about the original paper. Unit tests cover negligible and
above-target nondecaying spectra.

## Traveling Gaussian experiment

`ChalmersEvolution.py` generates inputs, launches `EvolveScalarWave1D`, and
postprocesses its native HDF5 output. The setup is deliberately nonuniform:

- `Psi = exp(-((x + 2 - t)/0.3)^2)`, with the corresponding analytic Pi and Phi.
- `RotatedIntervals` on `[-4,0]` and `[0,8]`, eight elements in each block.
  The left block starts with 16 points per element, the initially empty right
  block with five. The right block reverses logical orientation.
- The initial physical combined L2 error is `5.08e-10`, below the tightest
  requested target. The pulse travels from x=-2 to x=4 by t=6.
- Global Dormand-Prince fifth-order stepping, fixed `dt=0.0005`, adaptation
  every 0.02 time units, diagnostics every 0.2. No spectral filtering and no
  coarsening. Both boundaries use the analytic solution.
- All three fields are monitored. Absolute targets are `1e-4`, `1e-6`, `1e-8`;
  the relative target is zero. Compare p-only, hp with threshold 1, and hp with
  threshold 0.5. The p-only case changes only the h/p selection.
- Maximum refinement level 7; maximum polynomial degree 20 (21 LGL points).
  A pilot cap of level 5 was too restrictive and its failed run is retained.
- The fixed reference has 32 elements with 20 points each and the same time
  step. It achieves a final combined field error `2.75e-13` against the exact
  solution.

Errors below are **physical, unnormalized L2 norms**, integrated on independent
64- and 128-point Gauss grids in each element with the absolute mapping
Jacobian. The combined norm is `sqrt(||e_Psi||^2 + ||e_Pi||^2 + ||e_Phi||^2)`;
individual field errors and `||d_x Psi - Phi||` are also stored. The native
`ObserveNorms` output is nodal RMS and is not substituted for these integrals.
Doubling the independent quadrature changes the reported field norms by less
than `6e-16` in the resolved runs.

Results on the local Apple Silicon build:

<table>
<tr>
<th>
Target
</th>
<th>
Strategy
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
Combined L2 error
</th>
<th>
Evolution seconds
</th>
</tr>
<tr>
<td>
1e-4
</td>
<td>
p only
</td>
<td>
16
</td>
<td>
224
</td>
<td>
5.14e-5
</td>
<td>
13.1
</td>
</tr>
<tr>
<td>
1e-4
</td>
<td>
hp, threshold 1
</td>
<td>
55
</td>
<td>
421
</td>
<td>
8.14e-5
</td>
<td>
27.5
</td>
</tr>
<tr>
<td>
1e-4
</td>
<td>
hp, threshold 0.5
</td>
<td>
41
</td>
<td>
353
</td>
<td>
8.63e-5
</td>
<td>
22.9
</td>
</tr>
<tr>
<td>
1e-6
</td>
<td>
p only
</td>
<td>
stopped
</td>
<td>
stopped
</td>
<td>
degree limit reached
</td>
<td>
—
</td>
</tr>
<tr>
<td>
1e-6
</td>
<td>
hp, threshold 1
</td>
<td>
56
</td>
<td>
505
</td>
<td>
8.29e-7
</td>
<td>
29.2
</td>
</tr>
<tr>
<td>
1e-6
</td>
<td>
hp, threshold 0.5
</td>
<td>
56
</td>
<td>
505
</td>
<td>
8.29e-7
</td>
<td>
29.1
</td>
</tr>
<tr>
<td>
1e-8
</td>
<td>
p only
</td>
<td>
stopped
</td>
<td>
stopped
</td>
<td>
degree limit reached
</td>
<td>
—
</td>
</tr>
<tr>
<td>
1e-8
</td>
<td>
hp, threshold 1
</td>
<td>
102
</td>
<td>
864
</td>
<td>
9.82e-9
</td>
<td>
46.6
</td>
</tr>
<tr>
<td>
1e-8
</td>
<td>
hp, threshold 0.5
</td>
<td>
59
</td>
<td>
601
</td>
<td>
4.83e-9
</td>
<td>
30.7
</td>
</tr>
</table>

The two p-only failures have last recorded snapshots at t=4.2 and t=1.8,
respectively, with degree 20 reached. This is a limitation of that starting
mesh, estimator and degree cap, not evidence that p refinement generally cannot
solve this smooth problem.

At `1e-8`, threshold 0.5 produces 43 h splits instead of 86, uses 30% fewer
final
points and reduces the error by about a factor of two. Near the pulse it keeps
degree 9 on elements twice as wide as the original rule's degree-7 elements.
At `1e-6`, the two thresholds give identical solution and mesh histories. At
`1e-4`, p-only is both cheaper and more accurate in the combined norm. The
indicator is a local refinement target, not a guaranteed global error bound or
cost optimizer.

The fixed reference is also cheaper than either tight hp run in this small 1D
problem and substantially more accurate. These results support trying a lower
threshold, but do not demonstrate an overall AMR speedup over an appropriate
fixed grid. Timings use one Charm++ PE, an existing Debug build with effective
`-Oz -DSPECTRE_DEBUG -g0` flags, and the elapsed observer wall time between
first
and last observations. They include adaptation and diagnostics but exclude most
startup cost. They are local exploratory timings, not production scaling data.

Two further runs without volume output give 46.8–48.1 s for hp with threshold
1 at `1e-8`, 30.3–30.9 s for threshold 0.5, and 27.1–27.4 s for the fixed
reference. Thus the lower threshold is about 35% faster than the original rule
in this test; volume output is not responsible for the difference. Reductions
and the AMR domain checks remain enabled in these timing repeats.

Halving dt to 0.00025 at the same physical adaptation cadence preserves the
mesh histories for the fixed reference, p-only at `1e-4`, and both hp cases at
`1e-8`. Raw final-state differences are about `3.69e-11`. Floating-point slab
accumulation shifts the final time by `2.19e-12`; subtracting that known
analytic
phase shift reduces the differences to `0.95–1.23e-13`. Thus temporal error is
well below the measured adaptive spatial error. The script records both
comparisons. The largest final `dt/minimum_nodal_spacing` is 0.526 in hp at
threshold 1, below the Dormand-Prince stability factor used by SpECTRE.

`Fits.csv` records every sampled element's sigma, indicator and fit residual.
`Spectra.csv` stores Pi coefficients, analytic-interpolant coefficients on the
same mesh, and their difference. The plotted pulse-region spectra at `1e-8`
continue to decay through the retained modes; numerical modal errors there are
around `1e-10` or smaller. This checks the resolved baseline, not every possible
noise-contamination failure.

The separate under-resolved stress run starts with only six points per left
element. Its initial combined error is `3.79e-3` at the requested target
`1e-8`. With threshold 0.5 it requests h refinement beyond level 7 and stops
before the first t=0.2 diagnostic snapshot. This is a failed stress test, not a
successful recovery: it reinforces the need for an initial-resolution check
and a policy for unreliable spectra before more demanding evolutions.

The following 3D stage tests a resolved moving Gaussian on Cartesian cubes,
including axis-aligned and oblique propagation, with fixed-grid and p-only
controls. Curved cubed-sphere domains and black-hole evolutions remain later
stages; this branch has not added the criterion to generalized harmonic
executables.

## Reproducing the evolution study

Build `EvolveScalarWave1D`, `PyLinearOperators`, `Test_LinearOperators`, and
`Test_ParallelAmr`. Run the affected tests:

```sh
ctest --test-dir build-agents --output-on-failure \
  -R 'Unit\.Amr\.|Unit.LinearOperators.(Chalmers|Python.Chalmers)$'
ctest --test-dir build-agents --output-on-failure \
  -R '^InputFiles.ScalarWave.(Chalmers1D|PlaneWave1D).yaml.(parse|execute)$'
```

The new `tests/InputFiles/ScalarWave/Chalmers1D.yaml` is a short native
regression
run through t=1. The offline/C++ tests cover analytic exponential tails,
normalization, parity zeros, noise plateaus, tiny and large amplitudes, 3D
directions, Python/C++ agreement, factory parsing, serialization, joint flags,
topology restrictions, and the existing pure h/p scheduling behavior.
All 31 selected CTest entries pass (27 unit entries and four input-file
parse/execute entries). One earlier legacy parse invocation failed during
Charm++ startup; its immediate retry and the final complete selection passed.

For the full comparison (use an output directory that does not already exist):

```sh
for tolerance in 1e-4 1e-6 1e-8; do
  for strategy in p hp1 hp05; do
    build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersEvolution run \
      --executable build-agents/bin/EvolveScalarWave1D \
      --output build-agents/chalmers-evolution/reproduced/$strategy-$tolerance \
      --strategy "$strategy" --tolerance "$tolerance"
  done
done
MPLCONFIGDIR="$PWD/build-agents/matplotlib-cache" \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersEvolution analyze \
  build-agents/chalmers-evolution/reproduced/* \
  --output build-agents/chalmers-evolution/reproduced-results
```

The p-only limit failures intentionally return nonzero; the shell loop should
continue to the other cases. `--strategy fixed` runs the reference. Repeat with
`--dt 0.00025` for time convergence; compare the two directories with the
`difference FIRST SECOND --output comparison.json` subcommand. Use
`--no-volume` for timing repeats and `--underresolved` to reduce the initial
left-block resolution from 16 to six points for the separate stress test.

Each run retains its YAML input, executable log, run metadata and HDF5 data.
Postprocessing writes `History.csv`, `Fits.csv`, `Spectra.csv`, `Mesh.csv` and
`Analysis.json`. The combined output includes `Summary.json`,
`EvolutionComparison.png`, and `MeshAndSpectra.png`. Failed runs remain in the
summary with a reason and last observed mesh; they are excluded from final-error
plots, not counted as successful evolutions.

## Three-dimensional cube evolution

`ChalmersCube.py` runs the native `EvolveScalarWave3D` executable and compares
fixed grids, p-only refinement, and hp refinement at sigma thresholds 1 and 0.5.
The exact solution is
`Psi = exp(-((n.x - t + 1.8)/0.35)^2)`, with the corresponding analytic Pi and
three Phi components. The directions are `(1,0,0)` and `(1,2,3)/sqrt(14)`.

The domain is `[-1,1]^3`, initially eight affine blocks with no subdivisions.
Adaptive runs start at degree 13 in resolved directions. The downstream half
of the axis-aligned domain, or the positive octant of the oblique domain,
starts at degree 4. The two transverse degrees in the axis control are 2.
The initial combined L2 errors are `1.03e-7` and `8.05e-8`, respectively,
below the tightest target `1e-5`. The initially resolved tail enters through
analytic boundary data and evolves to t=2.

The axis symmetry control uses periodic transverse boundaries and analytic
boundaries in x (`--axis-periodic`). The oblique problem uses analytic data on
all six faces. Early axis runs with analytic data on all faces are retained
separately: their boundary discretization produces transverse errors and
therefore additional refinement, despite the transverse symmetry of the exact
solution. They are not included in the axis symmetry comparison.

Both fixed-step and CFL-controlled runs use global Dormand-Prince stepping,
physical adaptation times every 0.04, and diagnostics every 0.4. `StepToTimes`
lands exactly on that schedule. The fixed maximum step is 0.0005; landing steps
are sometimes shorter. There is no coarsening or filtering. Refinement stops
explicitly at level 3 (minimum element width 1/8) or polynomial degree 20.
A stop at these limits is retained as a failed experiment, not counted as a
successful error/cost result.

Errors use independent tensor-product Gauss quadrature with 24 and 32 points
per direction, including the physical Jacobian. The combined norm includes
all five fields; the reduction constraint is `||grad Psi - Phi||`. Cost plots
use the **largest observed** combined error at the six volume snapshots,
rather than assuming the requested local tolerance is attained globally.
This sampled maximum is not a bound between observations. Time-step comparisons
integrate the difference of the two numerical polynomials over the common
partition, including when the adaptive meshes differ.

### Additional multidimensional convention

The first 3D pilots exposed a problem with applying the stripe fit literally:
roundoff-sized transverse coefficient sequences can have nondecaying spectra,
whose extrapolated infinite tails trigger repeated refinement. The experimental
criterion now omits a whole stripe if its one-dimensional reference L2 norm is
at most `AbsoluteTarget / sqrt(2^(Dim-1) * number_of_stripes)`. Orthogonality of
the transverse Legendre modes bounds the combined norm of these omitted
stripes by `AbsoluteTarget` for that direction. This is an additional convention
of this implementation, not a prescription attributed to Chalmers et al.

The bound applies to the sampled interpolant, not the unsampled analytic
function or the global evolution error. The raw estimator remains available
with zero stripe budget, and the 1D criterion is unchanged. Active insufficient
or nondecaying fits are still handled explicitly; the guard does not silently
replace them with a preferred h/p decision. The original failed pilots and
subsequent guarded runs are retained under `build-agents/chalmers-cube/`.

### Cube findings

For the axis symmetry control at target `1e-4`, threshold 0.5 uses 2,340 final
points and 22.8 s, compared with 3,168 points and 29.2 s at threshold 1.
The maximum observed combined errors are `6.67e-5` and `6.18e-5`, respectively.
Thus lowering the threshold reduces the point count by 26% and elapsed time by
22% in this case. Both remain more expensive than the degree-13 fixed grid:
1,008 points, 9.50 s and maximum observed error `7.03e-5`. At `1e-3` the two hp
thresholds give identical mesh and solution histories. At `1e-5` both stop with
undetermined transverse fits near t=1.92.

All three oblique p-only runs complete. Targets `1e-3`, `1e-4`, `1e-5` give
final combined errors `1.19e-4`, `9.89e-6`, `1.45e-6`. At target `1e-4`, the
initially coarse positive octant ends with degrees `(8,10,12)`, consistent with
the wave-vector components increasing from x to z. In contrast, all six
oblique hp runs hit the level-3 cap, between t=1.0 and t=1.92. Their last saved
snapshots still have errors below the requested tolerance, so these failures
are not evidence of a visibly diverging evolution. They demonstrate excessive
refinement pressure and sensitivity to the chosen limits in this prototype.
Raising the cap has not been tested here.

The coefficient diagnostics distinguish the envelope plotted in
`CoefficientSpectra.png` from the worst individual stripe used for the decision.
For example, on the final downstream oblique p-only element, the analytic
interpolant itself gives minimum stripe sigmas about `(0.904,0.301,0.316)`.
The evolved field gives `(0.363,0.268,0.161)`, with maximum log-fit residuals
around 2.1. Thus small fitted sigma is not by itself evidence of a physically
nonsmooth solution. This is a reason to revisit fit windows and how stripes
are combined before applying the criterion to curved domains.

CFL stepping provides a second robustness check, using the same physical
adaptation times. Both axis hp runs at CFL factor 0.5 stop near t=0.96 because a
transverse fit is undetermined. At factor 0.25, threshold 0.5 completes with the
same final h mesh and point count as the fixed-step run, but threshold 1 stops
near t=1.96. These failures survive the whole-stripe norm guard: that guard
does not resolve every case with too few significant modes in a nearly
constant direction. They are recorded, not replaced by automatic h/p choices.

For successful axis runs, halving the CFL factor changes the final state by
`4.52e-11` on the degree-19 fixed grid and `1.50e-9` for p-only. The successful
hp threshold-0.5 run at CFL factor 0.25 differs from the fixed-step result by
`2.82e-10`. For the oblique wave, the corresponding CFL-halving differences
are `2.94e-9` for the fixed reference and `1.06e-8` for p-only. The p-only
CFL-0.25 result differs from fixed stepping by `3.96e-8`, below 0.5% of its
spatial error; its polynomial-degree pattern differs slightly, so this
comparison includes changes in the AMR path. Both oblique hp thresholds hit
similar h-refinement limits with fixed stepping and both CFL factors. The very
accurate fixed reference
at CFL factor 0.5 is time-error limited, but remains far below the adaptive
error targets. These successful comparisons do not establish time convergence
of the hp configurations that stopped.

The cube study therefore supports a limited benefit from lowering sigma in
the axis control, and **does not establish a robust or cost-effective 3D hp
method**. The next step is to fix and test the handling of nearly constant
transverse directions and the selection of a representative decay rate for
oblique data, then repeat these controls. Moving directly to cubed spheres or
black holes would obscure these already-reproducible problems.

At matched observed accuracy `1e-4`, the cheapest successful fixed-step cases
in the sampled sweep are:

| Wave | Fixed grid | p only | hp, sigma 1 | hp, sigma 0.5 |
| --- | ---: | ---: | ---: | ---: |
| Axis | 9.50 s | 9.76 s | 29.2 s | 22.8 s |
| Oblique | 44.5 s | 63.9 s | not attained | not attained |

These are discrete sampled comparisons, not interpolated optimal costs.
For the representative CFL-0.5 runs, the oblique degree-11 fixed grid costs
5.70 s at error `7.32e-5`, while p-only costs 11.9 s at error `1.09e-5`.
The fixed grids are cheaper at the shared `1e-4` accuracy requirement in
these controls. The hp CFL-0.5 cases do not complete. All timings remain
single-PE local experiments in the optimized Debug build, not production
performance measurements.

The main comparison records 64 attempts: 45 completed evolutions, 18 stops at
refinement limits or undetermined fits, and one Charm++ startup failure.
The startup failure was an auxiliary all-analytic-boundary axis case and its
successful retry is retained separately. The count excludes the earlier pilot
and diagnostic runs. `results/Summary.json` and `ExperimentResults.csv` retain
all attempts; `MatchedAccuracy.csv` excludes incomplete runs. The results
folder also contains `AccuracyCost.png`, `DirectionalRefinement.png`,
`CoefficientSpectra.png`, time-step comparisons, independent-quadrature checks,
source and executable hashes in `BuildMetadata.json`, and `CTest.log`.

### Reproducing the cube runs

Build `EvolveScalarWave3D`, `PyLinearOperators`, `Test_LinearOperators`,
`Test_ParallelAmr`, and `cli`. Run, for example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersCube run \
  --executable build-agents/bin/EvolveScalarWave3D \
  --output build-agents/chalmers-cube/example-axis \
  --orientation axis --axis-periodic --strategy hp05 --tolerance 1e-4
```

Repeat with `--orientation oblique`, strategies `p`, `hp1`, `hp05`, and
tolerances `1e-3`, `1e-4`, `1e-5`. Fixed controls use `--strategy fixed
--fixed-p DEGREE`, with degrees 7, 9, 11, 13, 15 and 19. For the CFL
comparison add
`--stepping cfl --cfl 0.5`, and halve the factor to 0.25 for time convergence.
Use `compare FIRST SECOND --output comparison.json` to compare final fields,
and `summarize RUN... --output RESULTS` to create accuracy/cost figures and
matched-accuracy tables. Run directories are never overwritten.

`tests/InputFiles/ScalarWave/Chalmers3D.yaml` is a short native regression
through t=1.2. The new tests also cover the stripe guard, its use by the 3D
criterion, physical quadrature normalization and polynomial differences on
nonmatching anisotropic partitions. All 35 affected CTest entries pass:
27 unit entries and eight input-file parse/execute entries, including the
existing 1D and 3D plane-wave inputs. Curved domains and generalized harmonic
black-hole evolutions remain untested; this branch still registers the new
criterion only in scalar-wave executables.


## Directional-energy follow-up

The cube failures were isolated with two changes tested separately. These are
implementation choices for this experiment, not additional rules attributed to
Chalmers et al.

1. **Guard the nonconstant part.** A large k=0 coefficient can coexist with
   negligible directional variation. The stripe guard now bounds only modes
   k>=1. It therefore recognizes a nearly constant direction even if roundoff
   leaves too few significant coefficients for a decay fit. Four previously
   failing axis cases now complete using this change alone: both sigma
   thresholds at target `1e-5` with fixed stepping, and at target `1e-4` with
   CFL
   factor 0.5. The two representative oblique hp cases still reach the h limit.
2. **Fit a directional L2 spectrum.** The new multidimensional variant combines
   transverse coefficients by their Legendre norms before fitting. In direction
   d, define

   ```text
   S_k^2 = sum_{j_perp} c_{k,j_perp}^2 product_{i!=d} 2/(2 j_i + 1).
   ||u||_reference^2 = sum_k 2 S_k^2/(2 k + 1).
   ```

   Both the error estimate and sigma now come from the same fit to S_k, instead
   of independently taking the largest stripe error and smallest stripe sigma.
   Transverse coefficients cannot cancel in this sum. The nonconstant norm
   `sqrt(sum_{k>=1} 2 S_k^2/(2k+1))` is compared with the absolute plus relative
   target before fitting. A significant unresolved direction still produces
   an explicit insufficient/nondecaying result; no arbitrary h/p fallback was
   introduced.

The energy variant uses a full reference-element L2 normalization. For an
axis-aligned wave in 3D this gives twice the corresponding 1D modal amplitude,
so equal local targets do not make the old and new criteria equally strict.
Comparisons should use attained physical errors. The fits still use all
significant modes with the original parity and terminated-tail conventions.
The variation guard is a bound on the sampled polynomial's variation, not a
bound on aliasing, the unsampled function, or the global evolution error.

In 2D and 3D the criterion now takes the explicit option
`UseDirectionalEnergy: true` or `false`. False selects individual stripes with
the improved nonconstant-part guard. `ChalmersCube run --aggregation energy`
is the default; `--aggregation stripes` selects that control. The option is
absent in 1D and its numerical behavior is unchanged. The raw
`chalmers_estimate`
Python binding remains available; `chalmers_energy_estimate` exposes the new
fit separately. Historical whole-stripe diagnostics are retained, and the
current analyzer refuses to silently replace them with fits from a different
algorithm.

The follow-up data are in `build-agents/chalmers-cube-v2/`. `variation-only/`
contains the six cases testing the first change alone. `energy/` contains the
three-tolerance CFL sweep, halved-CFL comparisons and representative fixed-step
runs. The original executable and kernel sources were preserved in `baseline/`,
and the variation-only executable was also saved. Each energy run records its
executable hash. The final reports retain failed control cases alongside the
successful energy cases.

New tests cover the full Legendre L2 normalization, a tiny transverse quadratic
on top of a large constant, preservation of significant unresolved variation,
parity, non-finite inputs, and invariance of directional energy under orthogonal
transformations of transverse modal coefficients. A mixed-spectrum criterion
test checks that the explicit selector and its serialized value choose p for
the energy fit and h for the stripe fit. The 3D input regression explicitly
selects the energy variant.


### Follow-up results

All **27 energy-variant evolutions complete**: 18 runs at CFL factor 0.5
(two orientations, three strategies, three tolerances), six runs at CFL factor
0.25 and target `1e-4`, and three fixed-step controls. The domain, initial data,
physical AMR cadence, degree cap and level-3 cap were retained. The largest
observed combined physical L2 error in each of the 18 main runs is below its
requested local target; this is an empirical result of these tests, not a
new error bound. The axis runs preserve the transverse degrees and widths.

At local target `1e-4` with CFL factor 0.5:

<table>
<tr>
<th>
Wave
</th>
<th>
Strategy
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
Maximum observed error
</th>
<th>
Evolution seconds
</th>
</tr>
<tr>
<td>
Axis
</td>
<td>
p only
</td>
<td>
8
</td>
<td>
1,080
</td>
<td>
5.82e-5
</td>
<td>
1.24
</td>
</tr>
<tr>
<td>
Axis
</td>
<td>
hp, sigma 0.5
</td>
<td>
24
</td>
<td>
2,052
</td>
<td>
5.33e-5
</td>
<td>
3.70
</td>
</tr>
<tr>
<td>
Axis
</td>
<td>
hp, sigma 1
</td>
<td>
44
</td>
<td>
3,168
</td>
<td>
3.61e-5
</td>
<td>
9.96
</td>
</tr>
<tr>
<td>
Oblique
</td>
<td>
p only
</td>
<td>
8
</td>
<td>
20,048
</td>
<td>
2.87e-5
</td>
<td>
11.1
</td>
</tr>
<tr>
<td>
Oblique
</td>
<td>
hp, sigma 0.5
</td>
<td>
18
</td>
<td>
30,120
</td>
<td>
2.17e-5
</td>
<td>
22.7
</td>
</tr>
<tr>
<td>
Oblique
</td>
<td>
hp, sigma 1
</td>
<td>
47
</td>
<td>
47,715
</td>
<td>
2.07e-5
</td>
<td>
33.3
</td>
</tr>
</table>

In the tighter oblique test (`1e-5`), sigma 0.5 uses 45,907 points and 34.2 s,
versus 90,841 points and 111 s for sigma 1. Their maximum observed errors are
`2.97e-6` and `2.54e-6`. The lower threshold therefore gives about 49% fewer
points and 69% lower elapsed time at similar accuracy in this case. These are
local exploratory timings; the point counts are the more portable comparison.
Fixed degree-13 grids and p-only refinement remain cheaper at the shared
`1e-5` accuracy requirement. This is progress in robustness, not evidence that
hp AMR is the optimal strategy for these plane-wave controls.

Halving the CFL factor preserves all six recorded mesh histories. Final-state
L2 differences range from `1.64e-9` to `1.84e-8`, below 0.07% of the respective
final errors. The three fixed-step comparisons with CFL factor 0.25 differ by
`7.31e-11` to `6.47e-10`, with identical final h meshes. The comparisons use
polynomial integration, not differences between scalar error norms. Independent
24^3/32^3 quadrature checks are also retained. Constraint norms remain recorded
separately and decrease with tighter tolerances; the field target does not
provide a separate guarantee for the derivative constraint.

All 35 affected CTest entries pass, including the 1D regressions and the
updated 3D regression. There are 17 Python test methods within the Python
CTest entry. A new selector test caught an initially unwired input option;
that was fixed before final validation. Two full native runs with the final
executable confirm that `UseDirectionalEnergy: false` reproduces the
variation-only control and `true` reproduces the energy pilot.

Reproduce a representative new run with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersCube run \
  --executable build-agents/bin/EvolveScalarWave3D \
  --output build-agents/chalmers-cube-v2/reproduced-oblique \
  --orientation oblique --strategy hp05 --aggregation energy \
  --tolerance 1e-4 --stepping cfl --cfl 0.5
```

Use `--aggregation stripes` for the improved stripe control. The result folder
contains the accuracy/cost and directional-refinement figures, `EnergyRuns.csv`,
`Ablation.csv`, `CflComparisons.json`, fixed-step comparisons, source/executable
hashes and the final CTest log. The next physics test can now be a scalar wave
on a curved cubed-sphere domain, retaining independent errors, constraint
monitoring, fixed-grid/p-only controls and step-size checks. Generalized
harmonic black-hole runs still require later validation.


## Curved cubed-sphere scalar wave

`ChalmersSphere.py` adds a native `EvolveScalarWave3D` experiment on three
concentric cubed-sphere shells, with radial boundaries at 1, 2, 3 and 4.
Each shell has six equiangular, genuinely curved hexahedral wedges. The
interior is excised and both spherical boundaries impose the exact analytic
solution. There is no central cube, moving map or black-hole system in this
test. The criterion and executable are unchanged from the final
`UseDirectionalEnergy: true` cube study.

The smooth regular spherical wave is

```text
F(s) = exp(-((s - 1.4) / 0.4)^2)
Psi(x,t) = (F(r-t) - F(-r-t)) / r
Pi = -partial_t Psi,   Phi_i = partial_i Psi.
```

The pulse center travels from radius 1.4 to 3 over `0 <= t <= 1.6`. Adaptive
runs start with `[14,14,18]`, `[14,14,16]` and `[14,14,8]` points in the three
shells, in angular/angular/radial order. Their initial physical combined L2
interpolation error is `1.23e-7`, checked before interpreting evolution errors.
This uses more initial resolution than a preliminary offline interpolation
trial, which had an unacceptable `8.78e-5` error in the second shell. The
preliminary trial was not evolved or included in the timing comparison.

The experiment retains the global Dormand-Prince 5 integrator, physical AMR
cadence 0.04, refinement-only policy, degree cap 20, h-level cap 3, and absence
of filtering. All runs land on the same scheduled times; independent volume
observations are at 0, 0.4, 0.8, 1.2 and 1.6. Fixed-grid controls have the same
18 blocks and uniform degrees; p-only uses the same estimator and initial
mesh as hp. Timings are single-PE local exploratory measurements, including
AMR and output but excluding the Python analysis.

Curved geometry requires different diagnostics from the Cartesian study:

- Read the exact serialized domain from the volume file and construct each
  element map. Check mapped evolution nodes against stored coordinates.
- Reconstruct each field's logical Legendre polynomial, evaluate it at
  independent 24^3 and 32^3 Gauss nodes, and weight the errors by the exact
  physical Jacobian determinant. Integrate over all wedges and shells. The
  reported field error combines the five physical L2 norms in quadrature;
  it is not volume normalized. Maximum errors are maxima over the five saved
  observations, not continuous-time bounds.
- Compute `partial_i Psi - Phi_i` with the inverse coordinate Jacobian, and
  check quadrature convergence of this derivative constraint as well as the
  fields. Also verify the integrated shell volume against `4*pi*(4^3-1)/3`.
- Compare CFL runs by integrating the difference between their reconstructed
  fields on a common block-logical partition, using the curved physical
  Jacobian. This supports different h and p meshes. Differences of scalar
  error norms are not used as estimates of time error.

Unit tests verify the analytic wave equation, physical gradients of `r^2` on
all six wedge orientations, shell volume, and comparisons between different
radial h meshes and polynomial degrees. `ChalmersSphere.yaml` is a shorter
native curved-domain parse/execute regression.

The criterion still estimates an **unweighted reference-element** spectral
error. It is not a physical-volume-weighted estimator, a global error bound,
or a separate bound on the derivative constraint. This symmetric radial test
exercises curved geometry and radial adaptation; it does not test angular
transport of a localized pulse. Those limitations matter before proceeding
to generalized harmonic or black-hole evolutions.

Reproduce a run with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersSphere run \
  --executable build-agents/bin/EvolveScalarWave3D \
  --output build-agents/chalmers-sphere/reproduced-hp05 \
  --strategy hp05 --tolerance 1e-4 --cfl 0.5
```

Use `--strategy p`, `hp1`, or `fixed --fixed-p 13` for controls, and
`--cfl 0.25` for the time-step check. The `analyze`, `compare` and `summarize`
subcommands reproduce the diagnostics and figures from saved runs. Inputs,
logs, volume data, mesh histories, directional fits and executable hashes
are retained under `build-agents/chalmers-sphere/`.


### Curved-shell results

All **20 study evolutions complete**: twelve adaptive runs (p-only and both
sigma thresholds at local targets `1e-3`, `1e-4`, `1e-5`, `1e-6`), four uniform
fixed-grid runs (degrees 9, 11, 13, 15), and four halved-CFL controls. No fit
failure or refinement-limit failure occurs. Both maximum observed field errors
and maximum observed constraint norms decrease as the local target is
reduced. The angular mesh stays at degree 13 without h refinement in every
adaptive run. All adaptation occurs radially.

At the tightest local target, `1e-6`, compared with the finest fixed control:

<table>
<tr>
<th>
Strategy
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
Maximum observed field error
</th>
<th>
Final constraint L2
</th>
<th>
Evolution seconds
</th>
</tr>
<tr>
<td>
p only
</td>
<td>
18
</td>
<td>
57,624
</td>
<td>
1.44e-6
</td>
<td>
2.40e-5
</td>
<td>
37.3
</td>
</tr>
<tr>
<td>
hp, sigma 0.5
</td>
<td>
24
</td>
<td>
64,680
</td>
<td>
2.48e-6
</td>
<td>
1.22e-5
</td>
<td>
41.7
</td>
</tr>
<tr>
<td>
hp, sigma 1
</td>
<td>
42
</td>
<td>
97,608
</td>
<td>
1.99e-6
</td>
<td>
1.34e-5
</td>
<td>
83.9
</td>
</tr>
<tr>
<td>
Fixed degree 15
</td>
<td>
18
</td>
<td>
73,728
</td>
<td>
8.33e-7
</td>
<td>
4.31e-6
</td>
<td>
39.1
</td>
</tr>
</table>

Sigma 0.5 uses about **34% fewer points and 50% less elapsed time** than sigma
1 at similar field accuracy in this tight case. It stops at radial h level 1;
sigma 1 reaches level 2. At the two intermediate targets (`1e-4`, `1e-5`) both
thresholds instead produce the same recorded meshes and error norms. The
lower threshold therefore helps at the tightest target, but does not alter
every refinement decision.

P-only is cheaper than hp in the tight comparison and has smaller field
error, but has a larger derivative-constraint error. Fixed degree 15 has
smaller field and constraint errors at a runtime comparable to the adaptive
p-only and sigma-0.5 cases. Fixed degree 13 also matches the `1e-5` local-target
hp runs' roughly `2e-5` physical field accuracy at about half the runtime.
There is no demonstrated hp cost advantage for this smooth radial problem.
The point counts and error norms are more portable evidence than the small
runtime difference between p-only and fixed degree 15. Coarsening remains
disabled, so resolution left behind the pulse is retained.

The requested local target is not a global physical bound: at target `1e-4`,
for example, the maximum physical field errors are `4.72e-4` (p only) and
`1.34e-4` (either hp threshold). `MatchedAccuracy.csv` compares only runs that
actually meet each physical error requirement; it does not equate local
criterion targets with achieved accuracy or interpolate between experiments.

The four CFL checks use factor 0.25 instead of 0.5 for the three adaptive
strategies at target `1e-4` and fixed degree 13. All four retain identical
recorded h/p mesh histories. Adaptive final-state L2 differences are about
`1.26e-10`; the fixed-grid difference is `1.57e-9`. The largest difference is
less than **0.013% of the respective final spatial error**. This supports
spatial-error dominance for those checked cases; the tighter-target runs were
not separately repeated at half CFL.

Across all 20 runs, changing independent quadrature from 24^3 to 32^3 changes
any reported field/constraint norm by at most `2.63e-15`. Mapped evolution
coordinates agree to `8.89e-16`, integrated shell volumes to `1.71e-13`, and
the independently coded analytic initial state agrees with native nodal data
to `6.67e-16`. These checks validate the measurement pipeline; they do not
turn the estimator into a rigorous error bound.

`results/AccuracyCost.png` compares attained errors against time and points.
`results/RadialRefinement.png` shows radial degrees and h splits in a
representative wedge at target `1e-6`. `Summary.json`, `RunSummary.csv`,
`CflComparisons.json`, `Validation.json` and `BuildMetadata.json` retain the
numerical results and provenance. `run_sweep.py`, `run_tighter.py` and
`finalize.py` under the study directory reproduce the local experiment matrix.

The next test should break spherical symmetry: propagate an off-center or
oblique smooth pulse across angular wedge interfaces with the same controls.
That would test angular h/p choices and transport between differently oriented
blocks before progressing to generalized harmonic and black-hole systems.

All seven affected CTest entries pass: the 1D, cube and curved-shell input
parse/execute tests, plus the Python entry containing 20 test methods. The
curved-shell regression takes about 19 seconds locally. Its first attempt
exceeded the harness's default 10-second timeout; the input now requests an
explicit 60-second timeout, and the successful rerun is in `results/CTest.log`.
The initial timeout log is preserved separately. Black, isort and
`git diff --check` also pass.


## Oblique transport across angular wedge interfaces

`ChalmersWedges.py` tests a nonspherical, smooth scalar wave on a single
cubed-sphere shell, `1 <= r <= 2`, made of six equiangular wedges. The analytic
solution is an oblique Gaussian plane pulse:

```text
n = (1, 2, 3) / sqrt(14)
u = n.x - t + 2.8
Psi = exp(-(u/0.8)^2)
Pi = -2*u/0.8^2 * Psi,   Phi_i = n_i * Pi.
```

The run covers `0 <= t <= 2.4`. Both spherical boundaries impose this exact
solution. The pulse is localized along its propagation direction and uniform
transversely; it is not a compact three-dimensional wave packet. Unlike the
radial study, propagation has a nonzero normal component across angular
interfaces. For example, the LowerY-to-UpperX face is `x=-y>0`, with unit normal
`(1,1,0)/sqrt(2)`, whose dot product with n is positive. A peak characteristic
`x(t)=(1.44,-0.72,0)+(t-2.8)n` crosses that face near `t=1.90` at radius about
1.84, inside the shell and away from the excision boundary.

The initial angular/angular/radial extents are `[10,10,14]`, `[18,18,14]`,
`[14,14,14]`, `[18,18,14]`, `[14,14,14]` and `[16,16,14]`, in Sphere's
UpperZ, LowerZ, UpperY, LowerY, UpperX, LowerX block order. Independent
interpolation checks give a combined initial physical L2 error of `1.48e-7`
with 19,544 points. These extents were selected before evolution, with the
same initial mesh for all adaptive strategies. Offline trials are retained
in `build-agents/chalmers-wedges/results/*ResolutionTrials.json`.

The single shell was chosen to isolate angular transport and keep the
experiments small enough for local controls. Preliminary interpolation trials
also considered the previous three-shell domain and broader pulses. No
preliminary trial was evolved or included in the cost comparisons. Results
from the radial and oblique studies concern different problems and should not
be compared as a common accuracy-versus-cost sweep.

Time integration, AMR cadence, refinement limits and the criterion are
unchanged: global Dormand-Prince 5, CFL 0.5, AMR every 0.04, degree cap 20,
h-level cap 3, no filtering and no coarsening. Seven volume observations
are saved, every 0.4. The native executable is the same as in the radial
study. Only the experiment setup and Python diagnostics were extended.

The volume analysis reuses the exact serialized maps and the derivative
constraint `partial_i Psi - Phi_i`, with independent physical 32^3/40^3
quadrature. The finest fixed-grid run exposed a `9.10e-12` difference with
24^3/32^3 quadrature that exceeded the convergence check. The saved fields were
reanalyzed at the higher orders; no native evolution was repeated for this
diagnostic change, and the original analyses and failure logs are retained.
It additionally measures both traces on the full LowerY/UpperX angular face:

- Physical surface L2 error of all five fields on each side, and the L2 norm
  of their jump.
- Physical energy flux `Pi * Phi.n` from each reconstructed trace, compared
  with the exact solution. This checks transport across the chosen face; it
  is not a conservation test of the DG numerical flux.
- Independent 24^2/32^2 quadrature on every rectangle in the common partition
  of the two adaptive faces. This resolves transverse h boundaries explicitly;
  a single quadrature grid over the whole face would not reliably integrate
  jumps between neighbouring elements.

The surface parametrization is verified against both saved block maps. Unit
checks cover different angular/radial partitions on the two sides, unequal
polynomial extents, the exact face area, and a deliberately discontinuous
trace on a radial half-face. A separate finite-difference check verifies the
analytic wave equation. The shorter radial and cube tests remain in place;
`ChalmersWedges.yaml` adds a native oblique-interface regression.

Reproduce a representative run with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersWedges run \
  --executable build-agents/bin/EvolveScalarWave3D \
  --output build-agents/chalmers-wedges/reproduced-hp05 \
  --strategy hp05 --tolerance 1e-6 --cfl 0.5
```

Use `--strategy p`, `hp1` or `fixed --fixed-p 17` for controls, and
`--cfl 0.25` for time-step comparisons. `analyze`, `compare` and `summarize`
regenerate diagnostics from saved data. Each run retains its input, log,
executable hash, volume data, mesh history, directional fits and interface
measurements. The local sweep driver is
`build-agents/chalmers-wedges/run_sweep.py`.


### Oblique-interface results

All **23 configurations complete**: twelve adaptive runs at local targets
`1e-3` through `1e-6`, seven fixed grids (degrees 9, 13, 14, 15, 16, 17, 20),
and four halved-CFL controls. Two additional attempts terminated with signal
11 during startup, before any evolution output; both succeeded on retry.
Their logs and separate failure records are retained. The cause of those
startup failures has not been diagnosed. There are no incomplete evolutions,
fit failures or refinement-limit failures among the successful launches.

Tightening the target decreases the maximum observed field error, derivative
constraint norm, interface trace error and interface jump for every adaptive
strategy. All h refinement is angular; the radial degree remains 13 in every
adaptive run. At target `1e-6`, sigma 0.5 makes one angular split, giving seven
elements, while sigma 1 reaches fifteen elements. Both reach angular h level
1, so the difference is how widely they split, not the deepest split.

At local target `1e-6`, with fixed degree 17 as a control:

<table>
<tr>
<th>
Strategy
</th>
<th>
Final points
</th>
<th>
Maximum observed field L2 error
</th>
<th>
Maximum constraint L2
</th>
<th>
Maximum interface jump L2
</th>
<th>
Evolution seconds
</th>
</tr>
<tr>
<td>
p only
</td>
<td>
25,144
</td>
<td>
2.24e-6
</td>
<td>
1.04e-5
</td>
<td>
2.05e-6
</td>
<td>
18.4
</td>
</tr>
<tr>
<td>
hp, sigma 0.5
</td>
<td>
26,124
</td>
<td>
2.14e-6
</td>
<td>
1.03e-5
</td>
<td>
2.03e-6
</td>
<td>
19.7
</td>
</tr>
<tr>
<td>
hp, sigma 1
</td>
<td>
44,744
</td>
<td>
1.61e-6
</td>
<td>
6.87e-6
</td>
<td>
1.91e-6
</td>
<td>
56.7
</td>
</tr>
<tr>
<td>
Fixed degree 17
</td>
<td>
34,992
</td>
<td>
4.34e-6
</td>
<td>
1.36e-5
</td>
<td>
9.96e-8
</td>
<td>
34.9
</td>
</tr>
</table>

Sigma 0.5 uses **42% fewer points and 65% less time** than sigma 1, with a
somewhat larger field error. P-only remains competitive: it is slightly
cheaper than sigma 0.5 at similar accuracy. At the three looser targets,
sigma 0.5 and p-only produce identical recorded meshes and fields.

Among the tested controls meeting a physical field-error requirement of
`1e-5`, p-only and sigma 0.5 beat fixed degree 17 in both error and cost.
They also beat the qualifying fixed degree 15 control for a `1e-4` requirement.
At `1e-3`, fixed degree 13 is cheaper. These comparisons concern the bulk
field norm: fixed degree 17 has approximately twenty times smaller interface
jumps than the tight p-only/sigma-0.5 runs. The criterion therefore does not
provide equally strong control of every diagnostic. Its local tolerance is
also not a global physical error bound; at target `1e-5`, the maximum field
errors are about `1.5e-5` to `1.8e-5`.

The four CFL checks repeat the adaptive strategies at target `1e-5` and fixed
degree 17 with CFL 0.25. Every recorded h/p mesh is unchanged. Final full-state
L2 differences are `1.63e-9` for p-only and sigma 0.5, `4.30e-11` for sigma 1,
and `9.15e-10` for fixed degree 17. Each is less than **0.022% of the respective
final spatial error**. This supports spatial-error dominance for the checked
cases; the `1e-6` adaptive runs were not separately repeated at half CFL.

Across all runs, the independent quadrature changes are at most `5.21e-17`
for volume norms and `1.56e-15` for interface measurements. Mapped coordinates
agree to `4.45e-16`, shell volumes to `7.11e-15`, and the analytic initial
state agrees with native nodal fields to `1.00e-15`. All maxima above refer to
the seven saved observations, not every time step.

Timings are single-PE local measurements, without repeated timing trials.
The fixed degree-14 evolution briefly overlapped a CMake configuration; its
accuracy and point counts are retained, but its runtime is excluded from cost
curves and rankings. `MatchedAccuracy.csv` selects only tested runs satisfying
each physical-error requirement, without interpolation. Offline resolution
selection, postprocessing and failed startup attempts are excluded from the
reported successful-evolution timings.

`results/AccuracyCost.png`, `DirectionalRefinement.png` and
`InterfaceTransport.png` show the cost comparison, angular adaptation and
transport diagnostics. The accompanying `Summary.json`, `RunSummary.csv`,
`CflComparisons.json`, `Validation.json` and `BuildMetadata.json` retain the
results and provenance. `run_fixed_fill.py`, `reanalyze.py` and `finalize.py`
under the study directory supplement `run_sweep.py` with fixed-grid controls,
quadrature reanalysis and final comparisons.

All **nine targeted CTest entries pass**: four scalar-wave input parse/execute
pairs and the Python entry containing 22 test methods. Black, isort and
`git diff --check` pass. The new native wedge regression takes about 28 seconds
locally. The C++ criterion and executable remain unchanged from the radial
study. At this stage generalized harmonic evolution and black holes had not
been tested; the following study checks their AMR/projector support and
establishes a stationary single-hole control.


## Stationary single-hole GH control

`ChalmersSingleHole.py` applies the criterion to a mass-one, nonspinning
Schwarzschild hole in Cartesian Kerr-Schild coordinates. The criterion is
registered in `EvolveGhSingleBlackHole`; the binary-hole executable is not
changed. The existing single-hole AMR machinery supplies:

- Domain, evolved-variable, mortar and time-step-history projectors for both
  p changes and h splits. Children are registered after projection.
- Filter projection and initialization/copying of the remaining element data,
  including control-measurement and interpolation state.
- `EvaluateAmrCriteria`, `AdjustDomain` and `CheckDomain` phase support, and
  AMR callback registration in the executable. Joint hp refinement uses these
  phases; the `RefineMesh` event only evaluates p-type criteria.
- Apparent-horizon observations that receive the current element IDs, meshes
  and fields after refinement.

The relevant implementations are
`Evolution/Initialization/{DgDomain,Evolution}.hpp`,
`Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp`,
`ParallelAlgorithms/Amr/Actions/InitializeChild.hpp` and
`Evolution/Executables/GeneralizedHarmonic/EvolveGhSingleBlackHole.hpp`
under `src/`. Conservative LTS steppers requiring past-step initialization
explicitly reject h refinement. Variable-order conservative LTS also requires
error-controlled restarting. These experiments use global Dormand-Prince 5,
so they do not validate either LTS path.

The static domain has twelve equiangular wedges in two shells,
`1.8 <= r <= 3` and `3 <= r <= 6`, with linear radial maps. Excision uses
`DemandOutgoingCharSpeeds`, which checks that no characteristic enters the
domain; it supplies no analytic inner-boundary data. The outer boundary and
gauge use the analytic Kerr-Schild solution. Control systems, moving maps,
filtering and coarsening are disabled. The horizon finder remains active.
The damping functions are retained from the repository's Kerr-Schild input.

The main controls evolve to `t=2M`, with AMR decisions every `0.1M` and five
volume/horizon observations every `0.5M`. All adaptive runs start with eight
points in each direction on every wedge, and monitor all components of
`SpacetimeMetric`, `Pi` and `Phi`. A separate `split` strategy uses
`DriveToTargetRefinementLevels` to force one angular and one radial split.
This tests projection mechanics independently of Chalmers' choices; it is
not another error-based strategy.

An important initialization detail differs from the scalar-wave experiments:
GH's `SetPiAndPhiFromConstraints` resets Pi and Phi to satisfy the discrete
initial gauge and three-index constraints. The metric agrees with analytic
nodal data to roundoff, while the reset fields have truncation error. The
diagnostics therefore check initial Phi against the numerical derivative of
the metric interpolant, and measure the full initial error separately.
Later refinement projects the evolved state; it does not replace it with
the analytic solution or recover information absent from the initial grid.

The independent analytic implementation uses
`g_ab = eta_ab + 2/r l_a l_b`, `l_a=(1,x_i/r)`, its exact spatial derivative
for Phi, and `Pi_ab = beta^i Phi_iab / alpha`. The analyzer evaluates the
saved numerical polynomials at independent 20^3/28^3 Gauss points using the
serialized native maps. It integrates field errors and
`partial_i g_ab - Phi_iab` with coordinate volume `d^3x`, counting both
off-diagonal tensor components. These are positive coordinate-component
norms, not invariant spacetime norms or proper-volume norms. Native gauge,
three-index and four-index constraint RMS measurements are saved separately;
their normalization and quadrature differ from the independent integral norm.
Horizon mass, relative area error and dimensionless spin are checked against
1, `16*pi` and 0.

For example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersSingleHole run \
  --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/chalmers-gh/reproduced-hp05 \
  --strategy hp05 --points 8 --tolerance 1e-4
```

Use `fixed --points 10`, `p`, `hp1` or `split` as controls. `--cfl 0.25`
provides the time-step check. The `analyze` subcommand regenerates diagnostics;
`compare` checks final-state differences and saved mesh histories for CFL
pairs with identical final meshes. The short regression input is
`tests/InputFiles/GeneralizedHarmonic/ChalmersSingleHole.yaml`.

The new runner defaults to six Charm++ processing elements (`--processes 6`),
with BLAS/OpenMP threading held at one. The initial short comparison batch
uses one processing element throughout. Core counts are recorded separately
so those timings are not mixed with the six-core follow-ups.


### Single-hole results and a worker-stack fix

The short matrix contains ten evolutions to `t=2M`: fixed degrees 7, 9 and 11,
Chalmers p-only and both hp thresholds at target `1e-4`, the forced angular/
radial split, and three halved-CFL controls. A separate fixed-grid pilot runs
to `0.5M`. All eleven complete. The common adaptive initial field error is
`1.77e-3`, already larger than the local criterion target.

<table>
<tr>
<th>
Strategy
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
Field L2 error at 2M
</th>
<th>
Maximum independent three-index constraint
</th>
<th>
Maximum horizon mass error
</th>
</tr>
<tr>
<td>
Fixed degree 7
</td>
<td>
12
</td>
<td>
6,144
</td>
<td>
2.32e-3
</td>
<td>
5.83e-3
</td>
<td>
4.68e-6
</td>
</tr>
<tr>
<td>
Fixed degree 9
</td>
<td>
12
</td>
<td>
12,000
</td>
<td>
9.39e-5
</td>
<td>
4.52e-4
</td>
<td>
2.34e-7
</td>
</tr>
<tr>
<td>
Fixed degree 11
</td>
<td>
12
</td>
<td>
20,736
</td>
<td>
3.32e-6
</td>
<td>
1.69e-5
</td>
<td>
1.13e-8
</td>
</tr>
<tr>
<td>
p only
</td>
<td>
12
</td>
<td>
7,776
</td>
<td>
1.21e-3
</td>
<td>
3.10e-3
</td>
<td>
5.74e-6
</td>
</tr>
<tr>
<td>
hp, sigma 0.5
</td>
<td>
12
</td>
<td>
7,776
</td>
<td>
1.21e-3
</td>
<td>
3.10e-3
</td>
<td>
5.74e-6
</td>
</tr>
<tr>
<td>
hp, sigma 1
</td>
<td>
48
</td>
<td>
24,576
</td>
<td>
1.32e-3
</td>
<td>
2.59e-3
</td>
<td>
7.16e-6
</td>
</tr>
<tr>
<td>
Forced angular/radial split
</td>
<td>
48
</td>
<td>
24,576
</td>
<td>
1.78e-3
</td>
<td>
3.35e-3
</td>
<td>
3.09e-6
</td>
</tr>
</table>

P-only and sigma 0.5 have identical saved meshes and fields. They increase
the angular degrees from 7 to 8; sigma 1 instead splits both angular
directions, retaining degree 7. No Chalmers run changes radial resolution.
The lower threshold uses 68% fewer final points than sigma 1 and has slightly
smaller field error, while sigma 1 has smaller gauge and three-index
constraints. The higher fixed resolutions improve every listed diagnostic
substantially. This is a projection/evolution control, not a demonstration
of hp efficiency at matched initial accuracy.

For sigma 0.5, sigma 1 and fixed degree 11, halving CFL leaves every saved mesh
unchanged. Final full-state L2 differences are respectively `3.34e-12`,
`2.86e-12` and `2.66e-13`, at most `8.03e-8` of the respective final spatial
error. Independent quadrature changes across the short matrix are below
`8.02e-16`. The initial metric agrees with analytic nodal data to `4.45e-16`;
initial Phi agrees with the numerical metric derivative to `7.52e-14`.

The six-core extension exposed a real implementation issue. Two sigma-1
attempts terminated with SIGBUS at the first h refinement, near `t=0.1M`.
Both macOS crash reports identify a worker-stack guard violation in the
history-projection call chain. Disassembly showed that
`ProjectTimeStepperHistory::apply` reserved approximately 150 KiB for a
temporary `History`, even on the global-stepping branch that does not use
it. Combined with the large AMR message frame, this exhausted the worker
stack. The fixed-grid and p-refining six-core runs did not hit that path.

`src/Evolution/Initialization/Evolution.hpp` now resets the history in place
by clearing its records and allocation caches and setting its integration
order. This removes the large temporary, without changing the projection
formula or requiring a modified Charm++ installation. The split/join unit
tests now also seed old steps and substeps and check that resetting clears
both. The previously failing sigma-1 evolution completes to `10M` on six
cores after this change. Failed runs and compact crash reports are retained;
they are not classified as startup failures or included as successful runs.

All three six-core controls have now completed through `10M` with the same
updated executable. Their common initial field L2 error is `1.77e-3`.

<table>
<tr>
<th>
Strategy
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
Field L2 error at 10M
</th>
<th>
Maximum independent three-index constraint
</th>
<th>
Maximum horizon mass error
</th>
</tr>
<tr>
<td>
Fixed degree 7
</td>
<td>
12
</td>
<td>
6,144
</td>
<td>
2.00e-3
</td>
<td>
6.52e-3
</td>
<td>
7.09e-6
</td>
</tr>
<tr>
<td>
hp, sigma 0.5
</td>
<td>
12
</td>
<td>
7,776
</td>
<td>
1.01e-3
</td>
<td>
4.01e-3
</td>
<td>
8.12e-6
</td>
</tr>
<tr>
<td>
hp, sigma 1
</td>
<td>
48
</td>
<td>
24,576
</td>
<td>
9.55e-4
</td>
<td>
3.38e-3
</td>
<td>
5.58e-6
</td>
</tr>
</table>

The two thresholds retain the same meshes as at `2M`. Sigma 0.5 uses 68%
fewer points, with a final field error about 6% larger than sigma 1; sigma 1
has smaller gauge and three-index constraints. These are short stationary
controls with an analytic outer boundary, not evidence of long-term or
binary-hole stability. Maxima are over the saved observations, every `2M`
in this extension, so they need not capture peaks visible in the denser
short-run output. The sigma-1 full state at `2M` is identical between the
original one-core run and the six-core run with the stack fix, evaluated
independently on the common mesh. Their complete saved mesh histories
differ because the observation cadences differ.

The study contains 16 successful evolutions (11 short controls/pilot, two
six-core controls before the stack fix, and three after it) and the two
diagnosed failures. Timings and executable hashes are retained per run;
the fixed degree-9 short-run timing is excluded from comparisons because
it briefly overlapped pilot diagnostics. No parallel-scaling claim is made.
After rebuilding the GH and scalar-wave executables and the initialization
tests, all 15 targeted CTest entries pass, including the Python entry with
23 test methods, the split/join history tests, and the new single-hole input.

Results are in `build-agents/chalmers-gh/`. The `results/StationaryControl.png`
and `results/SixCoreControl.png` figures show the short and six-core controls;
PDF versions are also saved. `Summary.json`, `RunSummary.csv`,
`CflComparisons.json`, `PostFixAt2M.json`, `Validation.json` and
`BuildMetadata.json` in `results/` preserve comparisons and provenance.
`CrashSummary.json` records the two stack failures. Each successful run
retains its input, native output, mesh history, field errors, constraints
and horizon diagnostics. Reproduce a six-core extension with the command
above plus `--end-time 10 --observation-interval 2 --processes 6`.

The next numerical step is to select sufficient resolution before evolution
and initialize the analytic GH data on that grid, then compare p-only and
both hp thresholds at controlled initial accuracy over a longer stationary
evolution. The large improvement from fixed degree 9 and 11 shows why this
should precede spinning or moving holes: projection during evolution cannot
recover information missing from the coarse initial data. The criterion's
local modal target is also not a guaranteed bound on the global field error.

## Controlled initial accuracy and p-coarsening

`Chalmers` now requires a `CoarseningFactor` option. Zero reproduces the
refinement-only experiments above. A value in `(0, 1)` enables p-coarsening
when the AMR policy also has `AllowCoarsening: true`. Refinement still uses
the original target and the sigma rule; coarsening uses a lower target to
provide hysteresis. `PRefineOnly` selects p instead of h for increases and
does not disable p decreases.

For a proposed reduction from degree `p` to `p-1` in direction `d`, form the
transverse-energy spectrum `S_k`. The candidate indicator is

```
E_coarse,d = sqrt(2/(2p+1)) S_p + E_Chalmers(S_0, ..., S_(p-1)).
```

The first term measures exactly the reference-L2 content removed by the
native modal-truncation projector. The second tests the retained spectrum,
including its new last mode and extrapolated tail. Requiring
`E_coarse,d <= CoarseningFactor * target / Dim` for every monitored component
budgets concurrent directional reductions conservatively. Refinement wins
if any component requests it, and any component can veto coarsening. A
negligible directional variation uses its full sampled norm instead of a
fit. Nondecaying or insufficient candidate fits veto coarsening.

The reducer preserves at least six points per direction, leaving three
modes of each parity available to fit when a quiet direction becomes active
again. Existing refinement limits still apply; configure a minimum no
higher than six, or use the existing limit-clamping policy if a larger
floor is desired. The tests retain `ErrorBeyondLimits: true` so unresolved
upper-limit requests remain visible. The estimator is an experimental local
indicator, not a global or derivative-error guarantee. In particular,
coarsening can increase GH constraint error even if the field indicator
accepts it.

Automatic h merging was not included in this p-coarsening study; the guarded
implementation and its validation are described in the next section. A small
child
spectrum does not certify a common parent: even two constants can disagree
at their interface. Existing joining machinery can enforce sibling agreement
on flags, but a reliable error decision also needs the sibling fields and a
check of the prospective parent projection. P reductions already release
resolution inside both unsplit and previously split elements.

The single-hole runner accepts `--coarsening-factor` (default zero), and its
`select-initial` command measures candidate initial grids after the native GH
initializer has reset Pi/Phi. It selects the first tested grid meeting a
specified full-state coordinate-volume L2 error at `t=0`, then comparison
runs independently initialize analytic data on that same grid. This is an
experiment driver, not a new production initial-AMR phase.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersSingleHole \
  select-initial --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/chalmers-initial-selection --target 1e-5 \
  --candidates 8 10 12 14 --processes 6

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersSingleHole run \
  --executable build-agents/bin/EvolveGhSingleBlackHole \
  --output build-agents/chalmers-coarsening-example --strategy hp05 \
  --points 12 --tolerance 1e-5 --coarsening-factor 0.1 --filter hypercube \
  --end-time 50 --observation-interval 5 --processes 6
```

The new `ChalmersSingleHoleCoarsening.yaml` regression starts at twelve
points and exercises native p reduction, including evolved variables and
RK-history projection. Unit tests cover discarded-mode loss, transverse
Legendre normalization, hysteresis, component vetoes, invalid options,
serialization, and the unchanged h/p refinement choices. The experiment
records AMR statistics every `0.1M`, independently of volume/horizon output,
so brief resolution changes are not hidden by the longer observation cadence.

### Travelling-pulse reversibility check

The first dynamic check uses the same width-0.3 Gaussian plane wave and
rotated one-dimensional intervals as the earlier scalar-wave experiment.
Both controls now initialize every element with 16 points, evolve to `t=6`
on six cores, and use the same local target `1e-5`. One disables coarsening;
the other uses factor 0.1. This isolates p changes before combining them
with h splits. AMR decisions occur every 0.02 and fields are saved every 0.2.

The refinement-only run ends with 265 points. P-coarsening ends with 128,
with a minimum of 120 during the evolution. Nine elements exhibit both
increasing and decreasing degree in their saved histories (45 observed
increases and 59 decreases across all elements). The common initial field
L2 error is `5.08e-10`; maximum observed errors are `4.74e-6` without
coarsening and `8.79e-6` with it. This is an accuracy/cost tradeoff, not an
equal-error comparison or a guarantee that a local target bounds global
error. The mesh map shows resolution following the pulse and being released
behind it.

Halving the fixed step from 0.0005 to 0.00025 produces the same saved meshes
and degree-change counts. Observation times differ by at most `2.20e-12`
from floating-point accumulation; mesh matching uses tolerance `1e-10`.
The maximum nodal full-state difference is `4.88e-11`, and independent
64/128-point quadrature differences are below `5.23e-16`. These comparisons
are recorded in `results/PulseResults.json` and
`results/PulseTimeStepComparison.json` under
`build-agents/chalmers-gh-coarsening/`. `TravellingPulse.png` shows cost and
error histories, and `PulseMesh.png` shows the moving degree distribution.

### Longer single-hole controls

Native initial-resolution probes give full-state initial errors `1.77e-3`,
`7.87e-5` and `4.79e-6` for 8, 10 and 12 points per direction. Thus the
initial error target `1e-5` selects 12 points. The comparison starts all
strategies from that same analytic initialization, with 20,736 points.

The first unfiltered fixed-grid attempt at `50M` fails near `27.98M` with a
non-finite CFL time-step estimate. Field and gauge-constraint errors grow
rapidly after about `15M`, before that assertion. Its last volume observation
at `25M` has field L2 error 1.38 and three-index constraint L2 error 1.16.
This run never adapts, so the failure cannot be caused by coarsening. Its
native log and partial diagnostics remain in `fixed-50M/`, explicitly marked
with nonzero return code and `incomplete: true`. The final time in partial
diagnostics is the last saved observation, not the requested end time.

The subsequent matrix restores the repository Kerr-Schild input's Hypercube
filter (HalfPower 64, enabled on every volume substep). Every strategy uses
this same filter and retains the same initial grid, target `1e-5`, CFL 0.5,
AMR cadence `0.1M` and five-M volume/horizon cadence. The CLI keeps
`--filter none` as its default for reproducibility of the earlier studies;
use `--filter hypercube` for this matrix. Filtering changes the evolved
spectrum and can affect the criterion's choices, so filtered results must
be compared against their equally filtered control. The unfiltered failure
alone does not distinguish spatial, temporal and boundary contributions to
the instability.

Full-state comparisons now integrate on a common block-logical partition
when adaptive meshes differ, using the native curved maps and both
20/28-point quadratures. Symmetric GH tensor components retain their full
multiplicities. A synthetic 50-component test checks this normalization
across different radial partitions. For identical meshes the comparator
subtracts nodal fields first, preserving exact zero differences.

A halved-CFL comparison with the Hypercube filter is a combined discretization
check: filtering is applied on every substep, so halving the time step also
changes how often the filter acts. It must not be described as isolating
pure time-integration error. Older archived input files require adding
`CoarseningFactor: 0.0` when replayed with the new executable; their original
inputs are retained unchanged as experimental records.

The travelling-pulse run can also be reproduced through the tracked CLI:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersEvolution run \
  --executable build-agents/bin/EvolveScalarWave1D \
  --output build-agents/chalmers-pulse-coarsening --strategy p \
  --tolerance 1e-5 --dt 0.0005 --end-time 6 \
  --initial-points 16 --coarsening-factor 0.1 --processes 6
```

Use factor zero for the refinement-only control, and the `analyze` subcommand
for independent diagnostics. The original defaults are preserved.

All four filtered primary runs completed `50M` successfully. Their common
initial field error is `4.79236e-6`; the table gives maxima over saved
observations, not bounds between observations. Field and three-index
constraint norms below are independent coordinate-volume integrals.

<table>
<tr>
<th>
Strategy
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
Max field L2 error
</th>
<th>
Max three-index constraint L2
</th>
<th>
Max horizon mass error
</th>
</tr>
<tr>
<td>
Fixed degree 11
</td>
<td>
12
</td>
<td>
20,736
</td>
<td>
7.585e-6
</td>
<td>
5.609e-5
</td>
<td>
4.629e-9
</td>
</tr>
<tr>
<td>
p only, coarsening factor 0.1
</td>
<td>
12
</td>
<td>
22,176
</td>
<td>
1.100e-5
</td>
<td>
8.826e-5
</td>
<td>
9.266e-8
</td>
</tr>
<tr>
<td>
hp, sigma 0.5, factor 0.1
</td>
<td>
36
</td>
<td>
36,092
</td>
<td>
1.230e-5
</td>
<td>
1.110e-4
</td>
<td>
8.722e-8
</td>
</tr>
<tr>
<td>
hp, sigma 1, factor 0.1
</td>
<td>
36
</td>
<td>
36,092
</td>
<td>
1.230e-5
</td>
<td>
1.110e-4
</td>
<td>
8.722e-8
</td>
</tr>
</table>

Both hp runs initially reduce to 18,144 points, peak at 38,686, and settle
at 36,092. P-only also initially reduces to 18,144, peaks at 22,782, and
settles at 22,176. The AMR statistics contain ten net point-decrease
observations in each adaptive run, versus nine increases for each hp run
and 28 for p-only. These are global net changes, not a count of every
element's decision. The hp runs split in angular directions; p reductions
also release points inside the resulting children. No h joins occur.

At factor 0.1, none of these adaptive single-hole runs improves the
cost--accuracy tradeoff over the equally filtered fixed grid. The pulse
demonstrates useful moving resolution, but that is not sufficient evidence
for a good GH coarsening policy. In particular, a field-based local
coarsening test does not certify the GH constraints, and the fixed grid
already represents this stationary solution well. `StationaryCoarsening`
and `FilteringControl` figures in the study's `results/` directory show
the complete histories. The primary 20/28-point quadrature changes are
below `5.70e-16`.

The same travelling-pulse setup also completes with both hp thresholds and
coarsening enabled. Sigma 0.5 ends with 44 elements and 280 points, with
maximum observed field error `4.65e-6`; sigma 1 ends with 57 elements and
353 points, with maximum error `6.32e-6`. Saved element histories contain
37 and 34 degree decreases, respectively, on children created by h splits.
Thus p reduction works after splitting as well as on the original mesh.
However, the element counts have grown from 16 and cannot fall under this
criterion. The six-point floor leaves a residual cost behind the pulse.
This makes prospective-parent error checks for h merging a concrete next
requirement, rather than assuming p-coarsening alone makes hp reversible.
Full metrics are in `results/HpPulseResults.json`.

### Coarsening strength and remaining limitations

Short hp controls at sigma 0.5 isolate the coarsening setting. All start
from the same twelve-point grid and use the same filter, target, CFL and
observation cadence. At the common time `5M`:

<table>
<tr>
<th>
Coarsening factor
</th>
<th>
Elements
</th>
<th>
Points
</th>
<th>
Field L2 error
</th>
<th>
Three-index constraint L2
</th>
</tr>
<tr>
<td>
0 (refinement only)
</td>
<td>
34
</td>
<td>
58,752
</td>
<td>
4.807e-6
</td>
<td>
2.606e-5
</td>
</tr>
<tr>
<td>
0.01
</td>
<td>
40
</td>
<td>
55,968
</td>
<td>
3.928e-6
</td>
<td>
2.219e-5
</td>
</tr>
<tr>
<td>
0.1
</td>
<td>
36
</td>
<td>
36,272
</td>
<td>
1.179e-5
</td>
<td>
9.360e-5
</td>
</tr>
</table>

Grid growth therefore also occurs without coarsening. Factor 0.01 uses
about 4.7% fewer points than refinement only and has smaller errors in this
short comparison. It also changes the h-refinement trajectory, so these are
different partitions, not a monotone convergence sequence on one mesh.
The factor-0.01 probe completes `10M`, ending with field error `3.669e-6`
and the same 55,968 points. This is a promising conservative setting for
further tests, not a demonstrated long-term optimum. Factor 0.1 releases
many more points relative to refinement only but raises both field and
constraint errors. The separate fixed grid has error `7.568e-6` at `5M`
with only 20,736 points. Thus comparison against a fixed grid at matched
accuracy remains necessary. `results/CoarseningAt5M.json` records these
common-time data; the 50M table above uses final point counts instead.

The filtered CFL 0.5/0.25 pair at `2M` ends on different meshes, with
33/18 elements and 36,772/23,736 points. Their maximum observed field
errors are `1.128e-5` and `1.137e-5`. The full-state difference integrated
on a common partition is `4.119e-6`, about 36% of the finer-step run's
field error, with quadrature change `1.15e-17`. This is a material
dependence on the combined time-step/filter setup. It does not establish
time convergence or identify which contribution dominates. It should be
resolved before extrapolating efficiency to moving holes or binaries.

The two primary hp thresholds have identical saved mesh histories and
exactly identical final nodal fields at `50M`. Their field-difference norm
is zero. P-only versus hp gives a full-state difference `8.814e-6` on the
common partition, with quadrature change `4.57e-17`. Each adaptive run's
saved element IDs and directional degrees remain unchanged from `10M`
through `50M`; the saved cadence cannot exclude changes between those
observations. These checks are in `results/hp1-vs-hp05.json`,
`results/p-vs-hp05.json`, and `results/RecordedMeshChanges.json`.

Validation comprises 19 passing targeted CTest entries, including 23
Python test methods. The final six-point-floor assertions and generalized
50-component comparison checks were rebuilt and re-run successfully.
Inputs, logs, analysis tables, source/executable hashes, and validation
metadata are retained under `build-agents/chalmers-gh-coarsening/`.
The study above tested p coarsening only. The subsequent h-joining work is
described below; a constraint-aware GH coarsening decision remains open.
The branch changes remain uncommitted.

## Guarded sibling joins and pulse exit

`HCoarsening: true` now enables automatic h reduction in the native AMR
phase. It requires positive `CoarseningFactor`, `PRefineOnly: false`, and
anisotropic policies. In multiple dimensions it also requires
`UseDirectionalEnergy: true`. The option defaults to false in the experiment
drivers; include it explicitly in YAML, including when replaying earlier
archived inputs with the new executable.

Before evaluating criteria, each exact face-sibling pair exchanges the nodal
values and meshes for every monitored tensor component. The criterion forms
the proposed parent with SpECTRE's native child-to-parent projection at the
maximum child extent in each direction. For every component it requires

```text
E_join = ||u_children - prolong(project(u_children))||_L2(reference parent)
         + sum_d E_Chalmers,d(project(u_children))
       <= CoarseningFactor * (AbsoluteTarget + RelativeTarget * max|u_parent|).
```

The first term integrates the difference between the represented polynomials
exactly, including a possible jump between the children. Unequal child orders
are evaluated at the larger extents, so parent modes absent on the smaller
child cannot be missed. The second term checks resolution on the prospective
parent, preventing a lossless projection from being mistaken for adequate
resolution. Any component can veto a join, and a local refinement request
takes precedence. This is an extension of the experimental criterion, not a
claim that this entire coarsening prescription was specified by Chalmers.

One logical axis is considered per AMR round, cycling through all dimensions.
Only exact face siblings can join, and a successful join retains polynomial
order. These restrictions make the checked two-child parent identical to the
parent actually constructed; independently accepted joins along several axes
would otherwise certify a different approximation. Existing sibling agreement,
2:1 balancing, and minimum h levels still apply. Nonconforming siblings must
first become an exact pair. The AMR round is inherited across splits and joins,
while temporary sibling data are discarded. Criteria with h coarsening disabled
retain their original communication path. DG nodegroups and retained coarse
grids are not supported by this new exchange.

The guard measures reference-element field error. It does not bound physical
derivative constraints, global accumulated evolution error, or unsampled
structure. Native validation below uses global Dormand-Prince stepping and
static maps; it does not establish LTS, moving-map, or GH h-coarsening accuracy.

The new tests check an interface jump between constant children, a resolved
parent represented with unequal child orders, all three join axes with mixed
transverse modes, and a lossless but insufficiently resolved parent. An
independent Gauss integration verifies projection loss. Action tests cover
both message-arrival orders, repeated exchanges, serialization, a second-field
veto, the minimum h level, rejection of nonconforming or cross-block candidates,
and inheritance of the AMR round. Existing history-projection tests include
unequal child meshes, derivative-only entries, time IDs, and integration order.
`ChalmersCoarsening1D.yaml` provides a short native join/p-reduction regression.

All scientific evolutions in `build-agents/chalmers-h-join/` use six Charm++
PEs with one BLAS/OpenMP thread. Native evolutions run sequentially. Mesh
observations are at every AMR interval (0.02 for the pulse, 0.04 for the shell),
so split/join counts come from parent-child ID transitions, not net element
counts. P-change counts include unchanged element IDs only. All raw inputs,
logs, mesh histories, errors, and failed attempts are retained. Wall times
are not a performance comparison; compilation and postprocessing can overlap.

### One-dimensional pulse exit

The pulse starts at x=-2 with width 0.3 on [-4,8]. All comparisons begin
with 16 elements and 16 points per element, monitor Psi/Pi/Phi, use sigma
threshold 0.5, and evolve to t=12 at dt=0.0005. The minimum h level is three
in each of the two blocks, so 16 elements is the configured floor. Independent
64- and 128-point Gauss integration measures the combined physical L2 error
at all 601 observations. The initial error is `5.078e-10` in every case.

At local target `1e-5`, with coarsening factor 0.1 when enabled:

<table>
<tr>
<th>
Policy
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
Splits / joins
</th>
<th>
Maximum field error
</th>
</tr>
<tr>
<td>
Refinement only
</td>
<td>
16
</td>
<td>
272
</td>
<td>
0 / 0
</td>
<td>
6.601e-6
</td>
</tr>
<tr>
<td>
P coarsening
</td>
<td>
57
</td>
<td>
342
</td>
<td>
41 / 0
</td>
<td>
5.990e-6
</td>
</tr>
<tr>
<td>
H+p coarsening
</td>
<td>
16
</td>
<td>
96
</td>
<td>
28 / 28
</td>
<td>
8.762e-6
</td>
</tr>
</table>

The hp run returns every element to six points after the pulse exits. It
records 105 p increases and 265 p decreases on persistent element IDs, with
no immediate p reversals. Three isolated split-then-join reversals occur at
consecutive AMR observations. Hysteresis therefore does not eliminate every
redundant topology change; no general chatter-free claim is made. Coarsening
can also change the refinement path: low-order elements subsequently split
where the refinement-only control retained higher order.

At the stricter target `1e-8`, h+p coarsening again ends at 16 elements and
96 points, with 65 splits, 65 joins, maximum error `6.282e-9`, and final error
`4.810e-11`. P coarsening alone ends at 59 elements / 354 points with maximum
error `7.912e-9`. The refinement-only run stops near t=1.84 on a request to
exceed polynomial degree 20; it is retained as a failure, not assigned a
final error. This is also the build's maximum polynomial degree.

Halving dt in the strict hp run gives identical meshes at all 601 observations.
The maximum full-state physical L2 difference is `9.660e-11`, with a maximum
nodal difference `1.277e-10`. Differences are measured between fields, not by
subtracting scalar error norms. The two isolated immediate h reversals are
the same at both time steps. A separate short, deliberately overresolved
pilot joins 32 initial elements down to 18 with maximum error `6.743e-10`.

### Curved spherical wedges

The spherical wave uses the previous 18-wedge, three-shell setup on 1<=r<=4,
now evolved to t=4.8. Its center crosses the outer boundary at t=2.6. All
main comparisons start on the same mesh with 49,392 points and initial
physical L2 error `1.234e-7`. The local target is `1e-4`. Every one of the 121
observations is independently integrated at 24^3 and 32^3 Gauss points with
the exact map Jacobian; mapped coordinates and the total shell volume are
also checked.

<table>
<tr>
<th>
Policy
</th>
<th>
CFL
</th>
<th>
Final elements / points
</th>
<th>
Maximum field error
</th>
<th>
Final field error
</th>
<th>
Final constraint norm
</th>
</tr>
<tr>
<td>
Refinement only
</td>
<td>
0.5
</td>
<td>
24 / 61,152
</td>
<td>
2.424e-4
</td>
<td>
2.410e-5
</td>
<td>
1.516e-3
</td>
</tr>
<tr>
<td>
P coarsening, factor 0.1
</td>
<td>
0.5
</td>
<td>
70 / 35,597
</td>
<td>
4.115e-4
</td>
<td>
3.046e-4
</td>
<td>
2.696e-3
</td>
</tr>
<tr>
<td>
H+p coarsening, factor 0.1
</td>
<td>
0.5
</td>
<td>
70 / 35,597
</td>
<td>
4.115e-4
</td>
<td>
3.046e-4
</td>
<td>
2.696e-3
</td>
</tr>
<tr>
<td>
H+p coarsening, factor 0.1
</td>
<td>
0.25
</td>
<td>
63 / 29,105
</td>
<td>
4.073e-4
</td>
<td>
2.339e-4
</td>
<td>
2.649e-3
</td>
</tr>
<tr>
<td>
H+p coarsening, factor 0.01
</td>
<td>
0.5
</td>
<td>
67 / 44,242
</td>
<td>
3.331e-4
</td>
<td>
2.380e-4
</td>
<td>
1.970e-3
</td>
</tr>
</table>

This long exit test does **not** demonstrate satisfactory curved-domain
cleanup. At factor 0.1, only two h joins are accepted, both immediately after
a split; at factor 0.01 and in the half-CFL run, none are accepted. A
substantial
Phi-dominated residual and `partial_i Psi - Phi_i` constraint remain after the
analytic pulse has left. Smaller coarsening factors improve the errors but
do not resolve the issue. Merely forcing joins would violate the guard's
resolution test.

The h-enabled and p-coarsening-only runs end on the same mesh and differ by
just `1.064e-6` in the combined physical L2 field norm, so h joining is not
the main source of their much larger residual. The half-CFL run follows a
different mesh path; its full-state difference from the standard-CFL run is
`1.859e-4`, about 80% of its final field error. There is no filter in these
scalar-wave runs. This establishes sensitivity of the coupled stepping/AMR
trajectory, not a clean temporal convergence result. Comparisons use a common
logical partition and converged curved-volume quadrature.

To test the join mechanism independently of that late residual, a short run
starts with each wedge divided radially once. It makes all 18 expected joins:
36 elements / 98,784 points become 18 elements / 21,000 points by t=0.24,
with 294 directional p reductions, no immediate reversals, and maximum field
error `2.198e-5`. Together with the unit tests and 1D exit runs, this validates
automatic h and p reduction, while leaving curved long-time accuracy open.

The next accuracy study should separate time integration from changing mesh
decisions, then test a derivative/constraint guard for coarsening. A frozen or
replayed AMR schedule would make the time-step comparison interpretable.
The present results do not justify advancing h coarsening to moving holes or
BBHs yet.

### Reproduction and validation

The local `run_pulse.py`, `run_sphere.py`, `mesh_changes.py`, comparison
scripts,
and `report.py` under `build-agents/chalmers-h-join/` preserve the exact study
setup. For a main spherical run the checked-in driver also supports:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersSphere run \
  --executable build-agents/bin/EvolveScalarWave3D \
  --output build-agents/reproduced-h-join-sphere \
  --strategy hp05 --tolerance 1e-4 --coarsening-factor 0.1 --h-coarsening \
  --end-time 4.8 --observation-interval 0.04 --cfl 0.5 --processes 6
```

Validation finishes with **45 passing CTest entries**: 31 selected unit/Python
entries and 14 input-file parse/execute entries, including the new 1D joining
regression and the existing GH inputs. This does not add a GH h-coarsening
evolution study. Sixteen scientific attempts are archived: 13 successful runs,
the tight refinement-only cap failure, and two Charm++ startup segmentation
faults before evolution. Both startup cases completed on fresh-directory
retries. A newly added unit-test fixture initially requested an unsupported
parent-to-child interpolation to fewer nodes; generating its exact polynomial
on a smaller source mesh fixed the fixture, and the final suite passes.

`results/Summary.json`, `results/CoarseningComparison.png`, and
`results/PulseMesh.png` collect the results. `results/BuildMetadata.json`
records source and executable hashes and the distinction between the initial
binaries and a later rebuild after formatting/help-text changes. Final source
snapshots, formatting checks, and validation logs are retained alongside them.

## Simple fixed-mesh and constraint-guard controls

`ChalmersConstraint.py` isolates two questions using only the flat-space
scalar-wave system: time convergence on a fixed mesh, and a prospective
p-coarsening constraint veto. There is no GH evolution, moving map, filter,
or h joining in this study. Hp refinement remains available in the adaptive
runs. All 12 scientific runs complete on six PEs, with one BLAS/OpenMP thread.
Their inputs, histories, comparisons, and source snapshots are retained in
`build-agents/chalmers-constraint/`.

The 1D pulse retains the previous [-4,8] domain, 16 initial elements, 16 points
per element, local field target `1e-5`, and end time 12. The sphere retains the
18 curved wedges on 1<=r<=4, but **all shells now start with [14,14,18]
points**.
This gives a common initial mesh for fixed and adaptive comparisons, with
63,504 points and initial physical L2 error `2.302e-8`. These spherical runs
are not numerically identical to the preceding section's nonuniform initial
orders. Their field target is `1e-4`, end time 4.8, and primary CFL 0.25.
Observations are every 0.2; maxima below are over those saved times.

### Fixed-mesh time convergence

All saved element IDs, orders, and coordinates are identical within each
time-step comparison. Differences are computed between the evolved fields,
not between scalar error norms. The spherical differences use converged
24^3/32^3 physical-volume quadrature at every saved time.

<table>
<tr>
<th>
Fixed-mesh comparison
</th>
<th>
Maximum combined field difference
</th>
<th>
Final field difference
</th>
</tr>
<tr>
<td>
Pulse, dt=0.0005 vs 0.00025
</td>
<td>
9.658e-11
</td>
<td>
3.456e-14
</td>
</tr>
<tr>
<td>
Sphere, CFL 0.5 vs 0.25
</td>
<td>
2.328e-9
</td>
<td>
2.837e-12
</td>
</tr>
<tr>
<td>
Sphere, CFL 0.25 vs 0.125
</td>
<td>
7.195e-11
</td>
<td>
1.602e-13
</td>
</tr>
</table>

The two spherical differences decrease by a factor 32.35, corresponding to
observed order 5.016, consistent with Dormand-Prince 5. Landing steps are
shortened to hit common scheduled times, so their sizes are not all related
by exactly a factor two. The fixed sphere's maximum field error is `4.075e-8`
and final constraint norm `1.549e-7` at all three CFL settings. The fixed
pulse's maximum field error is `2.807e-5`, despite its tiny time-step
difference: spatial error remains there. These controls establish temporal
accuracy on these resolved fixed meshes. They do not prove time convergence
on the much coarser meshes visited by an adaptive run.

### Small scalar-wave guard prototype

The opt-in `ScalarWave::ChalmersWithConstraintGuard` wraps the existing
field criterion. It first obtains the proposed flags. For every simultaneous
p reduction it constructs the discarded modal part, then evaluates

```text
C_i = partial_i Psi - Phi_i
delta C_i = partial_i(discarded Psi) - discarded Phi_i
max_over_old_nodes sqrt(sum_i delta C_i^2) <= ConstraintChangeTarget.
```

The sign is `C_old - C_new`; the norm is sign independent. The old physical
inverse Jacobian converts logical derivatives. Truncating the highest modal
slab in every reduced direction matches the native p projector, including
mixed terms when several directions coarsen together. If the check fails,
all proposed p decreases are vetoed, while refinement requests are retained.
Unchanged or increasing directions do not consume this budget. Existing
constraint violation is not included in the target and is not repaired.

This is deliberately a **per-event change** check, not a guarantee on the
accumulated constraint, its maximum between nodes, or the curl constraint.
It can also veto a projection that would reduce a pre-existing violation,
because it bounds the magnitude of the change rather than only an increase.
Pi does not enter the instantaneous one-index constraint, so its projection
is governed only by the original field criterion. Subsequent DG evolution
and interface errors are not controlled by this veto. Sibling joining is
explicitly rejected by its options (`HCoarsening` must be false), preventing
the p prototype from being mistaken for a validated h guard.

Unit tests use Psi proportional to a seventh Legendre polynomial and
Phi=partial_i Psi. Dropping P7 from Psi leaves degree-six Phi intact, injecting
the known endpoint constraint `28 * amplitude / length_scale`. Tests cover
1D/2D/3D, simultaneous directions and intersecting modal slabs,
inverse-Jacobian scaling, a nonlinear map,
a pre-existing low-mode violation, and a harmless constant. They also check
the actual criterion veto, acceptance at a looser target, factory creation,
serialization, and rejection of unsupported h coarsening. The native
`ChalmersConstraint1D.yaml` input provides a short parse/execute regression.

### Guard results

All adaptive runs use field coarsening factor 0.1. A smaller constraint-change
target is a stricter veto. The reported constraint is the independently
integrated physical L2 norm, not the nodal maximum used by the veto.

<table>
<tr>
<th>
Pulse policy
</th>
<th>
Minimum / final points
</th>
<th>
Maximum field error
</th>
<th>
Final constraint norm
</th>
</tr>
<tr>
<td>
No guard
</td>
<td>
122 / 342
</td>
<td>
5.616e-6
</td>
<td>
4.365e-5
</td>
</tr>
<tr>
<td>
Change target 1e-6
</td>
<td>
128 / 377
</td>
<td>
8.322e-6
</td>
<td>
6.116e-5
</td>
</tr>
<tr>
<td>
Change target 1e-7
</td>
<td>
130 / 366
</td>
<td>
8.613e-6
</td>
<td>
4.948e-5
</td>
</tr>
</table>

<table>
<tr>
<th>
Sphere policy
</th>
<th>
Minimum / final points
</th>
<th>
Maximum field error
</th>
<th>
Final field error
</th>
<th>
Final constraint norm
</th>
</tr>
<tr>
<td>
No guard
</td>
<td>
11,922 / 18,732
</td>
<td>
5.045e-4
</td>
<td>
1.305e-4
</td>
<td>
1.708e-3
</td>
</tr>
<tr>
<td>
Change target 1e-5
</td>
<td>
15,678 / 26,822
</td>
<td>
5.045e-4
</td>
<td>
1.305e-4
</td>
<td>
1.706e-3
</td>
</tr>
<tr>
<td>
Change target 1e-6
</td>
<td>
31,054 / 44,901
</td>
<td>
3.911e-4
</td>
<td>
3.237e-4
</td>
<td>
3.158e-3
</td>
</tr>
<tr>
<td>
Change target 1e-6, CFL 0.125
</td>
<td>
30,768 / 59,066
</td>
<td>
4.391e-4
</td>
<td>
4.367e-4
</td>
<td>
3.225e-3
</td>
</tr>
</table>

The prototype still permits substantial resolution reduction, but it is
**not an accuracy/efficiency improvement** in these tests. The loose spherical
guard first differs from the unguarded saved mesh at t=1.2, after substantial
error has accumulated. Tightening the guard retains more points, changes
the refinement path, and can worsen the final errors. The tighter guard's
CFL pair ends on different meshes, with full-state L2 difference `2.781e-4`
on a common curved partition. Thus the simple veto does not resolve adaptive
time-step sensitivity. It remains an experimental control, not a recommended
default.

A relevant property of this toy system is that `InitializeConstraints`
sets gamma2=0. The volume equations give
`partial_t C_i = -gamma2 C_i`, so the continuum constraint mode is undamped
and stationary here. Discretization, boundary/interface terms, and projection
can create violations that need not leave with the physical pulse. This is
a reason to test constraint damping next, not evidence that damping alone
will fix AMR. Another bounded follow-up is to compare a guard on *increases*
of the actual candidate constraint, and derivative/trace accuracy during
evolution, rather than only the immediate magnitude of a projection change.

Reproduce a guarded sphere run with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersConstraint \
  --model sphere --executable build-agents/bin/EvolveScalarWave3D \
  --output build-agents/reproduced-constraint-guard \
  --tolerance 1e-4 --end-time 4.8 --cfl 0.25 \
  --constraint-target 1e-6 --processes 6
```

Omit `--constraint-target` for the unguarded adaptive control. Use
`--fixed-grid` (without a constraint target) for fixed-mesh controls, and
`--model pulse --tolerance 1e-5 --end-time 12 --dt 0.0005` for the 1D model.
`results/Summary.json`, `results/TemporalOrder.json`,
`results/ConstraintGuard.png`, and `results/TimeAndFieldErrors.png` collect
the findings. `results/Validation.json` records the final test selection,
quadrature checks, formatting, and source/executable provenance.

Final validation passes 23 targeted unit-test entries (the scalar-wave suite,
Chalmers kernel/criterion, Python binding, and time-stepper-history regression)
and 12 scalar-wave Chalmers input parse/execute entries. C++ formatting,
Black, isort, and `git diff --check` pass. Both result figures were visually
inspected. All 12 scientific attempts finish successfully; no failed attempt
or retry is omitted.

## Filtered scalar waves on Kerr-Schild backgrounds

`ChalmersKerr.py` exercises three smooth problems: Schwarzschild with
`PureSphericalHarmonic` Y00 and Y20, and Kerr with spin a/M=0.5 and Y20.
The initial data have Psi=Phi_i=0 and
Pi=exp[-(r-10M)^2/(2M)^2] Y_l0. The background is stationary and analytic;
the curved scalar-wave executable uses gamma1=0 and gamma2 M=1.
The spherical domain has six angular wedges in each of six radial shells,
with boundaries at 1.7, 3, 6, 10, 16, 24, and 40 M. The inner boundary is
excised with `DemandOutgoingCharSpeeds`; the outer boundary uses
`ConstraintPreservingSphericalRadiation`.

All controls filter all evolved fields at every Runge-Kutta substep with
the exponential hypercube filter, HalfPower=256. This suppresses the
highest Legendre mode to roundoff while barely changing the next mode
at the allowed resolutions. `ChalmersIgnoringHighestMode` omits one highest
mode from each directional spectrum before fitting its decay. The omitted
mode becomes part of the unresolved tail; its filtered zero is not evidence
of a terminated polynomial. The normal `Chalmers` criterion remains
available as the control selected by `--include-highest`.

Coarsening also needs a filter-aware estimate. A p decrease discards the
old highest mode, then filtering on the smaller mesh discards the next
highest one. The estimate budgets both losses and the retained-spectrum
tail, and keeps at least six usable modes (seven grid points). For h joining,
it budgets native child-to-parent projection loss, the union of filtered
parent modal slabs, and the parent's directional tail estimates. These are
reference-element field indicators, not rigorous physical-error or
constraint bounds. The one-mode approximation is appropriate only for a
sufficiently steep filter; broader filtering requires a corresponding model
of the damped spectrum.

The curved scalar-wave executable now registers AMR initialization, phases,
events, criteria, and Charm++ refinement callbacks. Its projectors update
the domain, fields, time-stepper history, mortars, filters, and interpolation
point information. `ProjectBackground` reevaluates the analytic spacetime
at the new coordinates after p refinement, splitting, and joining, avoiding
interpolation error in the prescribed metric.

The main experiments in `build-agents/chalmers-kerr` evolve to t=20M with
Dormand-Prince 5, CFL 0.25, AMR every 0.4M, sigma threshold 0.5, absolute
field target 1e-4, and coarsening factor 0.1. Adaptive runs start at twelve
points in each direction and permit both h/p refinement and p coarsening.
Separate initially subdivided runs enable sibling joining and save every
AMR round, so h coarsening can be counted directly. Refinement levels up to
four are allowed. `NumPolynomialModes: [2, 19]` follows the current policy's
degree convention and permits up to twenty points; the criterion enforces
its stricter usable-mode floor.
The main runs save mesh and field observations every 2M, so their reported
minimum and maximum sizes are sampled values. The h-joining controls save
every 0.4M and can resolve individual adaptation rounds. The earlier Kerr
joining smoke test used a maximum refinement level of two; it did not request
a split beyond that limit during its t<=4M interval.

`ChalmersKerrAnalysis.py` independently reconstructs the solution from its
Legendre coefficients and evaluates proper-volume norms using the analytic
Kerr-Schild spatial metric. In particular, the one-index constraint uses
gamma^ij C_i C_j and sqrt(gamma) d^3x, and the curl norm uses
gamma^ik gamma^jl C_ij C_kl / 2. The state norm includes Psi, Pi, and Phi_i.
Gauss quadratures with 24 and 32 points per direction check integration
accuracy. Independent checks compare the analytic metric against the C++
Kerr-Schild solution and the Schwarzschild shell volume against a radial
quadrature. Waveforms are coordinate-sphere Y_l0 coefficients at r=12M.
The native `ObserveNorms` entries selected as `L2Norm` are nodal RMS values;
they are retained in `Norms.csv` but are not the physical norms used in these
comparisons. `PhysicalHistory.csv` contains the independent proper-volume
diagnostics.
Mesh comparisons integrate on the common refinement of both partitions.
There is no exact evolving solution in these tests: differences from fixed
sixteen-point runs retain finite-reference uncertainty.
Additional fourteen-point, half-CFL, and filter-strength controls cover only
the early t<=4M window. They do not establish long-time convergence of the
adaptive mesh history or the t=20M reference, and this pass does not test
late-time power-law tails.

### Three-case results at t=20M

All three main adaptive evolutions finish successfully. They begin with
62,208 points; the fixed sixteen-point reference has 147,456 points.
The following differences use the full state and proper-volume norm,
evaluated on a common partition with the sixteen-point reference:

<table>
<tr>
<th>
Background / initial mode
</th>
<th>
Smallest saved adaptive mesh
</th>
<th>
Final adaptive points
</th>
<th>
Fixed-12 difference
</th>
<th>
Adaptive difference
</th>
</tr>
<tr>
<td>
Schwarzschild Y00
</td>
<td>
28,563
</td>
<td>
70,546
</td>
<td>
0.242933
</td>
<td>
0.032760
</td>
</tr>
<tr>
<td>
Schwarzschild Y20
</td>
<td>
53,178
</td>
<td>
83,374
</td>
<td>
0.247845
</td>
<td>
0.036632
</td>
</tr>
<tr>
<td>
Kerr a/M=0.5, Y20
</td>
<td>
43,856
</td>
<td>
90,492
</td>
<td>
0.248411
</td>
<td>
0.036694
</td>
</tr>
</table>

The adaptive runs' final one-index constraint norms are respectively
0.01008, 0.009659, and 0.009690. The fixed-12 values are 0.09196, 0.10949,
and 0.10991; fixed-16 values are about 0.0177--0.0178. Coarsening nonetheless
creates transient constraint increases: the finest fixed grid has much
smaller constraints early in the evolution. Positive damping does not make
projection harmless or turn the modal target into a derivative-error bound.

Agreement in this volume norm does not imply improvement in every
observable. For example, the Kerr Y20 waveform at r=12M differs from fixed-16
by 4.39e-6 for AMR versus 7.88e-7 for fixed-12. The adaptive Schwarzschild
Y20 run also generates small spurious Y00/Y40 components, about -4.53e-6
and -3.39e-6 at t=20M. These are much smaller than the rotating case's
mixed modes: Kerr AMR gives Y00=5.48265e-4 and Y40=-2.64412e-4,
compared with fixed-16 values 5.52799e-4 and -2.61037e-4. A waveform-focused
accuracy study will need tighter targets or additional indicators.

### Explicit coarsening controls through t=4M

These runs start with 72 elements and 124,416 points (one radial subdivision
per block) and enable h joining as well as p coarsening:

<table>
<tr>
<th>
Background / initial mode
</th>
<th>
Smallest saved / final points
</th>
<th>
H split events
</th>
<th>
H joins
</th>
<th>
P decreases
</th>
</tr>
<tr>
<td>
Schwarzschild Y00
</td>
<td>
45,702 / 64,472
</td>
<td>
45
</td>
<td>
32
</td>
<td>
961
</td>
</tr>
<tr>
<td>
Schwarzschild Y20
</td>
<td>
59,080 / 59,080
</td>
<td>
28
</td>
<td>
48
</td>
<td>
702
</td>
</tr>
<tr>
<td>
Kerr a/M=0.5, Y20
</td>
<td>
56,502 / 64,252
</td>
<td>
42
</td>
<td>
35
</td>
<td>
647
</td>
</tr>
</table>

An h split can involve multiple axes and add more than one element. P counts
refer to changed directions on element IDs present in both adjacent
snapshots; they exclude simultaneous topology changes. The respective
full-state differences from fixed-16 at t=4M are 1.886e-4, 2.095e-4, and
2.299e-4. The different initial mesh makes these coarsening demonstrations,
not controlled wall-time comparisons with the main runs. Immediate h
reversals occur 1, 4, and 6 times, respectively. Their cost and possible
suppression by a dwell time or stronger hysteresis remain to be studied.

### Early convergence and sensitivity controls

At t=4M the full-state differences for resolutions 12/14/16 are:

| Background / initial mode | 12 versus 14 | 14 versus 16 |
| --- | ---: | ---: |
| Schwarzschild Y00 | 6.593e-4 | 5.275e-5 |
| Schwarzschild Y20 | 6.672e-4 | 5.347e-5 |
| Kerr a/M=0.5, Y20 | 6.679e-4 | 5.353e-5 |

The reduction is about a factor of 12.5 in each case. A fixed-12 Schwarzschild
Y00 control with CFL 0.125 instead of 0.25 differs by 2.482e-8 at t=4M,
well below the spatial-resolution differences in this early interval. This
does not check time-step sensitivity along changing adaptive mesh histories.
A fixed-16 Schwarzschild Y00 run with HalfPower=128 instead of 256 differs
by 1.754e-9 at t=4M. Both choices are extremely steep at this resolution;
this is an early fixed-grid sensitivity check, not validation of discarding
only one mode for arbitrary filters or later adaptive histories.

The full t=20M Kerr control with `--include-highest` also completes. Its final
mesh has 257 elements and 80,050 points, compared with 176 elements and
90,492 points for the filter-aware policy. Its full-state difference from
fixed-16 is 0.043108 instead of 0.036694, and its one-index constraint norm is
0.01860 instead of 0.009690. Thus the filter-aware policy improves these
diagnostics while using about 13% more points and fewer elements. It does
not improve every weak multipole: the ordinary policy's final Y00 component
is closer to the reference. No cost optimum or uniform accuracy advantage
is established. This comparison changes the tail interpretation, usable-mode
floor, and coarsening budgets together; it does not isolate only the fit.

Final saved modal coefficients confirm that the steep filter suppresses the
highest mode to roughly 1e-15 of the largest coefficient. Independent final
waveform quadratures (24/48/96 polar nodes and twice as many azimuthal nodes)
are retained in `WaveformQuadrature.json`; extraction uncertainty matters
when interpreting the smallest differences. Proper-volume norm quadrature
checks and every matched-mesh comparison are retained separately.

All 19 completed scientific runs and five failed development attempts are
listed in `results/Summary.json`. `results/KerrScalarWave.png`,
`results/WaveformDifferences.png`, and `results/ModeMixing.png` show the
main histories. `results/Validation.json` and `results/BuildMetadata.json`
record tests, quadrature checks, source snapshots, executable hashes, and
the limitations of the comparison. Native evolutions ran sequentially on
six cores; some builds and analyses overlapped, so recorded timings are
not clean performance benchmarks.

Reproduce one adaptive case with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersKerr \
  --executable build-agents/bin/EvolveCurvedScalarWaveKerrSchild3D \
  --output build-agents/reproduced-kerr-20 \
  --ell 2 --spin 0.5 --end-time 20 --processes 6
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersKerrAnalysis \
  build-agents/reproduced-kerr-20
```

Use `--spin 0 --ell 0` or `--spin 0 --ell 2` for Schwarzschild;
`--strategy fixed --points 16` for a fixed control; and
`--radial-level 1 --h-coarsening --end-time 4 --observation-interval 0.4`
for the short sibling-joining experiment. Every attempt retains its input,
native log, executable hash, and exit status, including failed attempts.

Validation passes 43 selected CTest entries: 22 Chalmers/curved-wave unit
entries, the existing flat-wave constraint-guard unit entry, and 20 input
parse/execute entries. The new curved-wave input exercises native splitting
and joining. Kernel tests cover filtered-tail fits, parity zeros, the
usable-mode floor, the additional filter loss after p coarsening, and exact
parent-filter losses including intersecting slabs in two and three dimensions.
The background projector is checked against direct analytic evaluation on
different grid sizes for Minkowski in 1D/2D/3D and Kerr-Schild in 3D.

The retained development failures exposed three issues before successful
long runs: the polynomial-limit option's degree convention, missing Charm++
refinement callback registration, and an h-refinement cap of two that was
too low as the wave reached the outer shell. These are corrected in the
current configuration and executable. Two pre-evolution startup SIGSEGVs
were retried in separate directories. The existing flat-wave 3D input also
needed a larger test timeout: its completed evolution was just at the old
10-second harness limit while other work was running. Original failing logs
remain alongside final passing tests.

## Kerr tolerance and coarsening study

The follow-up study in `build-agents/chalmers-kerr-convergence` keeps the
same Kerr background (a/M=0.5), Y20 initial data, gamma2=1, HalfPower=256
filter, sigma threshold 0.5, and 0.4M adaptation interval. All new evolutions
run to 20M, save fields at every adaptation time, and use six PEs serially.
Builds and heavy analyses do not overlap these timed evolutions.

The new adaptive runs start with 16 points in each direction. Independent
Gauss-24/32 proper-volume quadrature gives an initial Pi interpolation error
of 3.1492054e-4 at 12 points, versus 5.8532168e-7 at 16 points. Starting the
tighter tolerance sweep at 12 points would introduce an avoidable common
initial error. The two quadratures agree to 1.6e-15 in this check.

The completed matrix is:

| Case | Initial points per direction | Absolute target | Coarsening | CFL |
| --- | ---: | ---: | --- | ---: |
| Fixed reference | 18 | n/a | None | 0.25 |
| Fixed reference | 20 | n/a | None | 0.25 |
| Adaptive (three runs) | 16 | 1e-4, 3e-5, 1e-5 | p | 0.25 |
| Adaptive (three runs) | 16 | 1e-4, 3e-5, 1e-5 | h and p | 0.25 |
| Time-step control | 16 | 3e-5 | h and p | 0.125 |

Both adaptive families permit h and p refinement; only the permission to
join siblings differs. CoarseningFactor remains 0.1. The existing fixed-16
run supplies the lower-resolution member of the reference convergence check,
but its earlier timing is not a controlled performance comparator.

`Manifest.json` records the matrix and executable hash before running.
`run_study.py` retains failed attempts and retries only startup SIGSEGVs
before reduction output exists. `analyze_study.py` runs after all timed
evolutions and checks proper-volume norms, matched-time differences, mesh
changes, and coordinate-sphere multipoles. The local modal target is not a
bound on global solution error or derivative constraints.

For the smaller multipoles, `ChalmersKerrAnalysis --patch-waveforms` integrates
separately over the DG elements intersecting r=12M. It uses Gauss-24/32
quadrature on each spherical patch and checks complete 4pi solid-angle
coverage. This avoids applying a single angular quadrature across DG
interfaces. The measure remains coordinate solid angle, not Kerr proper
area. The extractor checks its assumptions of static spherical wedges,
linear radial maps, and a radial third logical axis.

`ChalmersKerrReversals.py` replays immediate split/join reversals. Observations
follow AMR, before the next evolution step. This is confirmed by the action
ordering and the changed mesh at the first 0.4M observation. For join-to-split
reversals, the replay uses the actual joined parent, applies the native steep
filter, and reconstructs the evolved parent from its new children at the
next observation. Reconstruction uses native projection matrices and is
accepted only if the reconstructed parent reproduces the children. A
simultaneous
one-point p decrease can be inverted to filter/roundoff accuracy by assuming
the old highest mode was suppressed by exp(-36) on the preceding RK substep;
use of this assumption is recorded per event. Larger decreases remain ambiguous.
A split absent from the local requests is consistent with mesh balancing;
neighbor messages are not part of the saved volume data. A reversal requires
overlap between the joined and split axes; opposite changes on different axes
are counted separately. The study does not introduce a cooldown without
evidence from this replay.

### Results of the tolerance sweep

All nine evolutions reached 20M on 2026-10-03. One h-joining attempt at
target 1e-4 had a startup SIGSEGV before any reduction output; an identical
retry completed. Every input, log, exit status, and executable hash is retained.
No C++ implementation changed during this study.

For the runs with both h and p coarsening, the final proper-volume diagnostics
and maximum waveform differences over 0--20M are:

<table>
<tr>
<th>
Absolute target
</th>
<th>
State difference from fixed-20
</th>
<th>
Final one-index constraint
</th>
<th>
Max waveform difference, l=0
</th>
<th>
l=2
</th>
<th>
l=4
</th>
</tr>
<tr>
<td>
1e-4
</td>
<td>
2.69464e-2
</td>
<td>
9.78442e-3
</td>
<td>
4.13547e-6
</td>
<td>
2.32173e-6
</td>
<td>
2.46569e-6
</td>
</tr>
<tr>
<td>
3e-5
</td>
<td>
7.59261e-3
</td>
<td>
1.16938e-3
</td>
<td>
1.28396e-6
</td>
<td>
1.27216e-6
</td>
<td>
7.83961e-7
</td>
</tr>
<tr>
<td>
1e-5
</td>
<td>
3.83085e-3
</td>
<td>
5.90142e-4
</td>
<td>
2.81098e-7
</td>
<td>
2.02739e-7
</td>
<td>
1.64870e-7
</td>
</tr>
</table>

Disabling h joining gives nearly identical error norms and waveforms. The
coarsening comparison is therefore mainly about mesh history and cost here:

<table>
<tr>
<th>
Target
</th>
<th>
H joining
</th>
<th>
Time-averaged points
</th>
<th>
Final points
</th>
<th>
Evolution seconds
</th>
<th>
H joins
</th>
<th>
P decreases
</th>
</tr>
<tr>
<td>
1e-4
</td>
<td>
Off
</td>
<td>
73,514
</td>
<td>
88,512
</td>
<td>
173.2
</td>
<td>
0
</td>
<td>
1,069
</td>
</tr>
<tr>
<td>
1e-4
</td>
<td>
On
</td>
<td>
72,485
</td>
<td>
90,784
</td>
<td>
192.4
</td>
<td>
10
</td>
<td>
1,012
</td>
</tr>
<tr>
<td>
3e-5
</td>
<td>
Off
</td>
<td>
77,952
</td>
<td>
101,982
</td>
<td>
221.2
</td>
<td>
0
</td>
<td>
911
</td>
</tr>
<tr>
<td>
3e-5
</td>
<td>
On
</td>
<td>
76,855
</td>
<td>
102,742
</td>
<td>
233.7
</td>
<td>
5
</td>
<td>
879
</td>
</tr>
<tr>
<td>
1e-5
</td>
<td>
Off
</td>
<td>
92,626
</td>
<td>
122,624
</td>
<td>
307.8
</td>
<td>
0
</td>
<td>
907
</td>
</tr>
<tr>
<td>
1e-5
</td>
<td>
On
</td>
<td>
87,185
</td>
<td>
117,056
</td>
<td>
284.6
</td>
<td>
5
</td>
<td>
841
</td>
</tr>
</table>

P decreases count changed directions on unchanged element IDs, excluding
simultaneous topology changes. Every observed h join locally reduces the
number of grid points. Time-averaged points use the post-AMR mesh on each
following 0.4M interval. H joining saves about 5.9% of average points and
7.5% of runtime at 1e-5; at the looser targets its small point savings do not
offset its overhead. These are single-run timings, including observations,
using the existing -Oz build with debug checks on a shared desktop; they are
not production performance benchmarks. All adaptive runs begin with 147,456
points, coarsen substantially, and subsequently refine as the wave moves.

The reference convergence check gives final state differences 2.90803e-2
for fixed-16 versus fixed-18 and 8.38302e-3 for fixed-18 versus fixed-20.
The latter is larger than the two tightest adaptive comparisons. Thus the
finest global error values should not be treated as errors against a fully
converged solution. A stronger radial reference is needed to establish them.
The waveform reference is much better resolved: over the full history the
18--20 differences are at most 1.263e-12, 3.268e-9, and 9.424e-13 for
l=0,2,4, respectively. The independent patch quadratures agree to about
6e-15; proper-volume norm quadrature changes remain below 1.14e-13.

### Half-CFL check

The 3e-5 h/p-coarsening control with CFL=0.125 completes in 425.4 seconds,
ending with 163 elements and 92,434 points. Its mesh first differs from the
CFL=0.25 run at 1.2M; only the first three of 51 saved meshes are identical.
The final full-state difference between the two adaptive runs is 1.33765e-4,
about 1.8% of their difference from fixed-20. Their final constraint norms
are 1.17117e-3 and 1.16938e-3. Maximum waveform changes are 9.45e-10,
1.72e-9, and 2.51e-8 for l=0,2,4. This is a combined time-step/AMR
sensitivity check, not a temporal convergence-order measurement on one mesh.

### Reversal diagnosis: sensitivity to symmetry zeros

The three full h/p-coarsening runs have 6, 4, and 0 immediate same-axis h
reversals at targets 1e-4, 3e-5, and 1e-5. The half-CFL run has 6. The earlier
initially subdivided Kerr 4M control also has 6. The saved joined parents do
not immediately request refinement when replayed before or after filtering.

The problematic fits mix a substantial, smoothly decaying parity sequence
with an almost absent opposite-parity sequence. Across 26 field/direction
triggers examined in these reversal replays, the weak parity's observed
logical L2 norm is at most 1.18e-9 times the requested target. Fitting the
strong parity separately gives estimated errors below 0.161 times the target.
Nevertheless, the combined fit can have a negative slope parameter and hence
an infinite extrapolated error because tiny symmetry coefficients lie just
above the numerical cutoff. These separate-parity fits are diagnostics only;
the native criterion has not been changed.

A direct sensitivity probe avoids any ambiguity from reconstructing a
pre-split field. On the actual saved parent in block 17 at 4.4M in the 3e-5
run, adding a constant 2e-15 to Phi_z changes the native decision from no
refinement to an h split in direction 1. The fitted error changes from
6.32e-11 to infinity. This perturbation is smaller than the roundoff changes
seen in native projection round trips. `RoundoffPerturbation.json` and
`probe_roundoff.py` retain and reproduce this check.

Reconstructed parent values reproduce their children to roughly 1e-14, but
that precision is comparable to the weak symmetry coefficients. Consequently,
the replay's local-criterion versus balancing classifications are indicative,
not a record of the original native flags. The direct saved-parent probe
independently establishes the fit's numerical fragility. A cooldown would
not address this fitting problem. The next implementation experiment should
make near-zero/parity handling robust, then repeat the 3e-5 joining and CFL
controls. The 1e-5 run is also a useful accuracy check because it showed no
immediate h reversals in this study.

`results/Summary.json` contains the complete diagnostics; `Validation.json`
and `BuildMetadata.json` in that directory preserve quadrature checks,
synthetic harmonic and projection checks, comparisons, source snapshots,
and provenance. The four inspected figures are `MeshAndConstraints.png`,
`WaveformConvergence.png`, `AccuracyCost.png`, and `ReversalSpectrum.png`.

## Separate-parity fitting experiment

`tools/Amr/ChalmersSymmetry.py` compares independent even/odd fits with
adjacent amplitude pairs. These are proposed modifications of this experiment,
not prescriptions attributed to Chalmers et al. The comparison and native
runs are retained in `build-agents/chalmers-parity`.

Both prototypes withstand 100 roundoff perturbations of an odd exponential.
Separate-parity fits also pass 75 known exponential-spectrum cases and bring
all 26 previously diagnosed reversal spectra below their requested targets.
The paired-envelope variant remains above target in two of those spectra.
Pairing uses only complete pairs for the slope, retains an incomplete final
coefficient in the error, and extrapolates a conservative pair-weight envelope.
It is nevertheless an indicator, not a bound on an unknown tail.

A mixed-spectrum control has even coefficients exp(-1.2 k) and odd
coefficients 1e-3 exp(-0.15 k). The combined and paired fits give sigma=0.664
and 0.648, while the separate fit correctly retains the slow sigma=0.15.
At sigma threshold 0.5 and target 1e-5, this changes the h/p choice. Low-order
coefficients alone would not establish the absence of this weak component.
For these reasons only the separate-parity candidate is implemented natively.

The opt-in factory name is `ChalmersParityIgnoringHighestMode` in the
curved-scalar executable. The original `Chalmers` and
`ChalmersIgnoringHighestMode` names preserve their estimator, including for
the control run. The new criterion requires directional energy and applies
one consistent fit to refinement and both coarsening candidates:

- Fit each parity against its original Legendre degrees, requiring three
  significant coefficients. A parity with two absent final modes is treated
  as a terminated polynomial; insufficient active fits veto coarsening and
  retain the existing diagnostic stop if refinement is required.
- Exclude a parity from regression if its resolved nonconstant reference-L2
  norm is at most 0.01 * budget / sqrt(2), but include that norm in the error.
  The two exclusions together consume at most 1% of the supplied budget in
  quadrature. The constant coefficient is omitted only from this norm test.
- Guard the last retained coefficient of each active parity and extrapolate
  its infinite stride-two Legendre tail. Combine the two errors in quadrature
  and use the smaller active sigma. This guards two highest usable modes when
  both parities are active, so it is also a change in the tail indicator.
- Continue accounting for exact removed slabs, join projection loss, and
  parent filtering loss in the coarsening candidates. Omit the intentionally
  filtered highest mode from fits, as in the preceding study.

This bounds omitted *resolved* variation, not unsampled function content or
physical derivatives. The 1% allocation is an experimental fixed choice.
The native saved-parent perturbation probe now makes the same no-refinement
choice for all seven constant offsets between -1e-14 and 1e-14; the old
criterion's 2e-15 offset triggered an h split.

The evolution matrix retains M=1, spin=0.5, Y20 data, gamma2=1, filter
HalfPower=256, initial 16 points per direction, adaptation every 0.4M, and
end time 20M. Five evolutions use ten PEs serially, with no concurrent build
or heavy analysis: separate-parity target 3e-5 with h joining on/off, the
h-joining run at half CFL, h joining at target 1e-5, and the legacy criterion
at target 3e-5 with joining. Earlier fixed-20 data provide a waveform/state
reference, retaining the global-reference accuracy limitation described above.

### Parity evolution results

All five primary runs reached 20M. The 1e-5 case had two pre-evolution
startup SIGSEGVs; the identical third attempt completed. The original failed
attempts remain in separate directories. At target 3e-5:

<table>
<tr>
<th>
Diagnostic
</th>
<th>
Legacy fit, h joining enabled
</th>
<th>
Separate parity, h joining enabled
</th>
</tr>
<tr>
<td>
Final elements
</td>
<td>
176
</td>
<td>
104
</td>
</tr>
<tr>
<td>
Final grid points
</td>
<td>
102,742
</td>
<td>
91,280
</td>
</tr>
<tr>
<td>
Time-averaged grid points
</td>
<td>
76,855
</td>
<td>
73,156
</td>
</tr>
<tr>
<td>
Evolution seconds, ten PEs
</td>
<td>
158.4
</td>
<td>
163.4
</td>
</tr>
<tr>
<td>
Final proper-L2 state difference from fixed-20
</td>
<td>
7.59261e-3
</td>
<td>
7.45532e-3
</td>
</tr>
<tr>
<td>
Final one-index constraint norm
</td>
<td>
1.16938e-3
</td>
<td>
4.38201e-4
</td>
</tr>
<tr>
<td>
Max waveform difference, l=0
</td>
<td>
1.28396e-6
</td>
<td>
4.77753e-7
</td>
</tr>
<tr>
<td>
Max waveform difference, l=2
</td>
<td>
1.27216e-6
</td>
<td>
5.31703e-7
</td>
</tr>
<tr>
<td>
Max waveform difference, l=4
</td>
<td>
7.83961e-7
</td>
<td>
2.81683e-7
</td>
</tr>
<tr>
<td>
Immediate same-axis h reversals
</td>
<td>
4
</td>
<td>
0
</td>
</tr>
<tr>
<td>
H joins
</td>
<td>
5
</td>
<td>
0
</td>
</tr>
<tr>
<td>
P decreases on unchanged element IDs
</td>
<td>
879
</td>
<td>
702
</td>
</tr>
</table>

The parity fit reduces average points by 4.8% and final points by 11.2%,
and improves these waveform differences by factors 2.4--2.8. The final
one-index constraint is smaller by a factor of 2.7. It does not provide a
measured runtime improvement in this control. The coarsening-disabled-in-h
parity case has the same mesh history and zero matched-time state difference
from the joining-enabled case: no join was accepted in either. Its 173.9 s
runtime versus 163.4 s is timing variability, not a demonstrated benefit from
h joining. Each makes 702 p decreases and no immediate p reversals.

At target 1e-5, the parity run ends with 106 elements and 105,956 points,
with 83,364 time-averaged points and 247.7 s evolution time. It makes 680 p
decreases, no h joins, and no immediate h or p reversals. Final state
difference from fixed-20 is 3.82087e-3 and the one-index constraint is
3.72037e-4. Maximum waveform differences for l=0,2,4 are 1.43326e-7,
3.10762e-7, and 8.58952e-8. Tightening the target improves all three over the
3e-5 parity case. Improvements over the legacy method are not universal:
the earlier six-PE legacy 1e-5 run had a smaller l=2 maximum difference,
2.02739e-7. Its timings are not a matched ten-PE comparison.

Halving CFL preserves all 51 saved meshes with the parity criterion. The
base/half-CFL state differences at 4,12,20M are 6.73091e-10, 2.45443e-9,
and 3.04335e-9, respectively. The half-CFL run has the same adaptation counts
and no immediate h or p reversals. This is a timestep/filter sensitivity
check at matching observed meshes, not a measured temporal convergence order.
The filter is applied at every RK substep.

The global-reference limitation remains: fixed-18 versus fixed-20 differs by
8.38302e-3 in the final state norm. The waveform references are much better
resolved, so the waveform improvements are more decisive than the slight
change in global state difference. The separate-parity fit is retained as an
opt-in experimental criterion; no cooldown is added.

### Explicit h-coarsening control and validation

Because none of the primary parity runs accepted an h join, an additional
4M run starts at radial refinement level one (72 elements), with initial
16-point extents, target 3e-5, and the same physics/filter/CFL. It completes
30 joins and 798 p decreases, with no splits or immediate reversals. The
mesh shrinks from 294,912 points to 60,644 and from 72 elements to 42. Every
join reduces point count; their combined immediate reduction is 82,320 points.
The final state difference from fixed-20 at 4M is 3.99446e-5 and the one-index
constraint norm is 7.07653e-5. This demonstrates actual h and p coarsening
with the new fit. Its finer initial radial mesh differs from the primary
matrix, so it is not a matched cost/accuracy comparison with those runs.

Nine distinct affected CTest entries pass: the four Chalmers kernel tests,
the Python test, both criterion tests, and the Kerr input parse/execute
checks. The final kernel/Python rerun additionally checks the h-join
projection-loss guard with parity fitting and prototype input validation.
Independent Python/native tail comparisons cover very slow exponential
decay, both terminal parities, unequal decay rates, and amplitude rescaling
between 1e-100 and 1e100. Native replays agree with the Python candidate on
all 26 saved spectra and 33 symmetry-breaking amplitudes. The 25 matched-time
state comparisons, proper-volume norm and patch waveform quadrature checks,
and three inspected figures are retained with provenance in `results/`.
`HCoarseningControl.json` records the extra short run separately from the
five-run primary summary.

To select the new criterion through the driver, add `--separate-parity`:

```sh
build-agents/bin/python-spectre -m spectre.tools.Amr.ChalmersKerr \
  --executable build-agents/bin/EvolveCurvedScalarWaveKerrSchild3D \
  --output build-agents/parity-example --ell 2 --spin 0.5 \
  --points 16 --tolerance 3e-5 --separate-parity --h-coarsening \
  --end-time 20 --observation-interval 0.4 --processes 10
```

The next physics test should break the exact initial symmetry, for example
with a nonaxisymmetric spherical harmonic, while retaining coarsening and
fixed-grid waveform controls. A stronger radial reference remains necessary
for certifying the global state error of the tightest runs.

## Nonaxisymmetric Kerr check: real Y21 initial data

The next study is in `build-agents/chalmers-nonaxisymmetric/`. It changes the
initial mode to `(ell,m)=(2,1)` and retains the preceding Kerr background,
gamma2=1, HalfPower=256 filter, ignored highest mode, sigma threshold 0.5,
coarsening factor 0.1, initial 16 points per axis, AMR interval 0.4M, and
CFL=0.25. Native source and executable are unchanged. All evolutions use ten
PEs, one run at a time. This breaks axisymmetry, but retains equatorial
reflection parity; it is not a test with all symmetries removed.

`ChalmersKerr.py --ell 2 --m 1` selects the initial data. The native initial
Pi is the radial Gaussian times **Re Y21**, without a sqrt(2) factor. The
extended patch extractor measures `integral Psi conjugate(Y_lm) dOmega` at
coordinate radius 12M, including both real and imaginary parts. It records
(2,1), (4,1), and (6,1), the negative-m conjugate, and forbidden control
modes. In particular, (4,3) can be generated by cubical-grid errors even
when reflection parity is preserved. Existing m=0 diagnostic columns and
default driver behavior are preserved. State comparisons now reject runs
with different m, interpreting older metadata without m as m=0.

The adaptive comparisons run to 20M: legacy and separate-parity fits at
target 3e-5, parity at 1e-5, and parity at 3e-5 with half CFL. A separate
4M control starts with one radial h-refinement. Matching fixed references
use 18 and 20 points per axis. All adaptive runs permit both h and p
coarsening. Element IDs and extents are recorded after every AMR event.

<table>
<tr>
<th>
Adaptive run
</th>
<th>
Final elements
</th>
<th>
Final points
</th>
<th>
h joins
</th>
<th>
p decreases on unchanged IDs
</th>
</tr>
<tr>
<td>
Legacy, 3e-5
</td>
<td>
158
</td>
<td>
87,692
</td>
<td>
4
</td>
<td>
974
</td>
</tr>
<tr>
<td>
Parity, 3e-5
</td>
<td>
112
</td>
<td>
84,608
</td>
<td>
0
</td>
<td>
718
</td>
</tr>
<tr>
<td>
Parity, 1e-5
</td>
<td>
138
</td>
<td>
110,896
</td>
<td>
2
</td>
<td>
678
</td>
</tr>
<tr>
<td>
Parity, 3e-5, half CFL
</td>
<td>
112
</td>
<td>
84,608
</td>
<td>
0
</td>
<td>
718
</td>
</tr>
<tr>
<td>
Initially refined radial control, 4M
</td>
<td>
42
</td>
<td>
59,880
</td>
<td>
30
</td>
<td>
802
</td>
</tr>
</table>

The half-CFL parity run preserves all 51 saved meshes. No immediate h or p
reversals appear in these parity runs. The radial control starts with 72
elements and 294,912 points, accepts 30 joins, and never splits. Thus the
symmetry treatment retains useful coarsening, including two h joins in the
tighter main evolution.

### Accuracy and cost for Y21

All seven runs completed without startup retries or evolution failures. The
table gives maximum **complex** waveform differences from fixed-20 over
0--20M, and proper-volume one-index constraint norms at 20M:

| Diagnostic | Legacy, 3e-5 | Parity, 3e-5 | Parity, 1e-5 |
|---|---:|---:|---:|
| Time-average grid points | 69,267 | 65,099 | 78,233 |
| Evolution seconds, ten PEs | 168.3 | 121.3 | 177.5 |
| Final constraint L2 | 2.29059e-3 | 3.47964e-4 | 1.03024e-4 |
| Maximum constraint L2 | 7.08138e-3 | 7.16190e-3 | 4.05992e-3 |
| Max difference, Psi21 | 6.16157e-7 | 4.91128e-7 | 1.52159e-7 |
| Max difference, Psi41 | 2.16948e-7 | 7.57413e-8 | 2.43938e-8 |
| Max difference, Psi61 | 5.94337e-8 | 1.42701e-8 | 3.14042e-9 |
| Max spurious Psi43 amplitude | 5.81574e-7 | 1.99551e-7 | 6.44796e-8 |
| Immediate same-axis h reversals | 6 | 0 | 0 |

At the same target, parity uses 6.0% fewer time-average points, with a
smaller improvement in final points (3.5%). The observed evolution takes
28% less wall time, but these are single shared-desktop measurements, not
a timing benchmark. The full norm and waveform analysis ran after the
timed evolutions. Only lightweight metadata checks and a subsecond
extraction test ran alongside the fixed references.

The dominant waveform improvement is modest (a factor 1.25 in maximum
complex difference); the weaker l=4,6 modes improve by factors 2.86 and
4.16. Spurious (4,3) leakage falls by a factor 2.91 and decreases further
with the tighter target. Modes forbidden by the remaining reflection
symmetries stay below 3.2e-13 in all runs. The real-field negative-m
conjugacy check agrees exactly at the reported precision.

These improvements are not uniform across diagnostics. The peak
constraint norm near 12M is essentially unchanged at target 3e-5 (about
1.1% higher for parity), despite the much smaller final norm. Maximum
dominant-mode phase difference, evaluated where the reference amplitude
exceeds 1% of its peak, is 2.76e-7 rad for parity versus 2.17e-7 for legacy.
The weak Psi61 signal has peak amplitude 2.20e-7: its maximum error is
still 6.5% of that peak for parity at 3e-5, falling to 1.4% at 1e-5.

Halving CFL leaves all 51 parity meshes unchanged. The proper-volume L2
difference between the final states is 4.70593e-9. Maximum waveform changes
for l=2,4,6 are 1.49329e-11, 4.15378e-12, and 9.73085e-13. These changes
are much smaller than the spatial differences from the fixed reference.
The native reversal replay classifies the six legacy reversals as four
local split-then-join events and two splits requested after intervening
evolution; the parity runs have none.

The radial coarsening control removes 82,320 points in its 30 h joins;
additional p decreases bring the final count to 59,880. Every accepted
join reduces the point count. Its state difference from fixed-20 at 4M is
2.29681e-5, its final constraint norm is 3.46359e-5, and its maximum Psi21
difference over 0--4M is 6.67787e-9. The tighter main run's two joins also
save points (686 in total), without immediate reversals.

### Reference limitations and validation

Fixed-18 and fixed-20 agree in the extracted l=2,4,6 waveforms to maximum
complex differences of 1.63403e-9, 4.23794e-13, and 5.95660e-14. This is
well below the adaptive differences. Their final global state difference,
however, is 5.92505e-3. The corresponding adaptive differences from
fixed-20 are 7.95364e-3 (legacy), 7.80479e-3 (parity 3e-5), and 4.14474e-3
(parity 1e-5). Thus the fixed references do not certify the tight run's
global error; stronger radial references are still needed. At 4M and 12M
the reference differences are much smaller, 2.57357e-7 and 9.74211e-6.

The affected `Unit.LinearOperators.Python.Chalmers` CTest passes. The new
analytic extraction test covers real and imaginary Y21 components,
negative-m signs, normalization, angular and radial element seams, and
unchanged real m=0 outputs. Actual native initial Pi agrees with its
closed-form profile to roundoff; the initial-16 extraction differs from
exp(-1)/2 by 2.08109e-10. Independent 24/32-point quadratures agree within
5.4e-14 for norms and 3.8e-15 for waveforms across all runs. Twenty-two
matched-time state comparisons pass their quadrature and domain-coverage
checks, and the different-m comparison guard rejects Y21 versus Y20.

`results/Summary.json`, `results/HalfCfl.json`, the three inspected PNG
figures, native inputs/logs, and source snapshots retain the results and
provenance. The build uses copied Python modules (`PY_DEV_MODE=OFF`), so
`RefreshPython.cmake` refreshes only the two changed pure-Python modules
with the project's `configure_file(... @ONLY)` rule. No native rebuild
was needed. `run_study.py`, `postprocess.py`, `analyze_study.py`, and
`report.py` reproduce the run matrix and analysis.

This extends the evidence for the opt-in parity criterion to a rotating
nonaxisymmetric scalar wave with both h and p coarsening. It does not yet
remove all reflection symmetries or test dynamical spacetime evolution.
