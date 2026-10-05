\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Geometric NP matching for a BBH worldtube

This experimental branch adds `PhysicalModel: QuadrupoleGeometric` and
`ThirdOrderGeometric` to the curvature-based `WorldtubeTypeD` boundary
condition. It is based on
`worldtube-radial-response-gauge` (`ba2e805ad6`), with per-object
spherical-harmonic shells and their radial-orientation, shape-transition and
distorted-frame fixes. It does not include the older `q8-worldtube-cluster`
affine gauge matcher or Gauss--Bonnet tracking system.

## Build and select the model

Fetch the branch from the `nikwit/spectre` fork, then use your usual cluster
SpECTRE environment and CMake configuration in a new build directory:

```sh
git fetch nikwit codex/geometric-np-bbh
git switch --track nikwit/codex/geometric-np-bbh
cmake --build build --target EvolveGhBinaryBlackHole -j 10
```

The directory `build` above is your configured cluster build. The model uses
GSL, already a SpECTRE dependency; it needs no Python runtime or external data
tape. A complete parse-only input is provided in
`tests/InputFiles/GeneralizedHarmonic/GeometricWorldtubeBinary.yaml`; replace
its placeholder initial-data path and choose numerical parameters before use.
Build from source on the cluster, not with the local macOS executable.

On the object's inner boundary use the following block, with `Mass` set to
that hole's mass in the simulation's units. The example value 0.2 corresponds
to the smaller mass in a nonspinning q=4 binary with total mass one:

```yaml
Interior:
  ExciseWithBoundaryCondition:
    WorldtubeTypeD:
      ConstraintPreservingSector: Bjorhus
      PhysicalSector: Bjorhus
      GaugeSector: SommerfeldAbsorbing
      PhysicalModel: QuadrupoleGeometric
      Mass: 0.2
      MomentRelaxationTime: None
```

`SommerfeldAbsorbing` is a moving-mesh-compatible example, not a newly
validated BBH gauge choice. Preserve an established compatible gauge setting
and moment-relaxation time when comparing physical models. `None` above
selects the instantaneous tidal fit. Changing only `QuadrupoleGeometric` to
`Quadrupole` gives the legacy second-order comparison. Do not change gamma2,
timestep or gauge simultaneously for that comparison. The single-hole
`ReferenceReplay` and `RadialResponse` gauge prototypes require a static mesh
and must not be copied to a moving binary run.

## Third-order matching

Use the following settings to add the documented third-order terms:

```yaml
PhysicalModel: ThirdOrderGeometric
Mass: 0.2                 # mass of this hole, in simulation units
MomentRelaxationTime: None
```

The complete example input selects this model. Replace it with
`QuadrupoleGeometric` for the second-order control, keeping the mass, gauge,
mesh and time-step settings fixed. Both models have the round inertial-sphere
restriction below. The third-order model currently requires forward time
evolution and `MomentRelaxationTime: None`; the older component-wise filter
is not a consistently transported third-order history.

At each RHS evaluation, the third-order model uses the same leading NP
registration, invariant radius/boost, screen eigenmap and polar-transported
dyad. It adds the electric and magnetic octupoles, the induction and near-zone
responses of the dotted quadrupoles, and the cut correction
`delta_t * (Edot + i Bdot)`. The near-zone profiles use the offline horizon
calibration, `e1 = -92/15`, `b1 = -76/15`. Only transverse model components
enter Psi0/Psi4; measured longitudinal channels remain in the registration.

The undotted fit has 24 real components (five electric and five magnetic
quadrupoles, seven of each octupole). It uses a column-scaled real SVD of
measured Psi4. The ten dotted components are supplied from history, not fitted
independently to the same sphere. Rank-deficient designs are rejected.

The causal construction is:

1. Fit quadrupoles and octupoles with the dotted terms temporarily zero.
   Retain these preliminary quadrupoles and the geometric directions at
   full-step RHS states. Differentiating their third-order bias changes the
   model only at fourth order under the slow-tide expansion.
2. Use the current sample and up to four strictly past full-step samples for
   a backward polynomial derivative on the actual, possibly nonuniform times.
   At the first sample the dotted terms are zero (quadrupole plus octupole
   only); the derivative order increases from one to four as history fills.
3. Reconstruct the no-screen flow `W = d/dT + V^A d_A`. Measure the rigid
   angular velocity from `W N`, where `N` is the inferred angular map.
   Subtract the rotation commutator `[Omega,H]` from both tensor derivatives
   before dividing by the model-time clock rate. This distinguishes rotating
   grid labels or auxiliary axes from a physically changing tide.
4. Reconstruct the model-time covector `tau = -u_flat/sqrt(1-2M/r)`.
   Its tangential potential is a small linear least-squares solve in the
   eight nonconstant `l<=2` scalar harmonics, equivalent to the offline
   linear/quadratic potential space. Fix its mean using geometric area
   weights. This gives `delta_t`; the weighted clock rate is `tau(W)` and
   must be positive at every face point.
5. With these dotted moments fixed, refit the 24 undotted components to
   Psi4, predict the incoming slot and transform back to the NR tetrad.
   The cached target enters the same `wminus/2` conversion and Bjorhus
   correction as second order.

The history and target live on the face-owning element and are PUP serialized
for migration and checkpoints made with this executable. Intermediate
integrator stages may evaluate the target but never enter the retained
history. Repeated full-step calls replace the endpoint, and backward-time
self-start resets or angular resolution changes clear the history. A fresh
volume-data start therefore has a short derivative warmup. Avoid AMR for the
first comparison; the implementation does not transfer history between faces.

The state retains the derivative order, clock rate, tilt-fit residual and
outgoing-fit residual/condition for inspection in a checkpoint. Existing
`ObserveWorldtubeMatching` output remains a passive legacy-model calculation;
it does not report the third-order target. Compare evolved curvature from
volume data. The potential residual measures truncation or nonclosure of the
measured time covector; no exact time-coordinate integrability is assumed.

Native tests compare every radial sector to independent offline Python
fixtures, recover held-out Psi0 for a manufactured time-dependent tide, and
check nonuniform causal history, rotating axes/grid labels, repeated calls,
rollback, serialization and the GH face-data integration. This is an
implementation check, not a demonstrated improvement in live BBH accuracy.
Temporal differentiation may amplify noisy fits; verify timestep and angular
resolution sensitivity in a controlled evolution. The underlying Kretschmann
time derivative retains its existing first-order backward estimate, so the
fourth-order moment stencil does not establish fourth-order accuracy of the
whole matching method.

## Required sphere and grid configuration

The worldtube must be one complete spherical-harmonic face. For ObjectB:

```yaml
# Entries to merge into DomainCreator: BinaryCompactObject:
ObjectA:
  UseSphericalHarmonics: false
ObjectB:
  UseSphericalHarmonics: true
InitialRefinement:
  ObjectBShell: 0
InitialGridPoints:
  ObjectBShell: [12, 8]  # radial points, angular L_max
```

These are fragments, not a complete input file. Keep all the other object,
envelope, outer-shell, initial-data and evolution settings. Both objects now
require an explicit `UseSphericalHarmonics` option. For the selected spherical
shell the refinement value is a scalar (radial only), and grid points are the
two entries `[N_r, L_max]`; the other blocks keep their existing layouts.
Angular L_max must be at least six for this BBH domain creator. Avoid AMR for
the first comparison. Keep `UseWorldtube: false` at domain level: that legacy
option controls the scalar-wave worldtube's functions of time and does not
select this GH boundary condition.

The sphere must remain round in **inertial coordinates**, though its center,
orientation and overall radius may change. Translation, rotation and uniform
expansion are compatible. The implementation checks the actual face's radius
and normal covectors before fitting and rejects nonspherical shape or skew
deformations. A straightforward first test uses `ShapeMapB: None` and
`SkewMap: None`. Shape maps on the other black hole remain independent. If a
spherical shell uses a shape/size map, its domain-map transition requires
`TransitionEndsAtCube: false`; this condition alone does not make a
nonspherical shape compatible with the geometric model.

## Outside-horizon worldtubes and tracking

An outside-horizon worldtube removes the smaller hole's apparent horizon
from the evolution domain. The standard BBH pipeline's horizon measurements
and controls therefore cannot simply be reused. In particular, disabling
only `ShapeB` and `SizeB` is insufficient: rotation, expansion, translation,
skew and ShapeA use the shared `BothHorizons` measurements in this executable.

For a short prescribed-motion test, use prescribed functions of time with
coverage of the full test interval, set all horizon-dependent controls to
`None`, and remove observation events that attempt to find/interpolate the
excised horizon. For example, the control entries are:

```yaml
ControlSystems:
  WriteDataToDisk: true
  MeasurementsPerUpdate: 4
  DelayUpdate: true
  Verbosity: Silent
  Expansion: None
  Rotation: None
  Translation: None
  Skew: None
  ShapeA: None
  ShapeB: None
  SizeA: None
  SizeB: None
```

This leaves motion prescribed by the domain maps; it is not a long-inspiral
tracking solution. Maintain a valid excision margin at the other hole during
the short test. A live BBH with Gauss--Bonnet worldtube tracking needs that
tracker ported separately from `q8-worldtube-cluster`; this branch does not
claim that integration. An old binary checkpoint from a different branch is
not a supported executable restart: start from compatible volume initial
data and choose `ElementsAreIdentical` according to the actual new mesh.

## Algorithm and checks

The measured NP scalars still determine the leading frame, invariant radius
and boost. The new part projects the coordinate-sphere tangents into that
frame's screen, constructs its two-metric, and solves a symmetric generalized
Laplace eigenproblem in real scalar harmonics. Its lowest nonconstant triplet,
aligned with the inertial angular labels, defines common angular directions.
Polar transport carries the complex dyad into that angular frame. The tidal
fit uses the map's angular-area weights and replaces the incoming NP slot.
The `wminus/2` conversion and physical Bjorhus factor two are unchanged.

The map uses `min(8, face L_max)` by default, hence an 81 by 81 eigenproblem at
L_max eight, solved at each matching evaluation. `NP_GEOMETRIC_LMAX` may
select another value between two and the face L_max; export it identically
to all ranks if used. The first test should use the default. No shared-file
research diagnostics or single-hole forcing hooks are included. Existing
`ObserveWorldtubeMatching` quadrupole columns remain passive legacy-model
calculations; they must not be interpreted as the geometric model's imposed
target. Compare curvature extracted from the evolved volume data.

The branch adds native tests of round-sphere eigenmodes, constant metric
rescaling, rotated collocation grids, an independent mixed electric/magnetic
tide with held-out Psi0, and translated/rotated/expanded face geometry. It
also includes the BBH domain and mortar consistency tests:

```sh
cmake --build build --target Test_NewmanPenrose Test_GeneralizedHarmonic \
  Test_DomainCreators Test_NumericalDiscontinuousGalerkin -j 10
ctest --test-dir build --output-on-failure \
  -R 'NP.*Tide|Worldtube|BinaryCompactObject|DG.MortarInterpolator'
build/bin/EvolveGhBinaryBlackHole --input-file Bbh.yaml --check-options
```

The single-hole studies established a regular geometric map and an accurate
Bjorhus prescription with donor Psi0. They did not establish a small live
model error: the geometric pulse error was 97.18% in the replay-gauge test,
versus 25.02% for the conditional offline fit. No BBH evolution or BBH
stability claim accompanies this branch. Treat the cluster run as a controlled
experiment, recording the source revision, exact input, finite reach,
constraints and the complex Psi0 error separately.
