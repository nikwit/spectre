\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Geometric NP matching for a BBH worldtube

This experimental branch adds `PhysicalModel: QuadrupoleGeometric` to the
curvature-based `WorldtubeTypeD` boundary condition. It is based on
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
  -R 'NP.GeometricTide|Worldtube|BinaryCompactObject|DG.MortarInterpolator'
build/bin/EvolveGhBinaryBlackHole --input-file Bbh.yaml --check-options
```

The single-hole studies established a regular geometric map and an accurate
Bjorhus prescription with donor Psi0. They did not establish a small live
model error: the geometric pulse error was 97.18% in the replay-gauge test,
versus 25.02% for the conditional offline fit. No BBH evolution or BBH
stability claim accompanies this branch. Treat the cluster run as a controlled
experiment, recording the source revision, exact input, finite reach,
constraints and the complex Psi0 error separately.
