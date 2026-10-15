\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Recording GH fields on a moving sphere {#example_gh_worldtube_recording}

`EvolveGhBinaryBlackHole` provides the `GhWorldtube` interpolation event. It
records `SpacetimeMetric`, `Pi`, `Phi`, their inertial spatial derivatives,
and `Kretschmann` and its inertial spatial gradient on spheres centered on a
selected domain excision center. There are 204 field components plus seven
geometry components per point. Only the center is mapped from grid to
inertial coordinates. The sampling radii are constant inertial
**coordinate** radii, even when the donor uses expansion, shape, size, or
skew maps. The angular axes remain fixed in the inertial frame. The donor's
domain, controls, and evolution are not modified by this event.

Add the target configuration under the existing `InterpolationTargets`
section:

```yaml
InterpolationTargets:
  GhWorldtube:
    LMax: 16
    Center: [*XCoordB, *YOffset, *ZOffset]
    Radius: 0.3
    AngularOrdering: Strahlkorper
  # Keep the other interpolation targets here.
```

This reuses the existing `Sphere` target and options. For `GhWorldtube`,
`Center` specifies the selected excision center in **grid coordinates**. Use
exactly the same values as the domain (the YAML anchors above are available
in the inspiral pipeline). The target maps this center to inertial
coordinates at each observation, while leaving the sampling radii fixed in
inertial coordinates. `Radius` also accepts a list, e.g. `[0.28, 0.30,
0.32]`, to record concentric spheres, stored in increasing radius order.
Choose positive radii entirely within the donor domain throughout the
recording interval. The example radius and resolution are illustrative, not
validated for a particular binary.

Existing binary input files need the new `GhWorldtube` options. Recording
only occurs when its event is triggered; leave the event out to disable
recording. The inspiral pipeline template and cylindrical binary test input
include an unused target configuration.

For a common physical-time cadence with local time stepping, append a dense
trigger to the existing `EventsAndDenseTriggers` list:

```yaml
EventsAndDenseTriggers:
  - Trigger:
      Times:
        EvenlySpaced:
          Interval: 0.01
          Offset: 1617.0
    Events:
      - GhWorldtube
```

This uses the evolution's dense output and spatial interpolation, rather
than forcing the donor to adopt the recording time step. Choose the cadence
and `LMax` by convergence tests. The offset anchors the time sequence; use
the run's usual stopping condition to bound the recording interval. A slab
trigger is also suitable if observations at slab boundaries are sufficient.

For later map replay, remove the `CleanFunctionsOfTime` event from the donor
input and arrange a volume observation at the end of the recording interval.
The existing volume writer saves the retained global functions-of-time
history. The surface output described below does not itself archive those
histories.

## Output

The callback writes to `Observers.SurfaceFileName` with the `.h5` extension,
under `/GhWorldtube.vol`, using the standard HDF5 volume format, version 0.
Each observation value is the absolute simulation time. Writes may arrive
out of chronological order, so readers must sort by observation value.

Each radius is a separate two-dimensional grid named, for example,
`Sphere0`, `Sphere1`, etc. The grid has extents `[LMax + 1, 2 * LMax + 1]`,
spherical-harmonic basis, and Gauss/equiangular quadratures. With
`AngularOrdering: Strahlkorper`, the theta index varies fastest.
`AngularOrdering: Cce` instead makes phi fastest and reverses the extents
and quadratures. There is no modal truncation after interpolation and no
conversion to float32.

Each grid contains:

- `SpacetimeMetric_*`, `Pi_*`, and `Phi_*`: all 50 independent components,
using
  the standard SpECTRE tensor component names and inertial Cartesian frame.
- `deriv(SpacetimeMetric)_*`, `deriv(Pi)_*`, and `deriv(Phi)_*`: 150
components
  of the inertial spatial derivatives, evaluated on the donor volume mesh.
  The metric derivative is computed independently of `Phi`, preserving the
  numerical three-index-constraint information.
- `Kretschmann` and `deriv(Kretschmann)_x/y/z`: the invariant and its
spatial
  gradient used by the worldtube matching models. Here `Kretschmann` means
  `C_abcd C^abcd = 8(E_ij E^ij - B_ij B^ij)`, equal to the Riemann
Kretschmann
  scalar in vacuum.
- `InertialCoordinates_x/y/z`: the actual target positions at that time.
- `InertialCenter_x/y/z` and `InertialRadius`: constant over each grid,
making
  the sampling geometry explicit without reconstructing the donor maps.

All spatial differentiation happens **before interpolation**, on each donor
DG element using its spectral differentiation matrices and the full
logical-to-inertial inverse Jacobian at the observation time. The volume
`deriv(Pi)` and `deriv(Phi)` feed the curvature calculation; the resulting
volume `Kretschmann` is then differentiated on that same volume mesh. Only
after all these operations are complete does the event interpolate to the
sphere. No spatial derivative is obtained by differentiating the recorded
spheres or finite differencing between the three example radii.

No time derivatives are stored. For this fixed-radius sphere with fixed
inertial angular axes, the sampled position is `x(t)=center(t)+offset`. A
time series derivative of K therefore gives `dK/dt = partial_t K +
dot(center)^i partial_i K`. Subtract the center-velocity term to obtain the
derivative at fixed inertial coordinates. The saved functions of time supply
that velocity.

The output callback rejects non-finite tensor values, including missing
points filled by the interpolation framework. Avoid placing sampling points
exactly on donor element interfaces when one-sided DG traces matter: this
recorder uses the interpolation framework's usual choice of the first
received trace at shared points, not a prescribed inner-side trace.

This target records reference fields in H5 format. The `WorldtubeReplay`
reader and boundary condition described below consume this output. It is
distinct from the existing static `FaceReplay` binary format.

## All-reference replay on one sphere

`WorldtubeReplay` is a three-dimensional GH ghost boundary condition
available in the existing BBH executable. Select exactly one radius from the
recorded `GhWorldtube` subfile:

```yaml
ExciseWithBoundaryCondition:
  WorldtubeReplay:
    H5Filename: /path/to/BbhSurfaces.h5
    SubfileName: GhWorldtube
    Radius: 0.30
```

The reader validates the double-precision fields, grid geometry, angular
ordering, constant radius, and uniform observation cadence. It preloads the
50 GH components of the selected sphere into immutable data shared by
boundary objects in each process. Loading happens during option construction
or boundary deserialization, before evolution I/O. It reconstructs each
inertial Cartesian component with spherical harmonics and local six-point
(degree-five) time interpolation. Both recorded angular orderings are
supported. There is no interpolation between radii, no rotation of tensor
components, and no time extrapolation beyond floating-point endpoint
tolerance. Repeated self-start calls evaluate the same time-dependent
prescription; this is interpolated reference data, not a replay of the
donor's discrete multistep history.

The receiver face must follow the same center and have the selected
**inertial
coordinate radius**. With uniform inner-region expansion, rotation,
translation, and identity skew, set all angular shape coefficients to zero
and
prescribe the size monopole as

\f[
\lambda_{00}(t)=\sqrt{4\pi}\left(R_{\rm grid}-R_*/a(t)\right).
\f]

The existing shape map supplies the spatial transition to the rest of the
mesh. The new `domain::FunctionsOfTime::FixedInertialRadius` owns a clone of
the archived expansion history and evaluates the monopole and its
derivatives analytically, including across expansion updates. It does not
accept control updates. For a recorded BBH, configure B as follows:

```yaml
ShapeMapB:
  LMax: Auto
  CoefficientTruncationLimit: 0.0
  TransitionEndsAtCube: true
  H5Filename: /path/to/BbhVolume0.h5
  SubfileName: ReplayVolume
  Replay: true
  FixedInertialRadius: 0.30
```

The grid inner radius is still supplied by `ObjectB.InnerRadius`. Choosing
`R_grid = R_*/a(initial_time)` starts with zero size deformation. The
fixed-radius option requires an archived scalar `Expansion` and the donor's
shape/size records. It assumes inner-region uniform expansion and identity
skew; the boundary rejects faces that do not lie on the selected sphere. It
also checks that the mapped receiver center agrees with the recorded centers
at the six sample times used for temporal interpolation. The exact mapped
center at the current time defines the angles and radius check. The center
itself is not temporally interpolated: coarse sampling of an accelerating
trajectory would otherwise spuriously fail the geometry check. The archived
map histories must cover the full interpolation stencil, including its
future sample times. The reader retains one evaluation, keyed by the exact
time, inertial point coordinates and mapped center. Repeated ghost and
model-correction queries reuse the interpolated fields; each query still
validates the recorded center trajectory. The cache is synchronized between
worker threads and is not part of serialized state. The physical model and
its sector projections also share one set of Bjørhus intermediate quantities
within each RHS evaluation.

For all remaining map imports, add `Replay: true`. This loads the **global**
functions-of-time history without reinitializing it. Omitting `Replay`
retains the ordinary numerical-restart behavior, using observation-local
functions and reinitializing them. Omit `FixedInertialRadius` on A to
preserve its archived shape and size. Disable all control systems
(`Expansion: None`, etc.), and remove horizon-finding events for the excised
object. Keep `UseWorldtube: false`: that domain option selects a different
worldtube evolution system.

Initialize the receiver from the donor's volume fields with
`ElementsAreIdentical: false`, since its maps and inner radial mesh have
changed. Retain the donor's gauge, damping, and numerical flux settings. The
ghost boundary passes the reference fields into native DG packaging and
lifting, which include the receiver mesh velocity and handle all
characteristic families. It does not replace any model sector. Spatial
derivatives in the tape remain available for subsequent model diagnostics
but are not required by this all-reference ghost boundary.

Stop before the last recorded surface time, allowing for any future times
requested by startup or the time stepper. A volume snapshot at the final
slab does not guarantee that the dense surface trigger fired at precisely
that endpoint. Establish the all-reference interpolation/discretization
error before introducing a physical-sector model; changing the excision mesh
prevents an expectation of discrete roundoff agreement with the donor.

## Receiver observations on the excision boundary

For replay diagnostics, use the existing `ObservationExcisionBoundaryB`
event
(and the corresponding A event for A). It locates the configured sphere in
**grid coordinates**, so roundoff in the inverse moving maps does not push
boundary points into the excised region. Set its center and radius to the
object's grid excision center and inner radius, and use
`AngularOrdering: Strahlkorper`. The target maps these same points forward
to
save their actual inertial coordinates at each observation time.

```yaml
EventsAndDenseTriggers:
  - Trigger:
      Times:
        EvenlySpaced:
          Interval: 0.01
          Offset: 1617.0
    Events:
      - ObservationExcisionBoundaryB
```

Remove the same event from slab triggers if using this dense trigger, to
avoid duplicate observations at shared times. `LMax: 16` gives 561 angular
points. This is interpolation from the exterior DG elements onto the
excision surface; angular points need not coincide with the DG face nodes.
At tangential element interfaces, the interpolation framework chooses the
first received trace.

The `/ObservationExcisionBoundaryB.vol` subfile in the surface H5 file
contains all 204 inertial GH/derivative/curvature components described
above, plus `Lapse`, the grid-frame `Shift`, and
`InertialCoordinates_x/y/z`. Spatial derivatives are computed on the volume
mesh before interpolation. The existing grid-frame surface harmonics are
still written to the reductions file. The same extensions apply to
`ObservationExcisionBoundaryA`.

For comparison with the reference tape, evaluate the donor data at the saved
**inertial positions and times**. The receiver's angular grid follows the
grid
frame and generally rotates relative to the donor's fixed inertial angular
axes, so matching array indices is incorrect. This diagnostic subfile is not
input for `WorldtubeReplay`: it uses grid angular coordinates and does not
contain the donor tape's inertial-center/radius metadata.

First check that every scheduled observation was written and all fields are
finite. Compare the receiver's incoming characteristic sectors and local
constraints with the donor at those positions. Repeat with higher polynomial
resolution near the excised object, keeping the same tape, observation
cadence, and angular resolution. A reduction in error establishes receiver
convergence; a plateau calls for separate time-step, tape-cadence, and
angular-resolution checks before introducing a model.

## Replacing only the physical sector by the order-zero type-D model

After validating all-reference replay, replace the boundary-condition name
`WorldtubeReplay` by `WorldtubeReplayTypeD`, keeping the same three options:

```yaml
ExciseWithBoundaryCondition:
  WorldtubeReplayTypeD:
    H5Filename: /path/to/BbhSurfaces.h5
    SubfileName: GhWorldtube
    Radius: 0.30
```

This uses the existing `EvolveGhBinaryBlackHole` executable and the same
reference sphere, initial volume data, prescribed maps and observations. It
requires `AveragedUpwindPenalty` and Gauss-Lobatto points in the direction
normal to the boundary. Moving maps and radial/angular element refinement
are supported; the lift uses the receiver element's inverse Jacobian and
normal polynomial extent.

The ghost data remain the recorded GH fields. An additional time-derivative
correction replaces only the physical projection of the incoming `u_minus`
reference-interface correction with the existing `WorldtubeTypeD` Bjorhus
correction. Schematically, if `C_ref` is the lifted replay correction and
`C_D` the model correction, the extra characteristic correction is
`P_phys(C_D_minus - C_ref_minus)`. It is converted back to evolved RHSs with
`delta g = 0`, `delta Pi = delta u_minus / 2`, and `delta Phi_i = -n_i delta
u_minus / 2`. Thus the reference gauge and constraint projections, and the
other characteristic families, are preserved. Points where `u_minus` is
outgoing receive no extra correction.

The model reads the **live receiver curvature**, computed from its volume
spatial derivatives at the face. It supplies the order-zero type-D incoming
Weyl mode (`Psi0 = 6 b^2 Psi2^K`) through the existing matching
implementation. It needs neither a prescribed mass nor the order-two
Kretschmann history. The recorded derivatives remain diagnostics; they are
not substituted for the receiver derivatives in the model.

Compare against an all-reference run at the same receiver resolution. Assess
the initial transient separately, then compare the incoming physical sector,
local constraints and exterior fields. This is a sector-isolation
experiment, with the constraint sectors supplied by replay by default.

To use live constraint preservation together with TypeD physics and replay
gauge, add `ConstraintPreserving: true` under `WorldtubeReplayTypeD`. This
uses the existing Bjørhus constraint terms for `u_psi`, `u_zero`, and the
constraint projection of `u_minus`. Each replacement acts only where that
field's moving-mesh characteristic speed is negative. The additional
correction subtracts the corresponding replay lift before adding the model
correction, so the two boundary prescriptions are not added together. The
gauge projection of `u_minus` and all of `u_plus` retain replay. Omitting
this option is equivalent to `ConstraintPreserving: false`. No relaxation is
introduced.

To replace the recorded gauge by the live inner-boundary Sommerfeld
condition, add `SommerfeldGauge: true` under `WorldtubeReplayTypeD`. This
selects the existing `WorldtubeTypeD` `SommerfeldAbsorbing` gauge terms,
using the current inertial distance from the mapped excision center,
including on moving meshes. It uses the inner-boundary sign, annihilating
`g(t+r)/r` with `(partial_t - partial_r - 1/r) f = 0`. The additional
correction replaces only the incoming gauge projection of `u_minus`,
subtracting its replay contribution first. It is independent of
`ConstraintPreserving` and introduces no relaxation. Omitting both
model-gauge options retains recorded gauge. `u_plus` and all outgoing or
zero-speed fields keep their reference contributions.

The experimental `QhatRadialResponse` option instead models the radial
derivative of `qhat_a = ell^b (Pi_ab - n^i Phi_iab)`, with `n` pointing into
the excised region. Specify three dimensionless coefficients `R K` for the
time, radial, and tangential channels, for example:

```yaml
QhatRadialResponse: [-2.654, -2.332, -2.106]
```

These example values come from one q=8 donor interval and are not
established as universal or stable coefficients. `None` disables the option,
which is mutually exclusive with `SommerfeldGauge: true`.

The model captures the receiver's initial `qhat` and its full radial
derivative and imposes `dr(qhat-qhat_initial) = K (qhat-qhat_initial)` in
all four gauge components. Spatial components of the reference rotate with
the `Rotation` function of time; angular labels must follow the same
rotation. The derivative includes the metric dependence of both null and
spatial normals. Only the incoming gauge correction is replaced; live GH
source terms remain. Changes to the characteristic basis caused by other
sectors' metric corrections are included. The reference is serialized and is
not recaptured during time-step self-start; a field-only restart captures a
new reference. No future donor gauge trace enters this model correction.

This prototype requires fixed-radius coordinate spheres on the lower radial
Gauss-Lobatto face of spherical-harmonic shells, with fixed mesh resolution
and no AMR. Initial radial refinement is supported, as are moving spheres
and time-dependent constraint damping. By default, `u_plus` and any sectors
not replaced by model options still use the replay interface.

To remove all dependence on recorded boundary fields, set
`ReplayBoundaryData: false`. The historical class name remains unchanged:

```yaml
WorldtubeReplayTypeD:
  ReplayBoundaryData: false
  Radius: 0.3
  PhysicalModel: Quadrupole
  Mass: 0.1111111111111111
  MomentRelaxationTime: None
  ConstraintPreserving: true
  QhatRadialResponse: [-2.654, -2.332, -2.106]
```

In this mode `H5Filename` and `SubfileName` may be omitted and are never
read. The ghost fields equal the receiver fields, so `AveragedUpwindPenalty`
contributes exactly zero. The time-derivative correction directly applies
the incoming physical, constraint, and gauge prescriptions. Outgoing and
zero-speed characteristics receive no boundary correction. This mode
requires constraint preservation and either gauge model; it returns an error
if any `u_plus` speed becomes negative because no incoming `u_plus` model is
supplied. Boundary replay remains enabled by default for existing inputs.
Coordinate maps may still be prescribed from recorded functions of time,
independently of this option, and initial fields may be imported.

When using spherical-harmonic shells around the objects with
`TransitionEndsAtCube: false`, restrict each distorted-frame
apparent-horizon search to its `ObjectAShell` or `ObjectBShell` block group.
The neighboring `ObjectACube` and `ObjectBCube` blocks have no distorted
frame in this setup. This restriction also applies to the control-system
horizon searches that share those options through YAML anchors. An
inertial-frame common-horizon search can still include the cube blocks.

The same physical-sector replay boundary can select the quadrupole model by
adding these options under `WorldtubeReplayTypeD` (the historical name is
retained):

```yaml
PhysicalModel: Quadrupole
Mass: 0.1111111111111111
MomentRelaxationTime: None
```

This example uses a hole of mass `1/9` in total-mass units. `None` imposes
the instantaneous tidal fit with no temporal relaxation. The live
Kretschmann scalar and its gradient are computed in the volume before
slicing onto the face. Its inertial time derivative uses the existing
backward difference and mesh-velocity correction. The BBH action updates
these data before the boundary condition is evaluated. The five tidal
moments are fitted on each element's excision face; on refined
tensor-product wedges this is a local patch fit, not a collective fit over
the whole sphere. To use tensor-product meshes, set `UseSphericalHarmonics:
false` explicitly under both `ObjectA` and `ObjectB`. The other sectors
still come from replay, exactly as in the TypeD experiment.

To isolate feedback through the physical model, add `PhysicalTarget:
RecordedModel` to `WorldtubeReplayTypeD`. The default, `Live`, retains the
behavior above. `RecordedModel` constructs the entire physical target from
the recorded metric, Pi, Phi, volume-computed derivatives of Pi and Phi, and
Kretschmann scalar and its volume-computed spatial gradient. This includes
the tetrad, rest-frame boost, Coulomb field and instantaneous tidal fit;
receiver fields and stored receiver moments are not inputs to the target.
Order-two recorded targets require a full spherical-harmonic face.

The recording has no stored time derivative of Kretschmann. The reader
differentiates its six-sample temporal interpolant at fixed recording angle
and subtracts the exact mapped center velocity contracted with the recorded
spatial gradient. The recording sphere has fixed inertial radius and angular
axes and translates with its center. No angular derivative is used to
replace a recorded spatial derivative, and no radial interpolation is
performed.

`PhysicalTarget: RecordedWeyl` provides a control: it prescribes the
recorded true incoming Weyl tensor instead of a model fitted to it. Both
recorded modes require `ReplayBoundaryData: true` and `MomentRelaxationTime:
None`. They still apply the physical correction with the receiver's live
Bjørhus operator, projections and characteristic speeds.
`ConstraintPreserving` and the gauge options keep their usual meanings;
`ConstraintPreserving: true`, `SommerfeldGauge: false` and
`QhatRadialResponse: None` select live CP and recorded gauge. These tests
isolate feedback through the physical target; they do not replace the
complete boundary RHS by donor values.

For the intermediate feedback test, select `PhysicalTarget: RecordedMoments`
with `PhysicalModel: Quadrupole`. This fits only the five complex tidal
moments from the donor (including the donor frame and boost needed for that
fit), then constructs the target using those moments with the receiver's
live frame, boost, curvature and retained Weyl scalars. Live stored moments
are ignored; there is no temporal relaxation. The receiver Kretschmann
history is still updated, so its boost uses the usual live backward time
difference and volume-computed spatial gradient. This mode requires boundary
replay and a full spherical-harmonic face, and also supports
`QuadrupoleCoulomb` donor fits. CP and gauge selection are unchanged.


See \ref tutorial_geometric_worldtube for the
geometric, third-order and longitudinal physical models.

## Derivative-input and electric/magnetic feedback experiments

For `PhysicalTarget: Live` with `PhysicalModel: QuadrupoleGeometric`, the
optional fit-only derivative substitution is:

```yaml
DerivativeFitReplay:
  LiveStrength: 0.0
  RampStart: 15.0
  RampDuration: 3.0
```

The numbers are example simulation times. This requires boundary replay and
`MomentRelaxationTime: None`. The strength begins at one and changes
smoothly to `LiveStrength`; one is the live control and zero selects
recorded `partial_i Pi_ab` and `partial_i Phi_jab` in the moment fit. These
are the inertial spatial derivatives supplied to the electric/magnetic Weyl
construction. Metric, Pi, Phi, the live frame and geometric map, and the
live Kretschmann data remain receiver quantities. The resulting moments
enter the usual live target assembly. CP, gauge selection and the physical
Bjorhus application retain their selected prescriptions. This is a
diagnostic mixed input, not a claim that the recorded derivatives are more
accurate.

A different diagnostic uses `PhysicalTarget: RecordedMoments` with the
legacy `PhysicalModel: Quadrupole`:

```yaml
MomentFeedback:
  Strength: [1.0, 0.0]
  StartTime: 15.0
  RampDuration: 3.0
  MomentsFile: control-moments.csv
  GeometryFile: control-geometry.bin
```

`Strength: [Electric, Magnetic]` independently controls the two
fitted-moment feedback gains. A scalar applies the same gain to both. The
imposed moments are

```text
Q_recorded + alpha_E Re(Q_live - R Q_control R^T)
           + i alpha_B Im(Q_live - R Q_control R^T) + R Q_pulse R^T
```

The reference files must come from the independently evolved
prescribed-moment control. `R` aligns its Cholesky axes with the receiver
using the proper polar rotations and face quadrature. The geometry recording
checks coordinate-map compatibility; extrapolation beyond reference times is
rejected. `[1, 0]` feeds back electric departures from the control; `[0, 1]`
feeds back magnetic departures. These are reference-centered perturbation
experiments, not removal of the other sector's recorded background. The
gains share a sin-squared ramp. `MomentPulse` can provide a controlled
forcing, and `MomentDiagnosticsInterval` requests full-step records of
donor, imposed and independently fitted moments.

## Experimental physical-trace boundaries

`WorldtubePhysicalGhost` prescribes the two incoming physical trace
components through the numerical ghost flux, with live CP and the selected
gauge correction. `WorldtubeEvolvedPhysicalGhost` instead evolves those two
amplitudes with the bulk time stepper. Its forcing may be prescribed
curvature, a quadrupole on a prescribed Schwarzschild frame, or recorded
incoming Weyl curvature. It does not perform a live NP moment fit.

These are experimental alternatives to the physical Bjorhus prescription.
They require gamma1 = -1; the evolved trace additionally requires global
time stepping and no surface AMR. Prescribed forcing uses a fixed sphere;
the recorded-Weyl mode supports moving maps with recorded gauge and live CP.
The `KerrSchildPhysicalGhost.yaml`, `KerrSchildEvolvedPhysicalGhost.yaml`
and `KerrSchildPrescribedNpGhost.yaml` inputs in
`tests/InputFiles/GeneralizedHarmonic/` exercise the supported
configurations.
