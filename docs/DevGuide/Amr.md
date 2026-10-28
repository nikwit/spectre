\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond

# Adaptive mesh refinement (AMR) {#dev_guide_amr}

\tableofcontents

### Introduction

SpECTRE implements a block-based anisotropic h-p adaptive mesh
refinement (AMR) algorithm.  Currently AMR can only occur at global
synchronization points, but it is planned to remove this restriction
in the future.

\warning AMR is still under active development. The elliptic
executables fully support AMR.  The evolution executables currently
only support p-adaptivity, and only when using a global time stepper,
and only when using a DG method.  The executable will abort with an
ERROR if local time stepping or h-refinement is done.

Here is an overview of what is described in detail in the sections below:

- \ref dev_guide_amr_algorithm "AMR Algorithm": An overview of the
  parallel phases, actions, and components that are used to do
  adaptive mesh refinement.
- \ref dev_guide_amr_criteria "AMR Criteria": Describes the
  expectations and restrictions on concrete classes implementing
  adaptive mesh refinement criteria.
- \ref dev_guide_amr_policies "AMR Policies": Describes restrictions
  and assumptions of the current algorithm, and what would be
  required to relax them.
- \ref dev_guide_amr_projectors "AMR Projectors": Describes the
  requirements on concrete classes which project data to the new or
  modified elements after mesh refinement.
- \ref dev_guide_amr_metavariables "AMR Metavariables": Describes the
  code that needs to be added to the metavariables of an executable in
  order to enable adaptive mesh refinement.
- \ref dev_guide_amr_input_options "Controlling AMR": Describes the
  options in the input file that control the behavior of adaptive mesh
  refinement.

### AMR Algorithm {#dev_guide_amr_algorithm}

The AMR algorithm requires communication via simple actions between
the elements of the array component that are responsible for the
Element%s of the Domain and a singleton amr::Component.  See \ref
dev_guide_parallelization_foundations
"Parallelization and Core Concepts" for an overview of phases,
actions, metavariables, and other SpECTRE parallelization concepts.

Currently the AMR algorithm is run over several phases when a phase
change is triggered by an executable.  See \ref
dev_guide_amr_input_options "Controlling AMR" for how to control when
the AMR phases are triggered.

#### Evaluating refinement criteria

When AMR is triggered, the executable should enter
Parallel::Phase::EvaluateAmrCriteria.  In this phase, the
amr::Component triggers the simple action
amr::Actions::EvaluateRefinementCriteria on each element of the array
component (specified by `metavariables::amr::element_array`).  Each
element evaluates (in order) the user-specified refinement criteria,
each of which returns a refinement decision in each logical dimension
of the element represented as an amr::Flag.  The overall refinement
decision in each logical dimension is obtained by taking the highest
priority decision in the dimension over all criteria, where
amr::Flag::Split has the highest priority and amr::Flag::Join has the
lowest.  See \ref dev_guide_amr_criteria "AMR Criteria" for more
details on implementing criteria.  Once the overall refinement
decision based on the refinement criteria is made, the decision may be
modified to satisfy the \ref dev_guide_amr_policies "AMR Policies".
Finally, the element sends its decision (and other information in
amr::Info) to its neighbors by calling the simple action
amr::Actions::UpdateAmrDecision.

In amr::Actions::UpdateAmrDecision, an element decides whether or not it
needs to update its own refinement decision (or other information in
amr::Info) based on decisions (or other information in amr::Info) sent
by its neighbor.  Possible reasons to update the refinement decision
are, e.g., if a sibling neighbor does not want to join an element that
had initially decided to join, or if any of the amr::Policies would
be violated.  If the refinement decision (or other information in
amr::Info) of the element changes, the element will send its
updated amr::Info to its neighbors.

Note that simple actions are executed asynchronusly, so both
amr::Actions::EvaluateRefinementCriteria and
amr::Actions::UpdateAmrDecision have to handle the possibility that
they are executed out of order.  For example,
amr::Actions::UpdateAmrDecision could be executed before
amr::Actions::EvaluateRefinementCriteria.  Eventually all elements
will reach a state where they do not need to send new information to
their neighbors, thus reaching a state of quiescence.  When this
happens Parallel::Phase::EvaluateAmrCriteria ends. At this point all
elements have reached a consensus on how the domain should be adapted,
but no changes have been made to the domain.  In order to adjust the
domain, the executable should next start
Parallel::Phase::AdjustDomain.

#### Adjusting the domain

At the beginning of Parallel::Phase::AdjustDomain, the amr::Component
calls amr::Actions::AdjustDomain on each element of the array
component (specified by `metavariables::amr::element_array`).  If the
element wants to split in any dimension, it will call the simple
action amr::Actions::CreateChild on the amr::Component in order to
create new refined elements that will replace the existing element.
If the element wants to join in any dimension, it will either call
amr::Actions::CreateParent on the amr::Component (if it is the joining
sibling with the lowest SegmentId in each dimension being joined) or
do nothing in order to create a single new coarsened element that will
replace multiple existing elements.  If the element neither wants to
split nor join in any dimension, then the element mutates items in its
db::DataBox.  First the Element, Mesh and neighbor Mesh are updated.
Then the element calls the
\ref dev_guide_amr_projectors "AMR Projectors" that will update all
the mutable items in the db::DataBox.
Finally the amr::Info of the element and its neighbors are reset and
cleared so they are ready for a future invocation of AMR.

When amr::Actions::CreateChild is executed by the amr::Component, it
creates a new element in the `metavariables::amr::element_array`
passing a Parallel::SimpleActionCallback that will be executed after
the new element is created.  The callback will be either to call
amr::Actions::CreateChild again for the next child to be created, or
amr::Actions::SendDataToChildren which will be invoked on the element
that is being split.  amr::Actions::SendDataToChildren will then
invoke amr::Actions::InitializeChild on each new child element,
passing along all of the mutable items in its
db::DataBox. amr::Actions::SendDataToChildren will then deregister and
delete the element being split. amr::Actions::InitializeChild will
create its Element, Mesh, and neighbor Mesh, and then call the \ref
dev_guide_amr_projectors "AMR Projectors" that will update all the
mutable items in its db::DataBox.

When amr::Actions::CreateParent is executed by the amr::Component, it
creates a new element in the `metavariables::amr::element_array`
passing a Parallel::SimpleActionCallback that will be executed after
the new element is created.  The callback will be
amr::Actions::CollectDataFromChildren which will be invoked on one
of the elements that are joining.
amr::Actions::CollectDataFromChildren will then invoke either itself
on another element that is joining or amr::Actions::InitializeParent
on the new parent element, collecting all of the mutable items in each
joining element's db::DataBox. amr::Actions::CollectDataFromChildren
will then deregister and delete the element that is
joining. amr::Actions::InitializeParent will create its Element, Mesh,
and neighbor Mesh, and then call the \ref dev_guide_amr_projectors
"AMR Projectors" that will update all the mutable items in its
db::DataBox.

When all of the actions finish, a state of quiescence is reached.
When this happens Parallel::Phase::AdjustDomain ends.  The executable
should now either return to the phase which triggered AMR, or enter
Parallel::Phase::CheckDomain.

### AMR Criteria {#dev_guide_amr_criteria}

When amr::Actions::EvaluateRefinementCriteria is executed, for each
element, it loops over a list of criteria that is specified in the
input file options
(See \ref dev_guide_amr_input_options "Controling AMR").
Each criterion must be an option-creatable class derived from
amr::Criterion (see its documentation for a list of requirements for
the derived classes).  In particular, each criterion should compute an
array of amr::Flag%s that represent the recommended refinemenent
choice in each logical dimension.  These choices can be computed from
any items in the db::DataBox of the element or additional compute
items specified by the criterion.  It is expected that a criterion
will compute something and then either choose to recommend refinement,
no change, or coarsening the element either in size (h-refinement) or
polynomial order (p-refinement).  A criterion does not need to handle
anything that is enforced by the amr::Policies such as worry about
bounds on the refinement level or number of grid points.

### Independent spectral extents in the GH SpEC-style controller

`gh::spec_amr::RefinementVariables` describes the independent refinement
variables of an element from its local mesh. A polynomial element has three
independent orders. A spherical-harmonic shell has a radial order and a single
angular order, with `m_max == l_max`. Its two angular mesh dimensions share one
spectral diagnostic and one history of constraint forecasts, target shifts and
cooldowns. The controller expands the angular decision into logical AMR flags
only on return. In its fixed-size diagnostic and history arrays, only the
first `RefinementVariables::size()` entries are active; unused entries do not
participate in decisions.

`gh::SpecSpectral` uses `gh::spec_amr::diagnose_spec` on both polynomial
(Legendre or Chebyshev) and `I1 x S2` elements. Polynomial and shell-radial
power monitors transform only the monitored direction, leaving transverse
points nodal. Shell-angular power uses tensor spherical harmonics in the grid
frame, with spin-dependent coefficient counts. The Pi/Phi group combines
squared coefficients and counts before normalization. It follows SpEC's
storage convention, including the zero imaginary m=0 scalar slot.

The SpEC diagnostic compacts exact zeros for weighted tails while retaining
original indices for slope fits. Projected constraints use SpEC's separate
element-to-grid and grid-to-inertial normalization, then its topological
aggregation: `sqrt(RMS_d)` for I1 and
`sqrt((RMS_theta + RMS_phi)/2)` for S2. These are deliberate changes to the
previous `SpecSpectral` polynomial diagnostic conventions, so previous
constraint tolerance values should not be interpreted as unchanged norm
thresholds. The older `ChalmersWithSpecCoarsening` controller and the legacy
`diagnose` / `diagnose_power` helpers retain their original diagnostics.

The element's actual `Filters::Tags::SpectralFilter` supplies its modal volume
response. The effective filtered count is `sum(1-response)`, rounded by the
spectrum diagnostic. The filter is also queried at each prospective lower
independent extent for the removal forecast. Disabled volume filters and
`Filters::None` give zero filtered modes; boundary filtering does not change
this count. The count describes a single filter application, not accumulated
damping over many timesteps. Missing filter metadata is an error instead of
an assumption of one filtered mode. A prospective mesh unsupported by the
local filter inhibits removal. `Hypercube`, `SphericalShell` and `None`
provide this query; other filters must implement it before use here.

Cartesian-to-tensor-harmonic matrices are cached per element, replaced when
l_max changes, and rebuilt after migration. They are not numerical history.
The minimum-point bounds remain in logical mesh coordinates; the two angular
bounds are combined into a single S2 floor. This is a diagnostic/controller
adapter, not an implementation of adaptive h refinement, adaptive filters or
the rest of SpEC's domain controller. Reference fixtures and their scope are
documented in
`tests/Unit/Evolution/Systems/GeneralizedHarmonic/SpecAmrReference/README.md`.

At a nonconforming polynomial / spherical-harmonic interface, p AMR preserves
the connectivity and equal-rate stepping policy. Each side uses its own face
mesh as its mortar; the wedge-side interpolator is rebuilt using both updated
meshes. The S2 side keeps a single self-keyed mortar rather than storing every
wedge in its fixed-size neighbor-mesh table. Logical h levels on opposite
sides are not comparable, so the two-to-one check is skipped there. Splitting
or joining an element adjacent to such an interface is not supported.
Resolution changes must run through `EvaluateAmrCriteria` / `AdjustDomain`
so the two sides exchange their proposed meshes before rebuilding the
interpolators. The legacy `RefineMesh` event does not provide that exchange.
Conservative LTS on other faces retains its history and projects only the
local geometric data. Volume time-stepper values and derivatives are
projected without restarting the integrator. The mixed-interface projector
regression checks both sides, changes on neighbors alone, and serialization
between AMR calls; the scalar-wave input
`tests/InputFiles/ScalarWave/NonconformingSphericalWaveAmrLts.yaml` exercises
repeated p increases and decreases during evolution with LTS.

### AMR Policies {#dev_guide_amr_policies}

The current AMR algorithm has a set of rules that are enforced after
the overall refinement choice is determined from the AMR criteria.

The following rules simplify the code:

- An Element cannot join if it is splitting in any dimension.  Instead
  any amr::Flag::Join is changed to amr::Flag::DoNothing. Relaxing
  this restriction would be a non-trivial change that would require
  more complicated logic in determining if neighboring elements can
  join, as well as an additional code path to handle the creation and
  initialization of the new elements. This should not be a significant
  restriction as in most cases when this would occur it is unlikely
  that the siblings of the Element would agree to a valid join.

- In two or three spatial dimensions, there must be a 2:1 balance
  (i.e. within one refinement level) between neighbors in the
  dimension(s) parallel to the face of the element.  Relaxing this
  restriction would be a significant change to the code as the
  code for communicating boundary corrections assumes 2:1 balance.  In
  order to maintain 2:1 balance with neighboring elements, a choice of
  amr::Flag::Join may be changed to amr::Flag::DoNothing, and
  amr::Flag::DoNothing may be changed to amr::Flag::Split.

The following policies are among those controlled by input file
options (see \ref dev_guide_amr_input_options "Controlling AMR" and
the documentation for amr::Policies):

- Whether to do isotropic or anisotropic refinement.  If isotropic
  refinement is done, the AMR decision in each dimension are changed
  to the decision of the dimension with the highest priority amr::Flag
  (i.e. amr::Flag::Split has the highest priority).

- Whether or not to enforce 2:1 balance for the refinement level in
  the direction perpendicular to the face.

- The minimum and maximum of the refinement level and the number of
  grid points.  The actual bounds on these will be the stricter of
  those specified by the input file and those set by the code.  For
  the number of grid points, the code bounds are set by
  Spectral::minimum_number_of_points and
  Spectral::maximum_number_of_points for the Spectral::Basis and
  Spectral::Quadrature of the Mesh.  For the refinement level, the
  minimum is zero, and the maximum is set by
  `ElementId::max_refinement_level`.

### AMR Projectors {#dev_guide_amr_projectors}

After the grid adapts, the mutable items in the DataBox of new or
existing elements must be updated.  The Element, Mesh, neighbor Mesh,
amr::Info, and neighbor amr::Info along with the items in
`Parallel::Tags::distributed_object_tags` are automatically updated by
the AMR actions, but all other items must be explicitly updated via
projectors that conform to amr::protocols::Projector.  If any mutable
item has not been explicitly handled by any of the projectors, a
`static_assert` will be triggered listing the tags for the items that
have not been projected.

When a new element is created, the items in its DataBox are default
constructed.  Their Element, Mesh, and neighbor Mesh are mutated to
their desired state.  Then the AMR projectors are called, passing
along the items in the DataBox(es) from their splitting parent element
or joining children elements.  Existing elements first mutate their
Element, Mesh, and neighbor Mesh, and then call the AMR projectors,
passing along the old Element and Mesh.

The `return_tags` of each projector lists the tags for the items that
the projector mutates.  From the pre-AMR state, the projector must
compute the post-AMR state of these items.  Typically during
Parallel::Phase::Initialization, a group of items will be initialized
by a specific action or mutator in the action
`Initialization::Actions::InitializeItems`.  Therefore it makes sense
to create a projector that will handle the same group of items.
However some items are initialized from input file options, while
others are left default initialized and not mutated by any
initialization action.  These items will still need to be handled by
some projector.  Two convenience projectors are provided:
amr::projectors::DefaultInitialize which value initializes the items
corresponding to the listed tags; and
amr::projectors::CopyFromCreatorOrLeaveAsIs that will either copy the
item from the parent (or children, asserting that all children agree)
of a newly created element, or leave the item unmodified for an
existing element for the items corresponding to the listed tags.

The list of projectors is specified in
`metavariables::amr::projectors` in the executable, and the projectors
are evaluated in the order in which they are specified.

### Enabling AMR for an executable {#dev_guide_amr_metavariables}

In order to enable AMR for an executable, the following changes need
to be made in the metavariables:

- The addition of a struct named `amr` such as the following example:
```
struct amr : tt::ConformsTo<::amr::protocols::AmrMetavariables> {
    using element_array = dg_element_array;
    using projectors = tmpl::list<::amr::projectors::DefaultInitialize<
        domain::Tags::InitialExtents<Dim>,
        domain::Tags::InitialRefinementLevels<Dim>,
        evolution::dg::Tags::Quadrature>>;
    static constexpr bool keep_coarse_grids = false;
  };
```

  where `element_array` specifies the array component on which AMR
  will operate and `projectors` is a type list of amr::projectors
  that govern how the items in the DataBox for the `element_array` are
  updated after refinement.

- The addition of amr::Component to the `component_list` of the metavariables

- The addition (or modification) of the following `tmpl::pair`s in
  `factory_creation::factory_classes`:
```
        tmpl::pair<amr::Criterion,
                   tmpl::list<LIST_OF_CRITERIA>,
        tmpl::pair<
            PhaseChange,
            tmpl::list<
                PhaseControl::VisitAndReturn<
                    Parallel::Phase::EvaluateAmrCriteria>,
                PhaseControl::VisitAndReturn<Parallel::Phase::AdjustDomain>,
                PhaseControl::VisitAndReturn<Parallel::Phase::CheckDomain>>,
```
where `LIST_OF_CRITERIA` should be a list of amr::Criteria.

- The addition of the following item in the list of PhaseActions for
  the component specified in `amr::element_array`:

```
          Parallel::PhaseActions<Parallel::Phase::CheckDomain,
                                 tmpl::list<::amr::Actions::SendAmrDiagnostics,
                                            Parallel::Actions::TerminatePhase>>,
```

- The possible addition of `PhaseControl::Actions::ExecutePhaseChange`
  to the action list for the appropriate Phase in the PhaseAction of
  the `amr::element_array`

- The addition of the mutator `amr::Initialization::Initialize` to the
  list of mutators in `Initialization::Actions::InitializeItems` in
  the Parallel::PhaseActions list for Parallel::Phase::Initialization.

- The appropriate includes for all of the above.

### Controlling AMR {#dev_guide_amr_input_options}

There are two places in the input file that control how and when AMR happens.

The "how" is controlled by the option group `Amr`.  Here you list the
amr::Criteria being used, the available options for the amr::Policies, and
the `Verbosity` of diagnostic messages that are printed.  Here is an
example:
```
Amr:
  Criteria:
    - TruncationError:
        VariablesToMonitor: [Psi]
        AbsoluteTarget: 1.e-6
        RelativeTarget: 1.0
  Policies:
    EnforceTwoToOneBalanceInNormalDirection: true
    Isotropy: Anisotropic
    Limits:
      RefinementLevel: Auto
      NumGridPoints: Auto
      ErrorBeyondLimits: False
    AllowCoarsening: True
  Verbosity: Verbose
```

Note that the values `Auto` for the Limits options choose the default
limits set by the code. Also note that if a criterion tries to refine beyond the
limits, whether or not the code should error is controlled by
`ErrorBeyondLimits`.

"When" AMR happens is controlled by specifying a Trigger and a list of
`PhaseChanges` in the top-level option `PhaseChangesAndTriggers`.  For
example:
```
PhaseChangeAndTriggers:
  - Trigger:
      Slabs:
        EvenlySpaced:
          Interval: 10
          Offset: 0
    PhaseChanges:
      - VisitAndReturn(EvaluateAmrCriteria)
      - VisitAndReturn(AdjustDomain)
      - VisitAndReturn(CheckDomain)
```

Both `EvaluateAmrCriteria` and `AdjustDomain` are required in order
for AMR to work.  `VisitAndReturn(CheckDomain)` performs diagnostics
and can be omitted.  To turn off AMR, omit the three phase changes above.

#### Initial p adjustment and junk-radiation timing for BBH

`EvolveGhBinaryBlackHole` can adjust its initial polynomial orders before
time-stepper self-start. With `SpecSpectral` as the sole AMR criterion, add
the following entry alongside the normal evolution AMR phase changes:

```yaml
PhaseChangeAndTriggers:
  - Trigger: {Not: Always}
    PhaseChanges:
      - InitialPRefinement:
          MaxIterations: 40
```

The controller runs automatically during startup; the false trigger prevents
evolution from requesting it again. Each iteration evaluates `SpecSpectral`,
applies the existing p-refinement policies and neighbor updates, and resamples
the original initial data on the resulting mesh. It then reapplies the initial
gauge constraints before evaluating again. Time, the h mesh and element
ownership stay fixed. Startup bypasses the criterion's time gates and cooldown,
permits coarsening on its first call, and does not advance the evolution
decision history. Unlike SpEC's multi-mode initial changes, each iteration
changes at most one point per direction.

Iterations stop when no permitted changes remain or `MaxIterations` is reached;
an element also stops requesting changes if its startup mesh cycles. The final
mesh always receives resampled data, including at the iteration limit. Output
reports changed elements, increasing directions blocked by policies, and
cycling elements. Stopping at a nonfatal point cap or iteration limit does not
mean that all error targets were met. Existing `Amr.Policies.Limits` settings,
including `ErrorBeyondLimits`, still apply. Omitting `InitialPRefinement`
preserves the previous startup behavior.

To suppress adaptation to outgoing junk radiation in the outer wave zone,
replace the scalar `StartTime` inside `SpecSpectral` with, for example:

```yaml
StartTime:
  Time: 0.0
  Blocks: [Envelope, OuterShell0, OuterShell1]
  RadiusPlusTime: 200.0
```

`Blocks` accepts domain block names and block-group names. Select the actual
wave-zone groups in the domain; the names above include the envelope of a
partitioned `BinaryCompactObject` domain. The envelope is included because
outgoing junk crosses it before reaching the outer shells. Each selected element freezes all
three p decisions until
`time >= max(Time, max_initial_grid_radius + RadiusPlusTime)`. The maximum
radius is cached from the element's initial grid points, so later p changes
do not move the threshold. Constraint and spectral diagnostics and their
history continue to update during this radial hold. Other blocks follow the
ordinary scalar `Time` gate. Both gates are bypassed during initial p
adjustment, which can also be used independently of the radial hold.

`RadiusPlusTime` is an absolute time offset. For a run starting at `t0`, use
`t0 + delay` if the desired release is `t0 + radius + delay`. SpEC's standard
bound-orbit SKS/SSphKS/CFMS setup uses `delay = max(200, 3 * initial_separation)`;
its SHK setup uses `2 * initial_separation`. Choose the offset to match the
initial-data family. The original scalar `StartTime: 0.0` remains supported
and imposes no radius-dependent restriction.


#### Redistributing BBH elements after p adjustment

The initial element distribution precedes initial p adjustment, so it can
become imbalanced as the polynomial orders change. `EvolveGhBinaryBlackHole`
provides `RedistributeElements` as a phase change. It uses the current mesh and
grid coordinates with the same `Parallelization.ElementDistribution` cost
model and block Z-curve partitioner used for initial placement. `RoundRobin`
also remains supported. Processors reserved by `ResourceInfo` remain excluded.

For example, redistribute after the first slab and every 100 slabs thereafter:

```yaml
PhaseChangeAndTriggers:
  # Keep the normal AMR entries above this entry, so coincident triggers first
  # adjust the mesh, then redistribute it.
  - Trigger:
      Slabs:
        EvenlySpaced:
          Interval: 100
          Offset: 1
    PhaseChanges:
      - RedistributeElements
```

Phase changes are checked after advancing time. Offset 1 avoids waiting for
slab 100 before the first redistribution. There is no separate startup option
or junk-traversal delay for redistribution; the `SpecSpectral.StartTime` hold
restricts mesh adaptation only.

Configure `RedistributeElements` once. Synchronized phases collect the
current costs, deliver the new destinations, migrate the halted elements,
then rebuild equal-rate time-stepping sections in `UpdateSections`.
Before migration, the old section cookies and multicast trees are discarded
at a quiescent slab boundary. Section creation does not reset the fixed LTS
ratios or mortar policies. Checkpoint restart also recreates these sections
before evolution resumes. The BBH p projectors retain the initial h refinement
levels so section membership can be reconstructed after mesh adjustment.
The element IDs, meshes, evolved fields, time-stepper and mortar histories,
and inboxes are retained. Migration calls the same element registration hooks
as the load-balancing phase, without invoking Charm's measured load balancer.
The output reports how many elements will move and the predicted maximum to
mean processor cost before and after redistribution. These are cost estimates,
not measured speedups.

This implementation supports ordinary element arrays with a fixed h mesh
(uniform h within each block), varying p, and equal-rate sections on mixed
polynomial / spherical-harmonic interfaces. Other active Charm array sections
and multigrid copies remain unsupported. It rejects unsupported layouts before
migration. It does not enable h
adaptation, change domain topology, or restart the evolution.
