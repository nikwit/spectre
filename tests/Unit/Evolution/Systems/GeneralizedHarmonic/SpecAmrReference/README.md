# SpEC AMR reference fixtures

These fixtures establish the behavior of the local SpEC checkout and check
the polynomial and spherical-shell adapters in `gh::SpecSpectral`. They do not
establish equivalence of the full AMR driver or mixed-interface evolution.

## Reproduction and provenance

Run from a configured SpECTRE checkout:

```sh
build-release-debug/bin/python-spectre \
  tests/Unit/Evolution/Systems/GeneralizedHarmonic/SpecAmrReference/GenerateReference.py \
  /path/to/SpEC --compiler clang++
```

The generator compiles unmodified arithmetic fragments from the supplied SpEC
source in a temporary directory. Its compatibility classes provide containers,
topology descriptions and a profiling stub. The numerical linear fit is also
extracted from SpEC, rather than replaced by SpECTRE's regression function.
Original SpEC source is not checked into this repository, and SpEC is not a
dependency of the unit tests.

`ReferenceData.json` records the SpEC revision and SHA-256 hashes of all source
files used. Hashes identify local modifications as well as committed content.
`ReferenceData.hpp` contains the same numerical results for the unit tests.
After regeneration, format the header with the repository's `clang-format`.

The fixtures cover:

- Tensor-component and transverse-point averaging of a directional polynomial
  power monitor, executing SpEC's `I1` summation. A nonconstant transverse
  profile also records the existing full-modal averaging difference.
- Error, prospective error after removing one mode, convergence and pile-up,
  with zero, one and three filtered modes; flat and zero spectra; exact parity
  gaps; a missing leading mode; resolved polynomials; and a roundoff tail.
- The final independent-topology grouping in
  `NormOfProjectedDerivConstraint`, given directional nodal RMS values, for
  `I1 x I1 x I1` and `I1 x S2`. This isolates aggregation, not the complete
  time-dependent coordinate transformation.
- Local add/drop permission expressions from `SetAddDropLimits`, with and
  without forced pile-up handling. These do not execute SpEC's full driver,
  neighbor smoothing, startup or adaptive filtering. The tests combine these
  permissions with the `LargeTrunc` condition for their stated simple cases.
- Effective filtered-mode accounting for a supplied smooth filter response.
  The generator executes SpEC's summation and rounds it as the diagnostic
  does. The local-filter tests additionally compare actual SpECTRE filter
  actions against their reported responses and query prospective lower meshes.
- Tensor-harmonic coefficient counts for the metric, Pi and spatial Phi
  sectors, executing SpEC's original storage iterator and spin-weight routine.
  Counts include its imaginary m=0 scalar slot. Radial points are counted
  separately by the adapter.

## Differences found before the refactor

The tests compare compatible cases directly with SpEC. Where semantics differ,
they retain measured pre-refactor polynomial outputs for the legacy
`diagnose_power` and explicitly test the difference. The new
`diagnose_spec_power` is checked against every case, including the gaps and
roundoff tails. No decision thresholds were retuned.

| Item | SpEC source in these fixtures | Legacy polynomial diagnostic |
| --- | --- | --- |
| Smooth spectra and genuine pile-up | Weighted tails and multi-fit slopes | Matches the reference cases, including multiple filter counts |
| Flat nonzero spectrum | Classified as degenerate | Classified as unresolved, with relative error one |
| Zero spectrum | Degenerate | Degenerate |
| Exact gaps and leading zeros | Removes zero entries before forming weighted tails; filtered-count cutoff also depends on the compressed size | Retains mode positions for weighted tails and fits |
| Roundoff and exact polynomial tails | Omits exact zeros; retains nonzero roundoff | Applies a relative roundoff cutoff and recognizes resolved tails |
| Projected constraint grouping | Takes an additional square root after directional nodal RMS, including for an `I1` extent; for `S2`, returns `sqrt((RMS_theta + RMS_phi)/2)` | Uses the directional nodal RMS directly; no spherical grouping |
| Transverse averaging | An individual basis transform leaves transverse coordinates nodal in the `PowerMonitor` path | Legacy polynomial diagnostics use the full tensor-product modal transform |

The extra square root in the constraint fixture was checked against SpEC's
`SimpleNorms::L2Norm`, which already returns a nodal RMS. `diagnose_spec`
now follows that source expression on **both** polynomial and spherical
meshes. It also adopts the nodal transverse averaging and exact-zero spectrum
conventions above. This changes the meaning of the old `SpecSpectral`
constraint targets; the hybrid controller and legacy diagnostic functions
retain their old conventions. An exactly straight fitted spectrum uses the
finite equal-weight limit if all regression errors vanish, avoiding SpEC's
possible 0/0 without an amplitude cutoff.

`Test_SpecAmrShell.cpp` tests static and sheared grid/inertial projections,
radial versus angular add/drop decisions from manufactured GH fields,
filter responses against real filtering, disabled/boundary-only filters,
resolution-dependent removal forecasts, cache invalidation, and wrapper
serialization. `Test_SphericalShellPowerMonitor.cpp` combines its existing
SpEC Kerr-Schild regression with the independently generated counts to test
the grouped Pi/Phi spectrum. These tests do not evolve through mixed
interfaces, validate time-stepper history projection, or test a BBH merger.

## Independent refinement variables

`gh::spec_amr::RefinementVariables` maps each element's local representation to
its independent extents. Polynomial elements have three variables. A spherical
shell has two: radial mode count and `l_max + 1`. The angular decision expands
to both mesh directions only when returning AMR flags. The second variable has
one constraint forecast, target shift and cooldown; the third storage entry is
inactive and ignored by the controller.

The layout and `diagnose_spec` support Legendre/Chebyshev polynomial
directions and an `I1 x S2` shell. `SpecSpectral` selects diagnostics from the
element's local mesh and selected filter. The legacy `diagnose` rejects shells
explicitly. Minimum theta and phi point counts are combined into one angular
floor; no angular h changes are requested.

History serialization now includes the variable layout. This is a checkpoint
format change for a newly built executable; it does not modify existing frozen
run executables or input files.

The transform cache is transient per-element state: copying a history
snapshot shares the cache, while unpacking clears it. It does not add fields
to the serialized numerical history.
