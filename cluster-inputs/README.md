\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Second-order worldtube validation sweep (grid frame)

Circular geodesic orbit at r_p = 8M, q = 1, M = 1, validating the
`second-order-grid-frame` branch against the convergence predictions of
arXiv:2304.05329.

The orbit is at r_p = 8M rather than the paper's r_p = 5M because the
particle orbit is numerically evolved on this branch (unlike the 2023 code,
which prescribed it analytically): circular geodesics inside the ISCO at 6M
are radially unstable and the orbit would exponentially depart from
circularity on an e-folding time of 1/(omega sqrt(6M/r - 1)).

Inputs: `Rp8_R<worldtube radius>_n<expansion order>.yaml` for worldtube radii
R = 1.6, 0.8, 0.4 at expansion orders n = 1, 2. Evolution runs to t = 4000M;
the paper found steady state between 3000M and 7000M. All inputs pin the
worldtube radius function to the grid excision radius (required by the
grid-frame scheme — the executable errors out otherwise).

## Validation targets

- Relative error scaling of the settled regular-field coefficients with
  worldtube radius: the error of Psi0 should scale as ~R^(n+1), i.e. ~R^2 at
  n = 1 and ~R^3 at n = 2 (the paper measured alpha = 2.07 and 3.08 at
  r_p = 5M). Psi0 is column 11 of `Reductions.h5` `/PsiTaylorCoefs.dat`
  (columns 2-10 are position/velocity/acceleration).
- Absolute comparison of the settled field along the z-axis against the
  analytic mode sum, Eqs. (53)-(54) of arXiv:2304.05329, which is valid at
  any orbital radius. The `PsiAlongAxis1` interpolation target in the inputs
  observes Psi along the z-axis for this purpose.

## Notes

- Smaller worldtube radii step slower (CFL) and may need more resolution;
  the paper increased resolution until the steady state stopped changing.
- The `WorldtubeSingleton` never terminates cleanly at Completion — the
  deadlock report at the end of the run is expected and harmless.
