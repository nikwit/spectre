# q=8 temporal ODE consistency test

This focused test uses all ten downloaded full-resolution volume observations
from \(T=950M_{\rm tot}\) through \(995M_{\rm tot}\), separated by
\(5M_{\rm tot}\).

At every time and worldtube radius:

1. the center is found independently from the Kretschmann dipole;
2. the derivative-assisted first-order \(V1+C3\) hierarchy is fitted;
3. the seven-parameter second-order subsystem consisting of symmetric
   spatial-strain rate and time acceleration is fitted;
4. the first- and second-order solves are alternated four times.

The first-order fit returns

\[
 a=\dot q^0,\qquad \Lambda_{(ij)}=L_{(ij)}/M_B,
\]

whereas the second-order response columns return

\[
 \alpha=\ddot q^0,\qquad S_{(ij)}=\dot L_{(ij)}.
\]

Because the volume observations are labeled by simulation time \(T\), while
the map is \(T=t+q^0(t)\), the tested relations are

\[
\alpha=(1+a)\frac{da}{dT},\qquad
S_{(ij)}=M_B(1+a)\frac{d\Lambda_{(ij)}}{dT}.
\]

The time derivatives are obtained from a seven-point cubic local-polynomial
fit. Centered differences and five- and nine-point local-polynomial
derivatives agree with it at the \(0.3\%\) level.

## Result

The current freely fitted second-order coefficients do **not** satisfy these
ODE relations.

- At the most favorable joint comparison,
  \(R=0.34M_{\rm tot}=3.06M_B\) after three alternations, the fitted
  acceleration is about \(277\) times larger than the derivative of the
  first-order clock rate.
- At the same point, the fitted strain-rate norm is about \(116\) times
  larger than the derivative of the first-order strain.
- Over all tested radii and alternations, the smallest strain-rate amplitude
  ratio is still about \(64\).
- The clock acceleration has nearly the same time dependence as the
  first-order derivative, but with a strongly radius-dependent amplitude.
  The strain-rate tensor is nearly orthogonal to, and often oppositely
  aligned with, the first-order derivative.
- Including the \(T\)-to-\(t\) conversion changes the result by less than one
  percent and cannot account for the discrepancy.

For example, at \(R=0.34M_{\rm tot}\) after three alternations,

\[
\ddot q^0_{\rm fit}\simeq-1.3\times10^{-2},
\qquad
(1+\dot q^0)\frac{d\dot q^0}{dT}
\simeq-4.6\times10^{-5},
\]

and

\[
\|\dot L_{\rm fit}\|_F\simeq8\times10^{-3},
\qquad
\left\|M_B(1+\dot q^0)\frac{d\Lambda}{dT}\right\|_F
\simeq7\times10^{-5}.
\]

A manufactured time-dependent data set constructed from exactly the same
first- and second-order response columns recovers both coefficient blocks and
all four local/integrated ODE comparisons to \(10^{-13}\)--\(10^{-14}\).
This verifies the factor of \(M_B\), the orthonormal symmetric-basis
conversion, the time-coordinate conversion, and the temporal differentiation
code.

The interpretation is therefore not that the \(5M\) cadence is too coarse.
Rather, the seven freely fitted second-order coefficients are currently
acting as effective residual-reduction directions. They are absorbing omitted
second-order coefficients and truncation error, so they cannot yet be
interpreted as the physical time derivatives of the recovered first-order
state.

The natural next test is to impose the ODE-predicted acceleration and strain
rate as known drives, then fit the remaining algebraic second-order
coefficients. This directly tests whether temporal information removes the
large degeneracy without sacrificing metric prediction.

## Direct per-slice metric derivative

The follow-up analysis constructs the metric time derivative independently on
every slice from the evolved generalized-harmonic variables,

\[
\partial_Tg_{ab}=\beta^i\Phi_{iab}-\alpha\Pi_{ab},
\qquad
\partial_Tg^{ab}=-g^{ac}g^{bd}\partial_Tg_{cd}.
\]

The velocity of the curvature-centered virtual worldtube is then included,

\[
D_Tg^{ab}=\partial_Tg^{ab}
v_{\rm center}^i\partial_i g^{ab}.
\]

The center velocity comes from the high-cadence control-system center plus
the derivative of the small curvature-center offset. As a direct
implementation check, this \(\Pi,\Phi\)-constructed derivative agrees with a
time derivative of the independently interpolated metric on the moving
spheres to about \(2\times10^{-4}\) at the smallest radius.

Fitting the first-order response matrix to \(D_Tg^{ab}\) confirms the
small physical rate scale. At \(R=0.20M_{\rm tot}\), after converting to
centered harmonic time,

\[
\begin{array}{c|ccc}
 & \Pi,\Phi\ {\rm direct} & {\rm finite\ difference}
 & {\rm free\ second\ order}\\ \hline
\operatorname{RMS}|\ddot q^0|
&3.95\times10^{-5}&4.42\times10^{-5}&4.08\times10^{-2}\\
\operatorname{RMS}_{ij}|M_B\dot\Lambda_{ij}|
&1.25\times10^{-5}&1.96\times10^{-5}&7.16\times10^{-3}
\end{array}
\]

Thus the direct derivative independently confirms that the free
second-order residual fit overestimates these rates by roughly three orders
of magnitude.

The direct derivative also exposes an important omission in the provisional
seven-parameter subsystem. The true vector-sector rates
\(\dot\beta_i\) and \(\ddot q^i\) have norms of a few
\(10^{-3}\), much larger than the true clock acceleration and strain rate.
They cannot be omitted while allowing the remaining seven rate columns to
fit freely.

The minimal-estimator comparison gives:

- V1 recovers \((\dot\beta_i,\ddot q^i)\) within about \(5\%\) of the
  finite-difference rates around \(R=0.28M_{\rm tot}\);
- the clock/strain sector is more estimator-dependent;
- C4 agrees with the finite-difference clock/strain rates to about \(12\%\)
  near \(R=0.48M_{\rm tot}\), but has only moderate held-out closure;
- C3 has the best clock/strain held-out closure, about \(18\%\) near
  \(R=0.58M_{\rm tot}\), with a \(25\%\) rate discrepancy.

So the direct metric derivative is already a strong measurement of all
thirteen evolution drives, particularly the vector sector. The
clock/strain rates should initially be treated as an estimator ensemble or
joint constrained fit rather than selected from a single block.

### Integrating back to the first-order velocities

The direct derivative can also be used to cross-check the first-order
velocities themselves. It does not determine their absolute values:
\(\Pi,\Phi\) measure \(d\dot q^\mu/dT\), so integrating requires one initial
value for every component. This is nevertheless an independent test of the
subsequent time dependence.

At \(R=0.20M_{\rm tot}\), anchoring only the first slice and integrating the
direct rates gives

\[
\frac{\|\Delta\dot q^0_{\rm fit}
      -\Delta\dot q^0_{\Pi,\Phi}\|}
     {\max(\|\Delta\dot q^0_{\rm fit}\|,
           \|\Delta\dot q^0_{\Pi,\Phi}\|)}
=10.4\%,
\qquad
\frac{\|\Delta\dot{\boldsymbol q}_{\rm fit}
      -\Delta\dot{\boldsymbol q}_{\Pi,\Phi}\|}
     {\max(\|\Delta\dot{\boldsymbol q}_{\rm fit}\|,
           \|\Delta\dot{\boldsymbol q}_{\Pi,\Phi}\|)}
=1.36\%.
\]

There is a genuinely absolute cross-check for the spatial velocity. With
\(T=t+q^0(t)\),

\[
\dot q^i=(1+\dot q^0)\frac{dz_{\rm center}^i}{dT}.
\]

The curvature-center velocity converted in this way agrees with the
algebraic metric fit to \(2.32\%\) at \(R=0.20M_{\rm tot}\), with a best
agreement of \(2.11\%\) at \(R=0.28M_{\rm tot}\). Omitting the time
conversion gives an error near \(10\%\), so the factor is resolved by the
data rather than being numerically negligible. Anchoring the integrated
spatial acceleration to the center velocity, rather than to the fitted
\(\dot q^i\), still reproduces the full fitted history to \(2.64\%\) at the
smallest radius.

Thus \(\dot q^i\) now has two independent checks: its absolute value from
center motion and its evolution from the direct metric derivative.
\(\dot q^0\) has an independent evolution check but still needs one
algebraic initial calibration; a time derivative alone cannot determine
its additive constant.

## Files

- `q8_temporal_ode_consistency.py`: full data extraction and fitting analysis
- `q8_temporal_ode_consistency.npz`: fitted coefficients and diagnostics
- `q8_temporal_ode_consistency_by_radius.png`: amplitude ratios and tensor
  alignments versus radius
- `q8_temporal_ode_best_timeseries.png`: representative time series
- `q8_direct_metric_time_derivative.py`: direct per-slice \(\Pi,\Phi\)
  derivative construction and rate fit
- `q8_direct_metric_time_derivative.npz`: interpolated derivatives, all
  thirteen recovered rates, and validation data
- `q8_direct_metric_time_derivative.png`: direct, finite-difference, and free
  second-order rate comparison
- `q8_direct_derivative_estimators.py`: all minimal block estimators applied
  to the metric time derivative
- `q8_direct_derivative_estimators.png`: rate agreement and held-out closure
  for V1--V5 and C1--C5
- `q8_first_order_kinematic_reconstruction.py`: integrates the direct rates
  and compares the spatial velocity with curvature-center motion
- `q8_first_order_kinematic_reconstruction.npz`: reconstructed histories and
  radius-dependent mismatches
- `q8_first_order_kinematic_reconstruction.png`: time-series and
  radius-dependent first-order velocity cross-checks
