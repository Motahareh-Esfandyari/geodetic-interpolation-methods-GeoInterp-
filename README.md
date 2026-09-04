# GeoInterp — Numerical Interpolation & Approximation Methods in Geodesy

A MATLAB collection of classical interpolation, approximation, and geostatistical
methods, with two real geodetic applications: **GPS satellite orbit interpolation
from precise (SP3) ephemerides** and **strain (dilatation) analysis of a geodetic
network**.

**Author:** Motahareh Esfandyari-Kaloukan

---

## Repository structure

```
GeoInterp/
├── 1-polynomial/            # 1D polynomial interpolation & approximation
│   ├── newton.m             # Newton divided differences
│   ├── lagrange.m           # Lagrange basis polynomials
│   ├── chebyshev.m          # Chebyshev polynomial approximation
│   ├── LagrangeDerivative.m # Evaluation of a tabulated function via Lagrange coefficients
│   └── CubicSpline.m        # Natural cubic spline (function)
│
├── 2-scattered-data/        # 2D scattered-data interpolation
│   ├── IDW.m                # Inverse Distance Weighting + support-domain study
│   ├── RBF.m                # Radial Basis Functions (Gaussian / MQ / IMQ) with derivatives
│   └── RMLS.m               # Recursive Moving Least Squares (quadratic basis)
│
├── 3-geostatistics/         # Kriging & variography
│   ├── KrigingDemo.m        # Ordinary kriging demo on a correlated random field
│   ├── kriging.m            # Ordinary kriging          [third-party, W. Schwanghart]
│   └── variogram.m          # Experimental variogram    [third-party, W. Schwanghart]
│
├── 4-applications/          # Geodetic applications
│   ├── OrbitInterpolation.m    # GPS orbit interpolation from SP3 (cubic spline)
│   ├── OrbitInterpLagrange.m   # Same, Lagrange variant (work in progress)
│   ├── Lagrange0.m             # Lagrange routine for the orbit variant (incomplete)
│   ├── sp3Cread.m              # SP3-a ephemeris reader  [original: LaQ, TU Delft]
│   ├── igs13730.sp3            # Sample IGS precise ephemeris file
│   └── DILITATION.m            # Dilatation of a 2D displacement field via RBF
│
└── 5-special-functions/
    ├── clpn.m               # Legendre polynomials Pn(z), complex argument [Jin & Zhang]
    └── Pnm.m                # Driver for clpn.m
```

## Highlights

- **Orbit interpolation** (`OrbitInterpolation.m`): reads a precise IGS ephemeris
  (15-minute sampling), rotates ECEF positions to an inertial frame using the
  Earth rotation rate, and interpolates the satellite position at any epoch with
  cubic splines. Plots the 3D orbit and the X/Y/Z components versus time.
- **Deformation analysis** (`DILITATION.m`): interpolates the displacement field
  of a scattered station network with RBFs (using their analytical derivatives)
  and compares the numerical dilatation (1/2)(du/dx + dv/dy) with the
  analytical solution.
- **Support-domain studies** (`IDW.m`, `RMLS.m`): how the interpolation error
  and coefficient convergence behave as the number of neighbours (k-NN support
  domain) grows.

## Requirements

- MATLAB (Symbolic Math Toolbox for `newton.m`, `lagrange.m`, `chebyshev.m`;
  Statistics Toolbox for `knnsearch`; Image Processing Toolbox for the kriging
  demo's `imfilter`/`fspecial`)
- [`variogramfit.m`](https://www.mathworks.com/matlabcentral/fileexchange/25948-variogramfit)
  from the MATLAB File Exchange (needed by `KrigingDemo.m`; not included here)

## Usage examples

```matlab
% GPS orbit interpolation at a given epoch
OrbitInterpolation        % prompts for satellite PRN and time [h m s]

% Cubic spline at a query point
yInt = CubicSpline(x, y, x0);

% Local RBF interpolation with derivatives (Gaussian kernel)
[u0, du_dx, du_dy] = RBF(x, y, u, x0, y0, 10, 10, 'G');

% Ordinary kriging demo (needs variogramfit.m)
KrigingDemo
```

## Status / known issues

- `OrbitInterpLagrange.m` and `Lagrange0.m` are work in progress: the Lagrange
  routine is unfinished and the plotting section has an inconsistent variable
  name (`yInt` vs `pInt`).
- `Pnm.m` still contains two leftover lines from the original Fortran
  conversion that must be removed before running.

## Credits

- `kriging.m` and `variogram.m` — Wolfgang Schwanghart (MATLAB File Exchange),
  included unmodified as dependencies; full documentation and license terms in
  the file headers.
- `sp3Cread.m` — original implementation by LaQ, MGP – TU Delft (2003); adapted
  here for the orbit-interpolation application.
- `clpn.m` — MATLAB conversion of the CLPN routine from Jin & Zhang,
  *Computation of Special Functions*.

All remaining scripts and functions were written by
**Motahareh Esfandyari-Kaloukan**.
