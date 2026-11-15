# Numerical reproduction of the canonical predictor

The program computes the canonical maximum from scratch for ten square board sizes. Separately, it solves the equations for the profile of mass one and evaluates the Nyström discretization. The reference integers for the count come from OEIS A000170. The finite-activity reference and the canonical reference for n=5 are kept in separate files.

## Running the package

Node.js 20 or later and Python 3 with NumPy are required. Install the dependencies in your environment with `python -m pip install -r requirements.txt`. From this directory, `node reproducir.mjs` creates a new directory for all outputs and timings. On Unix you can also use `sh reproducir.sh`. Set `PYTHON` to select the interpreter. The calculations use the CPU and double precision.

To check the supplied values without running the solvers again, use `node comprobar.mjs`. In this distribution, the default is the English article in `../paper`. To check the approved Spanish version instead, pass the path to its `articulo` directory. Each language has the same 90 checks. They cover the individual entries in both tables, the fits, both grid parities in the quadrature, the decimal precision in the abstract and introduction, and the empirical checks in the discussion. The process returns a nonzero exit code for a discrepancy. The full reproduction also runs `barrera_4_5.py`, saves its output and checks its two boundary transitions. The tabulated canonical value for n=5 is obtained separately from its reference at mass one.

## Measured running times

The CPU times below were measured in the JavaScript V8 environment used for the Spanish version and rounded to milliseconds. The Python time depends on the interpreter and is written to the timing file when the full reproduction is run.

| Calculation | Cases | Seconds |
|---|---:|---:|
| Canonical maximum, pseudodeterminant and gap | 10 | 3.251 |
| Profile and Nyström discretization | 50 pairs of grids | 19.541 |
| Finite activity | 13 × 61 | 2.506 |
| Alternative barrier solver using symmetries, JavaScript | 3 × 41 | 0.182 |

## Results and files

`canonico.csv` records extrema, residuals, iterations and gaps before rounding. The ten values of log Q_pred agree with the reference table to seven decimal places and the gaps to five. The gap test uses inverse iteration with a residual check and shifted Cholesky factorization on the complement of the zero mode. The checks of positive definiteness and interiority are numerical.

| n | log Q_pred | min p | max p | 1/4 − max p | min h | gap |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 24.3667849 | 0.02473820 | 0.06915949 | 0.18084051 | 0.17854349 | 0.14177800 |
| 27 | 39.9860346 | 0.01787269 | 0.05152626 | 0.19847374 | 0.17749523 | 0.14130756 |
| 40 | 73.4844735 | 0.01177350 | 0.03486465 | 0.21513535 | 0.17636131 | 0.14083134 |
| 60 | 132.9170795 | 0.00771761 | 0.02329709 | 0.22670291 | 0.17560694 | 0.14047948 |
| 80 | 199.0846876 | 0.00573905 | 0.01749172 | 0.23250828 | 0.17522429 | 0.14029676 |
| 100 | 270.2725639 | 0.00456768 | 0.01400209 | 0.23599791 | 0.17499290 | 0.14018495 |
| 130 | 384.3475269 | 0.00349688 | 0.01077685 | 0.23922315 | 0.17477810 | 0.14008030 |
| 160 | 505.3805043 | 0.00283273 | 0.00875917 | 0.24124083 | 0.17464325 | 0.14001418 |
| 200 | 675.3684475 | 0.00226030 | 0.00700937 | 0.24299063 | 0.17452600 | 0.13995644 |
| 240 | 853.3899390 | 0.00188032 | 0.00584225 | 0.24415775 | 0.17444764 | 0.13991772 |

The central reference is α₀=1.944001. The two outer columns below use 1.944000752 and 1.944001082, respectively. The least-squares fits use Householder transformations. The aggregate endpoints printed in the article are rounded.

| Range of n | Degree in 1/n | C at lower α | C at α₀ | C at upper α |
|---|---:|---:|---:|---:|
| 40–130 | 2 | 0.94011554 | 0.94017293 | 0.94019191 |
| 40–130 | 3 | 0.94005476 | 0.94013367 | 0.94015977 |
| 40–240 | 2 | 0.94007420 | 0.94015487 | 0.94018155 |
| 40–240 | 3 | 0.94002363 | 0.94013272 | 0.94016878 |
| 60–240 | 2 | 0.94004855 | 0.94014309 | 0.94017435 |
| 60–240 | 3 | 0.94000433 | 0.94013200 | 0.94017421 |
| 20–240 | 2 | 0.94036043 | 0.94042196 | 0.94044231 |
| 20–240 | 3 | 0.93968175 | 0.93976215 | 0.93978873 |

`ny_resultados.csv` contains the 50 grids. The pairs 192/384 and 193/385 give 0.94014395698 and 0.94014395711, respectively. The nodal profile factors at m=385 are in `perfil_m385.json`. Extrapolating I₂ and log det₂ from the even pair and then exponentiating log det₂ gives I₂≈1.02638534 and det₂≈0.005162694. The difference between the two parities measures numerical stability. The estimate has no certified error bound.

`referencias/anunciados.json` records the value of Q(28) reported in November 2026 (Cherepanov, Zenodo; OEIS A000170). It is excluded from the tables and the figure. Its only use here is to check the informative value of C′₀(28) in Section 9. `actividades_793.csv` documents the activity grid. `referencias/barrera_archivada.txt` contains the historical output of the Python program. `barrera_alternativa.json` records the 123 points from an independent JavaScript implementation: for n=4 and 5, the maximum first reaches the threshold 0.24999 at Λ=100 and Λ=10,000, respectively. For n=6 at Λ=10⁹, the maximum is 0.2165479. The threshold refers to the grid tested. These finite-activity results and the estimated canonical value for n=5 have different interpretations.

`comprobar.mjs` compares the reference values with the LaTeX text of the sections reporting numerical results. The supplied verification uses the archived Python output for the barrier calculation and the recomputed canonical values for Table 1. The full reproduction replaces the finite-activity diagnostic with values from the new Python run.

The code is distributed under the MIT license. See `CITATION.cff` for authorship information.
