# Graded-window barrier certificates

## Project overview

This repository accompanies Mohamed Fradj's preprint **Thermodynamic saturation and window barriers for graded spaces in a subexponential 3-SAT algorithm**, scientific version **1.0**. The accompanying certificate software is version **1.0.0**. It studies the dimensions of prescribed graded windows motivated by the Appendix B construction in Anthropic's 3-SAT algorithm. The graded spaces arise from the vacuum and nonzero-coset characters of D4 and additive isotropic codes over the Klein four-group.

The repository provides reproducible numerical evidence for the three sign conditions in Appendix A, exact checks of finite coefficient identities, and the data displayed in its interval table. The manuscript supplies the analytic arguments connecting these computations to the window bounds.

## Main mathematical result

Fix a compact density interval $J=[x_1,x_2]\subset(0,\infty)$. Theorem 6.5 gives constants $\eta_J>0$ and $C_J>0$, depending only on $J$ and the fixed character alphabet, such that

$$
\rho_{\mathrm{win}}(V_S;h,\lambda)
\geq \frac{1}{\sqrt{2}}+\eta_J-C_Jm^{-1/2}.
$$

The quantifiers are uniform over integers $m\geq1$, additive isotropic codes $S\subset(\mathbb F_2^2)^m$, admissible positive levels $h\in\frac12\mathbb Z$ with $h/m\in J$, and $\lambda\geq1$ satisfying $2^\lambda\leq\dim(V_S)_h$. The grading and the real windows are those defined in the manuscript. Admissibility means that the level is present. Smaller lengths are absorbed by enlarging $C_J$; the strict-margin argument first applies beyond a threshold depending on $J$. For integer-degree windows, use the even-weight subcode as described in the manuscript.

The theorem states the displayed strict fixed-density bound on each fixed compact interval. The scripts do not compute an optimal value of $\eta_J$, and the estimate is not uniform as $h/m$ approaches zero or infinity.

This is a bound on the **cost of prescribed windows**. It does not establish a universal lower bound on intrinsic matrix or factorization ranks.

## What this repository certifies

The canonical Arb program uses 256-bit ball arithmetic, exact integer coefficients, and explicit analytic truncation majorants. It verifies the numerical conditions used for:

1. $(0,4]$: the modular point bound $R(\pi^2)<4\pi^2/3$;
2. $[4,5]$: negativity of the equation (A.3) upper bound on ten rational boxes;
3. $[5,\infty)$: positivity of the large-inverse-temperature bracket at $\beta=5$.

The point tests imply the stated domains through the analytic monotonicity and tail arguments in Appendix A. The program does not sample or directly enumerate an infinite domain. Every failed sign condition raises an error or exits nonzero, and coefficient invariants remain active under Python optimization.

The separate `mpmath.iv` program uses a different interval implementation with the same analytic tail majorants. This is an implementation cross-check, not an independent proof of those majorants or independent human review. The exact checker verifies identities for every retained half-step exponent strictly below 200; the infinite identities also use the classical arguments in the manuscript. No script is a Lean proof of the complete theorem.

## Repository structure

```text
graded-window-barrier-certificates/
├── README.md
├── LICENSE
├── LICENSE-paper
├── CITATION.cff
├── requirements.txt
├── .gitignore
├── SHA256SUMS.txt
├── certify_arb_independent.py
├── supplement/
│   ├── series_exact.py
│   ├── certify_mpmath_iv_crosscheck.py
│   └── regenerate_interval_table.py
├── results/
│   ├── arb_execution_final.log
│   ├── series_exact.log
│   ├── mpmath_iv_crosscheck.log
│   └── interval_table.csv
└── manuscript/
    ├── Fradj_Thermodynamic_Saturation_and_Window_Barriers_v1.0.tex
    └── Fradj_Thermodynamic_Saturation_and_Window_Barriers_v1.0.pdf
```

| Script | Scientific role |
| --- | --- |
| [`certify_arb_independent.py`](certify_arb_independent.py) | Canonical Arb certificate for the three Appendix A numerical conditions. |
| [`supplement/series_exact.py`](supplement/series_exact.py) | Exact integer checks of the Jacobi quartic identity, parity, initial coefficients, positivity, character relation, and finite triple products. |
| [`supplement/certify_mpmath_iv_crosscheck.py`](supplement/certify_mpmath_iv_crosscheck.py) | Separate 80-decimal-digit interval cross-check of the same three numerical conditions. |
| [`supplement/regenerate_interval_table.py`](supplement/regenerate_interval_table.py) | Regenerates the eleven-row CSV after running the interval checks. Decimal formatting uses exact rational representations of interval endpoints. |

The CSV includes ten-decimal approximations for $A$, $\Delta_A$, and $\mathcal C$, four-decimal upper bounds for $U$, and outward-rounded 30-decimal enclosures. Its last row has no outgoing box and therefore has empty $U$ fields. A separate Arb calculation checked every displayed entry and exported enclosure against the manuscript table. The table generator itself uses `mpmath.iv`.

## Requirements

The recorded executions used:

| Component | Version or setting |
| --- | --- |
| CPython | 3.14.3, 64-bit Windows |
| python-flint | 0.9.0 |
| FLINT | 3.6.0 |
| Arb | Integrated into FLINT; no separate Arb version is exposed by this installation |
| mpmath | 1.3.0; `iv.dps = 80` |
| Arb working precision | 256 bits |
| Manuscript compiler | Tectonic 0.17.0 |

Install the pinned Python packages in [`requirements.txt`](requirements.txt). The exact checker needs only the Python standard library. The commands below also apply on Linux/macOS, but these operating systems have not been executed as part of the recorded validation.

## Reproduction instructions

Clone the versioned source and run from its root:

```sh
git -c core.autocrlf=false clone --branch manuscript-v1.0 --depth 1 https://github.com/MFRADJ/graded-window-barrier-certificates.git
cd graded-window-barrier-certificates
```

On Linux/macOS:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip check
python certify_arb_independent.py
python supplement/series_exact.py
python supplement/certify_mpmath_iv_crosscheck.py
python supplement/regenerate_interval_table.py --output interval_table.local.csv
```

On Windows PowerShell:

```powershell
py -3.14 -m venv .venv
$python = '.\.venv\Scripts\python.exe'
& $python -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed' }
& $python -m pip check
if ($LASTEXITCODE -ne 0) { throw 'Dependency check failed' }
foreach ($script in @('certify_arb_independent.py', 'supplement/series_exact.py', 'supplement/certify_mpmath_iv_crosscheck.py')) {
    & $python $script
    if ($LASTEXITCODE -ne 0) { throw "Verification failed: $script" }
}
& $python supplement/regenerate_interval_table.py --output interval_table.local.csv
if ($LASTEXITCODE -ne 0) { throw 'Table generation failed' }
```

Compare the regenerated table byte for byte, on either platform:

```sh
python -c "from pathlib import Path; assert Path('interval_table.local.csv').read_bytes() == Path('results/interval_table.csv').read_bytes(); print('TABLE MATCH: PASS')"
```

The generator defaults to `results/interval_table.csv` if `--output` is omitted. Recorded logs include UTC start/end times, dependency versions, source SHA-256 values, stdout, stderr, and exit codes. Re-execution produces the same mathematical indicators; timestamps in a new execution record will differ.

Verify the complete distributed file manifest before changing any file (cross-platform):

```sh
python -c "from pathlib import Path; import hashlib; entries=[line.split('  ',1) for line in Path('SHA256SUMS.txt').read_text(encoding='utf-8').splitlines() if line]; assert all(hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest for digest,name in entries); print('SHA-256: PASS')"
```

On Linux, `sha256sum -c SHA256SUMS.txt` is an alternative; on macOS use `shasum -a 256 -c SHA256SUMS.txt`. Text files are committed with LF line endings. Clone with `git -c core.autocrlf=false clone ...` if a local Git configuration would otherwise rewrite line endings and prevent byte-exact verification.

## Expected results

Each program must return exit status zero. The success markers are:

```text
==> APPENDIX A: CERTIFIED
==> EXACT SERIES CHECKS: PASS
==> MAIN (0,4]/[4,5]/[5,infty) INTERVAL CROSS-CHECK: PASS
==> INTERVAL TABLE: PASS (11 rows; exact endpoint rounding)
```

The Arb log reports approximately $R(\pi^2)=10.2925020821$, a largest box upper endpoint of $-0.161844090586$, and a positive large-inverse-temperature bracket of $42.9378683703$. These decimals summarize interval enclosures; the decisions use the full intervals.

## Manuscript

The public manuscript is **v1.0** (13 November 2026). The filename and the version printed in the manuscript agree. The `manuscript-v1.0` tag records the synchronized manuscript and documentation; `v1.0.0` remains the original certificate-software tag. The assets attached to the software release include the current manuscript and an updated reproducibility archive.

Read the [preprint PDF](manuscript/Fradj_Thermodynamic_Saturation_and_Window_Barriers_v1.0.pdf) or its [LaTeX source](manuscript/Fradj_Thermodynamic_Saturation_and_Window_Barriers_v1.0.tex). The PDF is compiled from this source. With [Tectonic](https://tectonic-typesetting.github.io/en-US/) 0.17.0 installed:

```sh
tectonic --untrusted manuscript/Fradj_Thermodynamic_Saturation_and_Window_Barriers_v1.0.tex
```

Tectonic obtains the required LaTeX packages on first use. A conventional LaTeX distribution with the packages declared in the source can also compile it; that alternative build has not been validated here. Rebuilding may change PDF metadata or typesetting-library output, so byte identity of rebuilt PDFs is not asserted. The distributed PDF's exact bytes are covered by the manifest.

## Limitations

This work neither proves nor refutes ETH, establishes a universal bound on intrinsic ranks, nor constructs the missing algebraic operators. It concerns the specified windows at fixed density. Analytic inequalities, tail majorants, reductions, and hypotheses must be assessed in the manuscript; software execution alone does not certify the entire theorem. Independent qualified human review remains outstanding, and no peer-reviewed publication status is claimed.

## Citation

Use [`CITATION.cff`](CITATION.cff) for machine-readable metadata. Cite the manuscript as:

> Mohamed Fradj. *Thermodynamic saturation and window barriers for graded spaces in a subexponential 3-SAT algorithm*. v1.0, preprint, 2026.

Cite the accompanying software as version **1.0.0**, *Graded-window barrier certificates*, Mohamed Fradj, at the [versioned release](https://github.com/MFRADJ/graded-window-barrier-certificates/releases/tag/v1.0.0). The [repository](https://github.com/MFRADJ/graded-window-barrier-certificates) contains the versioned sources. No manuscript or software DOI has been assigned.

For Arb, cite Fredrik Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Transactions on Computers 66 (2017), 1281-1292, [doi:10.1109/TC.2017.2690633](https://doi.org/10.1109/TC.2017.2690633).

## License

The software and repository documentation are licensed under the [MIT License](LICENSE), copyright Mohamed Fradj. The manuscript source and PDF are **all rights reserved**, as specified in [`LICENSE-paper`](LICENSE-paper). The manuscript has no Creative Commons license. Installed dependencies retain their own licenses and are not bundled with this repository.
