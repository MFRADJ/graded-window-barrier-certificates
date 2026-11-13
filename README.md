# Graded-window barrier certificates — software v1.0.0

This repository distributes certificate software, numerical execution records and reproduction instructions for specified graded-window inequalities. The graded spaces involve D4 characters and additive isotropic Klein-four codes.

The manuscript PDF and LaTeX source are not distributed in this certificate-only snapshot. Mathematical hypotheses, analytic reductions and tail estimates remain separate proof obligations.

## What this repository certifies
The canonical Arb program uses 256-bit ball arithmetic, exact integer coefficients, and explicit analytic truncation majorants. It verifies the numerical conditions used for:

1. $(0,4]$: the modular point bound $R(\pi^2)<4\pi^2/3$;
2. $[4,5]$: negativity of the equation (A.3) upper bound on ten rational boxes;
3. $[5,\infty)$: positivity of the large-inverse-temperature bracket at $\beta=5$.

The point tests imply the stated domains through the analytic monotonicity and tail arguments in Appendix A. The program does not sample or directly enumerate an infinite domain. Every failed sign condition raises an error or exits nonzero, and coefficient invariants remain active under Python optimization.

The separate `mpmath.iv` program uses a different interval implementation with the same analytic tail majorants. This is an implementation cross-check, not an independent proof of those majorants or independent human review. The exact checker verifies identities for every retained half-step exponent strictly below 200; the infinite identities also use the classical arguments in the manuscript. No script is a Lean proof of the complete theorem.


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

Install the pinned Python packages in [`requirements.txt`](requirements.txt). The exact checker needs only the Python standard library. The commands below also apply on Linux/macOS, but these operating systems have not been executed as part of the recorded validation.

## Reproduction instructions
Clone the versioned source and run from its root:

```sh
git -c core.autocrlf=false clone --branch v1.0.0 --depth 1 https://github.com/MFRADJ/graded-window-barrier-certificates.git
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


## Limitations
This work neither proves nor refutes ETH, establishes a universal bound on intrinsic ranks, nor constructs the missing algebraic operators. It concerns the specified windows at fixed density. Analytic inequalities, tail majorants, reductions, and hypotheses must be assessed in the manuscript; software execution alone does not certify the entire theorem. Independent qualified human review remains outstanding, and no peer-reviewed publication status is claimed.


## Citation

Use [CITATION.cff](CITATION.cff) to cite Mohamed Fradj, *Graded-window barrier certificates*, software version 1.0.0. No DOI is assigned.

For Arb, cite Fredrik Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Transactions on Computers 66 (2017), 1281–1292, [doi:10.1109/TC.2017.2690633](https://doi.org/10.1109/TC.2017.2690633).

## License

Software and repository documentation are licensed under the [MIT License](LICENSE), copyright Mohamed Fradj. Dependencies retain their own licenses and are not bundled.
