# Graded-window barrier certificates — software v1.1.0

[Français](README_FR.md)

Reproducible numerical and integer-series checks for prescribed graded-window barriers. Software **v1.1.0** adds the high-density F4/F5 programs and a unified **14-suite** reproduction runner. A fresh Windows execution passed **14/14 suites**, with no skipped tests; the separate auxiliary F3 checker passed **37/37 checks**. Successful execution is computational evidence, not qualified human review or a formal proof of the complete theorems.

This repository distributes certificate software, test data, execution logs and reproduction documentation. Manuscript PDFs and LaTeX sources are not included. The v1.0.0 and v1.1.0 tags and downloadable packages provide certificate-only snapshots.

## Reproduce software v1.1.0

```sh
git clone --branch v1.1.0 https://github.com/MFRADJ/graded-window-barrier-certificates.git
cd graded-window-barrier-certificates
```

Windows PowerShell:

```powershell
py -3.14 -m venv .venv
$python = '.\.venv\Scripts\python.exe'
& $python -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed' }
& $python -m pip check
if ($LASTEXITCODE -ne 0) { throw 'Dependency check failed' }
foreach ($mode in @('verify', 'all', 'f3')) {
    & $python reproduce.py --mode $mode --timeout 3600
    if ($LASTEXITCODE -ne 0) { throw "Verification failed: $mode" }
}
```

Linux/macOS (instructions provided; these platforms were not executed for this release):

```sh
set -e
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip check
.venv/bin/python reproduce.py --mode verify
.venv/bin/python reproduce.py --mode all --timeout 3600
.venv/bin/python reproduce.py --mode f3
```

The validated environment is **CPython 3.14.3, 64-bit Windows, python-flint 0.9.0, FLINT 3.6.0, mpmath 1.3.0**. Dependencies are pinned in [`requirements.txt`](requirements.txt). Arb is integrated in FLINT. Use Python without `-O` or `PYTHONOPTIMIZE`; the historical programs use assertions.

## Runner behavior

| Mode | Scope | Success condition |
| --- | --- | --- |
| `verify` | SHA-256 source and distribution manifests | Every listed file matches |
| `standard` | Six standard-library / mpmath suites | 6 PASS, no other status |
| `arb` | Four suites requiring python-flint | 4 PASS, no other status |
| `all` | Standard + Arb + four larger audits | 14 PASS, no other status |
| `f3` | Separate auxiliary arithmetic checker | 37/37 auxiliary checks |

Each invocation creates a fresh directory under `results/local_runs/`, copies the verified Python sources, and runs them in isolation. It writes individual logs, an interval CSV where applicable, generated audit artifacts, `CERTIFICATION_RESULT.json`, and `SHA256_EXECUTION.txt`. Existing source files and historical result logs are not rewritten. New runs are ignored by Git. An external output parent can be selected with `--output`; output inside the repository is restricted to `results/local_runs/`.

`--timeout 3600` is a limit in seconds **per program**. `SKIPPED` (missing dependency), `FAIL`, `MISSING`, or `TIMEOUT` cannot produce a successful full-run status. Source corruption stops the runner before execution. Use the exit code and machine-readable summary, not a success marker copied from an old log.

## Sources and evidence

| Location | Contents |
| --- | --- |
| [`certify_arb_independent.py`](certify_arb_independent.py), [`supplement/`](supplement/) | The four original v1.0.0 programs, unchanged |
| [`certificates/high_density/`](certificates/high_density/) | High-density interval and constant checks, historical S2 audit |
| [`certificates/review/`](certificates/review/) | Arb chain rule / constants, exact coefficient and large-m audits |
| [`tests/`](tests/) | Small-case, threshold, and auxiliary F3 checks |
| [`certificates/CERTIFICATE_MATRIX_14.csv`](certificates/CERTIFICATE_MATRIX_14.csv) | Suite IDs, methods, source hashes, historical mathematical context |
| [`results/v1.1.0_2026-11-15/`](results/v1.1.0_2026-11-15/README.md) | Fresh release executions, logs, environment and hashes |
| [`results/f5_imported_2026-11-15/`](results/f5_imported_2026-11-15/) | Original delivered F5 summary: **6 PASS, 8 SKIPPED**, preserved as such |
| [`docs/PROVENANCE_v1.1.0.md`](docs/PROVENANCE_v1.1.0.md) | Archive hashes, source selection, licenses and historical release baseline |

The original records `results/arb_execution_final.log`, `series_exact.log`, `mpmath_iv_crosscheck.log`, and `interval_table.csv` retain their exact bytes. New evidence is stored separately. [`SHA256SUMS.txt`](SHA256SUMS.txt) covers every other distributed file. `.gitattributes` preserves distributed file bytes on checkout, including on Windows.

The four original programs can still be run individually using the v1.0 instructions. For new runs, the unified runner prevents historical outputs from being overwritten. Its fixed-density Arb suite checks the three Appendix A numerical conditions; its separate interval checker uses the same analytic tail majorants with another interval implementation.

## Scientific scope and limits

The checks concern graded spaces constructed from D4 characters and additive isotropic Klein-four codes, at fixed and high density. The matrix's theorem numbers identify the F4/F5 development context. The mathematical text is distributed separately from this repository; software execution alone does not establish its hypotheses or analytic arguments.

The S2-bis suite combines **31 Arb checks and three numerical checks**. Several integer-coefficient audits use floating logarithmic comparisons. The complete set is therefore not described as 14 fully rigorous interval certificates. The exploratory `certificates/review/fixed_density_limit_V2-S4.py` is included as source, outside the 14 suites, and is not run by `--mode all`.

Analytic tail bounds, reductions, and hypotheses remain mathematical proof obligations. These computations do not constitute a Lean proof of the complete theorems, establish universal intrinsic-rank bounds, construct missing algebraic operators, or prove or refute ETH. Independent qualified human mathematical review remains outstanding; no peer-reviewed publication status is claimed.

## Citation and license

Use [`CITATION.cff`](CITATION.cff) to cite software **1.1.0**. No DOI is assigned.

The software and repository documentation use the [MIT license](LICENSE), copyright Mohamed Fradj. Dependencies retain their own licenses and are not bundled. See the [changelog](CHANGELOG.md) for the software changes.
