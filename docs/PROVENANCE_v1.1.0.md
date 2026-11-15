# Software v1.1.0 provenance

This certificate software was released on 15 November 2026. This distribution includes the calculation programs, evidence and documentation. Manuscript PDF and LaTeX files are distributed separately.

## Delivered inputs

| Input | SHA-256 |
| --- | --- |
| `ETH_alphaC_F5_DOSSIER_VALIDATION_CONDITIONNELLE.zip` | `f2b39d2b924efcc0e195921183336fc087109d38ae1423eb19e87af936c4c1bc` |
| Nested `F4_ARCHIVE_ORIGINALE_NON_MODIFIEE.zip` | `0f1ac488742e0b7316456249a2ffc75e0d97f9c02371677d454a30084b744c47` |

All 27 F5 manifest entries and all 135 F4 manifest entries were checked before selecting the software files. The imported Python source bytes are unchanged. `certificates/SOURCE_SHA256SUMS.txt` records all 16 immutable Python inputs, including the four historical programs and the auxiliary F3 checker. `certificates/CERTIFICATE_MATRIX_14.csv` maps the 14 suites to their methods, source hashes, and historical manuscript references. Those references identify the F4/F5 development context; they do not rename or revise the existing public manuscript v1.0.

The source-frozen software license in F4 is MIT, copyright Mohamed Fradj, consistent with the existing repository `LICENSE`. Dependency packages are installed separately and retain their own licenses.

## Three distinct sets of evidence

1. `results/arb_execution_final.log`, `series_exact.log`, `mpmath_iv_crosscheck.log`, and `interval_table.csv` are the existing v1.0 records, retained byte for byte.
2. `results/f5_imported_2026-11-15/` contains the delivered F5 summary and matrix: **6 PASS, 8 SKIPPED**. The missing dependency in that environment was `python-flint`. These files remain historical records of that environment.
3. `results/v1.1.0_2026-11-15/` contains fresh executions made for this software release on Windows, CPython 3.14.3, python-flint 0.9.0 / FLINT 3.6.0, and mpmath 1.3.0. These runs are distinct from the historical reviewer logs. The new runner also re-executed all 14 suites and the 37 auxiliary F3 checks using isolated copies.

The delivered F4 runner was first executed successfully against its complete original manifest. The published v1.1.0 runner was then exercised against the selected software sources. Source hashes and per-log hashes identify the actual executed bytes. Personal local path prefixes in the published evidence, if any, are replaced with `<REPOSITORY_ROOT>` / `<PYTHON>`; mathematical output is retained. Summaries identify any such replacements.

## Certificate-source preservation

`docs/history/V1_0_0_SCRIPT_BLOBS.json` records the unchanged Git blob IDs and SHA-256 hashes of the four historical certificate programs. All sixteen Python inputs, the runner, dependency pins and original mathematical execution logs are retained byte for byte. The published history and software tags have been rewritten to exclude manuscript files; earlier commit IDs and package digests no longer identify the current distribution.

## Interpretation

A successful suite means that the distributed program completed successfully. The S2-bis suite combines 31 Arb checks and three numerical checks. Several large integer-coefficient audits use floating logarithmic comparisons; they are not wholly rigorous interval certificates. The optional `fixed_density_limit_V2-S4.py` heuristic is included as source but is outside the 14-suite matrix and is not run by `--mode all`.

Numerical execution does not supply independent qualified human mathematical review, prove the complete theorem formally, or prove or refute ETH. No new manuscript license or manuscript publication is asserted by software v1.1.0.
