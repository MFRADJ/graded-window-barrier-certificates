# Fresh software-v1.1.0 executions — 15 November 2026

| Execution | Result |
| --- | --- |
| Original unmodified F4 runner, all suites | 14 PASS, 0 SKIPPED, 0 FAIL, 0 TIMEOUT |
| Distributed v1.1.0 runner, all suites | 14 PASS, 0 SKIPPED, 0 FAIL, 0 TIMEOUT |
| Distributed auxiliary F3 checker | 37/37 PASS (one auxiliary suite) |
| Negative execution checks | Missing dependency, timeout, changed source, optimized Python and protected output correctly rejected |

Environment: CPython 3.14.3, Windows 64-bit, python-flint 0.9.0, FLINT 3.6.0, mpmath 1.3.0. `pip check` found no broken requirements. The summaries give UTC times and the full platform string.

`all/` and `f3/` contain the new runner's summaries, stdout/stderr logs and generated data. Each result identifies its source SHA-256 and log SHA-256. The F3 JSON contains all 37 check results. The generated integer-audit TSV is included. `original_f4/` retains the separate first replay of the original archive's software runner; it is an execution of the original F4 code, not a manuscript compilation. Its relative reproduction-log paths are preserved.

7 local path-prefix occurrences were replaced in public evidence with `<REPOSITORY_ROOT>`, `<LOCAL_WORKSPACE>` or `<PYTHON>`. These replacements affect launcher metadata and generated-file locations, not mathematical output. `SHA256SUMS.txt` in this directory covers the published evidence bytes. For each new runner log that contains such a prefix, the summary preserves `raw_log_sha256`, records `local_path_prefix_redacted: true`, and sets `log_sha256` to the published bytes. Every published log digest is checked during packaging.

The delivered F5 6/14 result remains separate in `../f5_imported_2026-11-15/`. No skipped suite is counted as a success. These new runs do not supply independent qualified human review or a formal proof of the complete article. Manuscript files are not part of this certificate-only distribution.
