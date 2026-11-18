# Changelog

## Certificate-only cleanup — 19 November 2026

- Remove the manuscript directory, PDF/LaTeX files, paper-specific license and paper citation metadata from all published history and software snapshots.
- Delete the manuscript-specific tag; retain v1.0.0 and v1.1.0 as certificate-only software tags.
- Replace release archives and checksum files with certificate-only packages, and remove standalone manuscript assets.
- Preserve all certificate source bytes, dependency pins and mathematical execution records.

## 1.1.0 — 15 November 2026

- Add the F4/F5 high-density interval programs, Arb constant and chain-rule checks, integer-series audits, and the 14-suite source/method matrix.
- Add an isolated `reproduce.py` runner with source verification, per-suite timeouts, UTF-8 output, dependency detection, fresh logs, source hashes, and machine-readable results. A skipped, missing, failed, or timed-out suite cannot count as a successful full run.
- Add a new local 14/14 execution and 37/37 auxiliary F3 execution, separate from the imported incomplete F5 record (6 PASS, 8 SKIPPED).
- Add release-specific provenance, distribution and source checksums, and bilingual reproduction guidance.
- Update the software version in `CITATION.cff` to 1.1.0.

The four historical certificate scripts, existing mathematical result files, dependency pins and software license retain their exact bytes. This distribution contains certificate software and documentation.

## 1.0.0 — 13 November 2026

Original public certificate-software release with four calculation programs and their execution records.
