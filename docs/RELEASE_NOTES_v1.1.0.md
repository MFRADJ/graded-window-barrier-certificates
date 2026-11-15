# Certificate software v1.1.0 — 14-suite reproduction, certificates only

Software v1.1.0 adds high-density F4/F5 programs and an isolated reproduction runner. The recorded Windows execution using CPython 3.14.3, python-flint 0.9.0 / FLINT 3.6.0 and mpmath 1.3.0 passed **14/14 suites, 0 skipped**. The separate F3 checker passed **37/37 auxiliary checks**.

The package includes source hashes, a suite/method matrix, individual logs, machine-readable summaries, SHA-256 manifests, and English/French reproduction instructions. Every new reproduction uses a fresh workspace. Altered inputs, missing dependencies, timeouts and program failures cannot produce a successful complete run.

The four original certificate programs, additional Python inputs and mathematical execution records retain their original bytes. The original F5 result remains recorded separately as 6 PASS / 8 SKIPPED.

On 19 November 2026, at the author's request, manuscript PDF/LaTeX files were removed from the published history, tags and both release packages. The manuscript-specific tag was deleted. This archive replaces the earlier v1.1.0 package; archive and distribution checksums have changed. The software version remains 1.1.0 because the calculation programs are unchanged.

```sh
python -m pip install -r requirements.txt
python -m pip check
python reproduce.py --mode verify
python reproduce.py --mode all --timeout 3600
python reproduce.py --mode f3
```

Use Python without optimization. Windows is the recorded validation platform; Linux/macOS instructions are provided without claiming a new execution there.

Assets: `reproducibility_v1_1_0.zip`, `reproducibility_v1_1_0.zip.sha256`, and `SHA256SUMS_v1_1_0.txt`. No manuscript asset is distributed. Software/documentation: MIT.

Successful numerical execution is not independent qualified human review or a formal proof of complete theorems. S2-bis mixes 31 Arb checks and three numerical checks; several integer-coefficient audits use floating logarithmic comparisons. No claim to prove or refute ETH is made.
