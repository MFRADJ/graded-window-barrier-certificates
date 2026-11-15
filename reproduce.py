#!/usr/bin/env python3
"""Run the v1.1.0 software checks in a fresh workspace; never edit source inputs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parent
INPUT_MANIFEST = ROOT / "certificates/SOURCE_SHA256SUMS.txt"
# The four v1.0.0 scripts retain their historical paths and exact bytes.
STANDARD = [
    ("fixed_iv_3regimes", "supplement/certify_mpmath_iv_crosscheck.py", "mpmath", "interval"),
    ("fixed_interval_table", "supplement/regenerate_interval_table.py", "mpmath", "interval"),
    ("fixed_series_exact", "supplement/series_exact.py", None, "exact_integer"),
    ("highdensity_iv_32", "certificates/high_density/certify_highdensity_mpmath_iv_S5_4.py", "mpmath", "interval"),
    ("highdensity_15_small", "tests/ETH_alphaC_V2-S5_2_check_independent.py", None, "integer_coefficients_plus_floating_comparisons_small_cases"),
    ("highdensity_m20_h300", "tests/ETH_alphaC_V2-S5_2_check_threshold.py", None, "integer_coefficients_plus_floating_comparisons"),
]
ARB = [
    ("fixed_density_arb_256", "certify_arb_independent.py", "flint", "Arb"),
    ("highdensity_s2bis_arb_200", "certificates/high_density/certify_constants_arb_V2-S2-bis.py", "flint", "mixed_31_Arb_3_numeric"),
    ("highdensity_s4_arb_256", "certificates/review/certify_constants_V2-S4.py", "flint", "Arb"),
    ("highdensity_chainrule_arb_256", "certificates/review/arb_chain_rule_C1_V2-S5_5.py", "flint", "Arb"),
]
HEAVY = [
    ("s2_old_full_audit", "certificates/high_density/verify_v2s2.py", "flint", "historical_computational_audit"),
    ("s4_exact_6118_reviewer", "certificates/review/audit_exact_V2-S4.py", "flint", "integer_coefficients_plus_floating_log_comparisons"),
    ("s4_large_m", "certificates/review/audit_large_m_V2-S4.py", "flint", "large_m_numeric_audit"),
    ("threshold_mh4000_26", "certificates/review/exact_hd_mh4000_V2-S5_5.py", "flint", "integer_coefficients_plus_floating_log_comparisons"),
]
F3 = [("f3_auxiliary_37", "tests/F3_CONTROLES_INDEPENDANTS.py", "mpmath", "auxiliary_numeric_spot_checks")]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(manifest, root):
    """Return verified paths. Missing, altered or escaping paths are errors."""
    paths = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        expected, relative = line.split("  ", 1)
        path = (root / relative).resolve()
        if not path.is_relative_to(root.resolve()) or not path.is_file():
            raise ValueError(f"Missing or unsafe input: {relative}")
        if sha256(path) != expected:
            raise ValueError(f"SHA-256 mismatch: {relative}")
        paths.append(relative)
    if not paths or len(paths) != len(set(paths)):
        raise ValueError("Empty manifest or duplicate paths")
    return paths


def environment():
    info = {"python": platform.python_version(), "platform": platform.platform()}
    for name, key in (("flint", "python_flint"), ("mpmath", "mpmath")):
        if importlib.util.find_spec(name) is not None:
            module = __import__(name)
            info[key] = getattr(module, "__version__", "unknown")
            if name == "flint":
                info["flint"] = getattr(module, "__FLINT_VERSION__", "not exposed")
        else:
            info[key] = None
    return info


def save_json(path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("verify", "standard", "arb", "all", "f3"), default="all")
    parser.add_argument("--timeout", type=int, default=3600, help="positive time limit per suite, in seconds")
    parser.add_argument("--output", type=Path, default=ROOT / "results/local_runs")
    args = parser.parse_args(argv)
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if sys.flags.optimize:
        parser.error("Run without -O/PYTHONOPTIMIZE; historical checks use assertions")
    try:
        inputs = verify_manifest(INPUT_MANIFEST, ROOT)
        if args.mode == "verify":
            distributed = verify_manifest(ROOT / "SHA256SUMS.txt", ROOT)
            print(f"SOURCE INTEGRITY: PASS ({len(inputs)} inputs)")
            print(f"DISTRIBUTION INTEGRITY: PASS ({len(distributed)} files)")
            return 0
    except (OSError, ValueError) as exc:
        print(f"INTEGRITY ERROR: {exc}", file=sys.stderr)
        return 3
    tasks = {"standard": STANDARD, "arb": ARB, "all": STANDARD + ARB + HEAVY, "f3": F3}[args.mode]
    missing = {task[1] for task in tasks} - set(inputs)
    if missing:
        print(f"INTEGRITY ERROR: unmanifested sources: {sorted(missing)}", file=sys.stderr)
        return 3
    output = args.output.resolve()
    # Allow the default ignored output folder, or a folder outside the repository.
    allowed = (ROOT / "results/local_runs").resolve()
    if output.is_relative_to(ROOT) and not output.is_relative_to(allowed):
        parser.error("Inside the repository, --output must be under results/local_runs")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    run = output / f"{stamp}_{args.mode}_{uuid.uuid4().hex[:8]}"
    work = run / "sources"
    logs = run / "logs"
    logs.mkdir(parents=True)
    for relative in inputs:
        dest = work / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, dest)
    # Verify the copied bytes as well as the distributed inputs.
    verify_manifest(INPUT_MANIFEST, work)
    env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8", PYTHONOPTIMIZE="")
    results = []
    for name, source, dependency, method in tasks:
        record = {"name": name, "source": source, "source_sha256": sha256(work / source), "method": method}
        if dependency and importlib.util.find_spec(dependency) is None:
            record.update(status="SKIPPED", reason=f"missing dependency {dependency}")
        else:
            command = [sys.executable, "-X", "utf8", str(work / source)]
            if name == "fixed_interval_table":
                command += ["--output", str(run / "interval_table.csv")]
            start = time.monotonic()
            log = logs / f"{name}.log"
            with log.open("wb") as stream:
                try:
                    process = subprocess.run(command, cwd=run, env=env, stdout=stream,
                                             stderr=subprocess.STDOUT, timeout=args.timeout)
                    record.update(status="PASS" if process.returncode == 0 else "FAIL", exit_code=process.returncode)
                except subprocess.TimeoutExpired:
                    record.update(status="TIMEOUT", reason=f"exceeded {args.timeout} seconds")
                except OSError as exc:
                    record.update(status="FAIL", reason=str(exc))
            record.update(elapsed_s=round(time.monotonic() - start, 3), log=f"logs/{name}.log", log_sha256=sha256(log))
            if name == "f3_auxiliary_37" and record["status"] == "PASS":
                try:
                    auxiliary = json.loads((work / "tests/F3_CONTROLES_RESULTATS.json").read_text(encoding="utf-8"))
                    if (auxiliary["suite_count"], auxiliary["pass_count"], auxiliary["fail_count"]) != (37, 37, 0):
                        raise ValueError("Expected 37/37 auxiliary checks")
                    record["auxiliary_checks"] = 37
                except (OSError, ValueError, KeyError) as exc:
                    record.update(status="FAIL", reason=str(exc))
        results.append(record)
        print(name, record["status"], flush=True)
    counts = {status: sum(r["status"] == status for r in results)
              for status in ("PASS", "SKIPPED", "FAIL", "MISSING", "TIMEOUT")}
    ok = len(results) == len(tasks) and counts["PASS"] == len(tasks)
    payload = {"software_version": "1.1.0", "mode": args.mode,
               "runner_sha256": sha256(Path(__file__)), "input_manifest_sha256": sha256(INPUT_MANIFEST),
               "started_utc": stamp, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "environment": environment(), "expected": len(tasks), "observed": len(results),
               "all_pass": ok, "status_count": counts,
               "note": "PASS means successful script execution, not human review or a complete formal proof. SKIPPED is never PASS.",
               "results": results}
    save_json(run / "CERTIFICATION_RESULT.json", payload)
    # Input copies may contain generated JSON; hash every artifact after execution.
    manifest = [f"{sha256(path)}  {path.relative_to(run).as_posix()}"
                for path in sorted(run.rglob("*")) if path.is_file()]
    (run / "SHA256_EXECUTION.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print("COUNTS:", counts)
    print("EXECUTION:", run)
    print("RESULT:", "PASS" if ok else "INCOMPLETE/FAIL")
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
