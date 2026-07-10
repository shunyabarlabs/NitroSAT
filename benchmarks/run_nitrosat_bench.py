#!/usr/bin/env python3
"""
NitroSAT Benchmark Runner

Runs the nitrosatv2 binary on a directory of CNF files,
collects results (from JSON output + timing), and produces a summary.

Usage:
    python3 run_nitrosat_bench.py --cnf-dir ./cnf_suite --solver ../src/c/v2/nitrosatv2 --timeout 300

Outputs:
    - results.csv
    - summary.txt
"""

import argparse
import json
import subprocess
import time
import csv
import os
import sys
from pathlib import Path
from typing import Dict, Any, Optional

def run_solver(solver_path: str, cnf_path: str, timeout: int, exact: bool) -> Dict[str, Any]:
    """Run nitrosat on a CNF and return parsed results + metadata."""
    result = {
        "cnf_file": os.path.basename(cnf_path),
        "full_path": str(cnf_path),
        "status": "ERROR",
        "satisfaction_rate": None,
        "satisfied": None,
        "unsatisfied": None,
        "variables": None,
        "clauses": None,
        "wall_time_s": None,
        "return_code": None,
        "error": None,
        "raw_json": None,
    }

    cmd = [solver_path, cnf_path]
    if exact:
        cmd.append("--exact")

    start = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        wall_time = time.time() - start
        result["wall_time_s"] = round(wall_time, 3)
        result["return_code"] = proc.returncode

        stdout = proc.stdout.strip()
        stderr = proc.stderr.strip()

        # Try to parse JSON from stdout
        try:
            data = json.loads(stdout)
            result["raw_json"] = data
            result["status"] = data.get("status", "UNKNOWN")
            result["satisfied"] = data.get("satisfied")
            result["unsatisfied"] = data.get("unsatisfied")
            result["variables"] = data.get("variables")
            result["clauses"] = data.get("clauses")
            result["satisfaction_rate"] = data.get("satisfaction_rate")
            if result["satisfaction_rate"] is None and result["satisfied"] is not None and result["clauses"] is not None:
                result["satisfaction_rate"] = float(result["satisfied"]) / float(result["clauses"])
        except json.JSONDecodeError:
            # Fallback: try to parse basic info from stderr if JSON failed
            result["error"] = "Failed to parse JSON output"
            if "NitroSAT" in stderr:
                result["status"] = "PARTIAL" if "PARTIAL" in stderr else "ERROR"

    except subprocess.TimeoutExpired:
        result["wall_time_s"] = timeout
        result["status"] = "TIMEOUT"
        result["error"] = f"Exceeded {timeout}s timeout"
    except Exception as e:
        result["error"] = str(e)

    return result


def main():
    parser = argparse.ArgumentParser(description="Benchmark NitroSAT on CNF files")
    parser.add_argument("--cnf-dir", type=str, required=True, help="Directory containing .cnf files")
    parser.add_argument("--solver", type=str, default="../src/c/v2/nitrosatv2", help="Path to nitrosatv2 binary")
    parser.add_argument("--timeout", type=int, default=300, help="Per-instance timeout in seconds")
    parser.add_argument("--output", type=str, default="results.csv", help="Output CSV file")
    parser.add_argument("--exact", action="store_true", help="Run solver with --exact flag")
    args = parser.parse_args()

    cnf_dir = Path(args.cnf_dir).resolve()
    solver = Path(args.solver).resolve()

    if not solver.exists():
        print(f"ERROR: Solver not found at {solver}")
        sys.exit(1)

    cnf_files = sorted(list(cnf_dir.glob("*.cnf")))
    if not cnf_files:
        print(f"No .cnf files found in {cnf_dir}")
        sys.exit(1)

    print(f"Found {len(cnf_files)} CNF files in {cnf_dir}")
    print(f"Using solver: {solver}{' --exact' if args.exact else ''}")
    print(f"Timeout per instance: {args.timeout}s")
    print("-" * 60)

    results = []
    for i, cnf in enumerate(cnf_files, 1):
        print(f"[{i}/{len(cnf_files)}] Running on {cnf.name} ...", end=" ", flush=True)
        res = run_solver(str(solver), str(cnf), args.timeout, args.exact)
        results.append(res)

        sat = res.get("satisfaction_rate")
        t = res.get("wall_time_s")
        status = res.get("status")
        print(f"{status} | sat={sat} | {t}s")

    # Write CSV
    output_path = Path(args.output)
    fieldnames = ["cnf_file", "status", "satisfaction_rate", "satisfied", "unsatisfied",
                  "variables", "clauses", "wall_time_s", "return_code", "error"]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print("-" * 60)
    print(f"Results written to {output_path}")

    # Simple summary
    successful = [r for r in results if r["status"] in ("SATISFIED", "PARTIAL", "UNSATISFIABLE")]
    timeouts = [r for r in results if r["status"] == "TIMEOUT"]
    errors = [r for r in results if r["status"] == "ERROR"]

    avg_sat = None
    if successful:
        sats = [r["satisfaction_rate"] for r in successful if r["satisfaction_rate"] is not None]
        if sats:
            avg_sat = sum(sats) / len(sats)

    print("\n=== SUMMARY ===")
    print(f"Total instances: {len(results)}")
    print(f"Successful:      {len(successful)}")
    print(f"Timeouts:        {len(timeouts)}")
    print(f"Errors:          {len(errors)}")
    if avg_sat is not None:
        print(f"Average satisfaction (successful): {avg_sat:.4f}")

    # Also write a small summary file
    summary_path = output_path.with_suffix(".summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"NitroSAT Benchmark Summary\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Solver: {solver}{' --exact' if args.exact else ''}\n")
        f.write(f"Instances: {len(results)}\n")
        f.write(f"Successful: {len(successful)}\n")
        f.write(f"Timeouts: {len(timeouts)}\n")
        f.write(f"Errors: {len(errors)}\n")
        if avg_sat is not None:
            f.write(f"Avg satisfaction: {avg_sat:.4f}\n")

    print(f"Summary written to {summary_path}")


if __name__ == "__main__":
    main()
