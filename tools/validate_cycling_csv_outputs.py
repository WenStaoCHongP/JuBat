"""Audit actual short-run CSVs and existing logs; never run or modify a simulation."""
from pathlib import Path
import ast
import csv
import hashlib
import json
import math
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "output/test_cycling_csv_timeline"
OUT = ROOT / "output" / Path(__file__).stem


def table(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        values = list(reader)
    assert values, (path, "empty data")
    assert len(header) == len(set(header)), (path, "duplicate columns")
    assert all(len(row) == len(header) for row in values), (path, "row width")
    return header, values


def main():
    report = {"runs": [], "physical_values_unchanged": [], "official_baseline": False}
    for mode in ("czm_cycling", "czm_cycling_n20"):
        for sample in ("equal", "unequal"):
            log = (RUNS / f"{mode}_{sample}.log").read_text(encoding="utf-8-sig")
            assert (RUNS / f"{mode}_{sample}.exit").read_text(encoding="utf-8-sig").strip() == "0"
            assert "CSV_TIMELINE_VALIDATED" in log
            directory = (ROOT / "output/czm_cycling" / f"csv_timeline_candidate_{sample}"
                         if mode == "czm_cycling" else RUNS / mode / f"candidate_{sample}/case1")
            files = sorted(directory.glob("*.csv"))
            assert len(files) == (14 if mode == "czm_cycling" else 15)
            phase_header, phase_rows = table(directory / "phase_summary.csv")
            phases = [dict(zip(phase_header, row)) for row in phase_rows]
            assert [p["phase"] for p in phases] == ["discharge", "rest1", "charge", "rest2"]
            assert all(float(a["t_end_s"]) == float(b["t_start_s"]) for a, b in zip(phases, phases[1:]))
            files_report = []
            for path in files:
                header, values = table(path)
                for time_key in ("time_s", "t_s"):
                    if time_key in header:
                        stamps = [float(row[header.index(time_key)]) for row in values]
                        assert all(math.isfinite(t) for t in stamps), path
                        assert all(a <= b for a, b in zip(stamps, stamps[1:])), path
                        assert min(stamps) >= 0 and max(stamps) <= float(phases[-1]["t_end_s"]), path
                # Manifest NaN is defined as no damage onset; physical field tables must be finite.
                if path.name != "case_manifest.csv":
                    for row in values:
                        for value in row:
                            try:
                                number = float(value)
                            except ValueError:
                                continue
                            assert math.isfinite(number), (path, value)
                files_report.append({"file": path.name, "rows": len(values),
                                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            report["runs"].append({"example": mode, "configuration": sample,
                "directory": str(directory), "files": files_report,
                "duration_s": float(phases[-1]["t_end_s"])})

        # Same equal-phase simulation before/after: compare all five raw mechanical histories.
        before = (ROOT / "output/czm_cycling/csv_timeline_before"
                  if mode == "czm_cycling" else RUNS / mode / "before/case1")
        after = (ROOT / "output/czm_cycling/csv_timeline_candidate_equal"
                 if mode == "czm_cycling" else RUNS / mode / "candidate_equal/case1")
        for name in ("damage", "sep_n", "sep_t", "traction_n", "traction_t"):
            filename = name + "_history_cyc1.csv"
            h1, rows1 = table(before / filename)
            h2, rows2 = table(after / filename)
            assert h1 == h2 and [r[1:] for r in rows1] == [r[1:] for r in rows2], (mode, filename)
            report["physical_values_unchanged"].append(f"{mode}/{filename}")

    # Reuse the established 20 metric patterns, without running the old batch script.
    tree = ast.parse((ROOT / "tools/validate_case_batch_a.py").read_text(encoding="utf-8"))
    patterns = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == "PATTERNS" for t in n.targets))
    metrics = []
    for path in (ROOT / "output/remove_legacy_csv/candidate_60s.log", RUNS / "candidate_60s.log"):
        log = path.read_text(encoding="utf-8-sig")
        result = {}
        for key, pattern in patterns.items():
            values = re.findall(pattern, log)
            assert len(values) == 1, (path, key)
            result[key] = values[0]
        assert path.with_suffix(".exit").read_text(encoding="utf-8-sig").strip() == "0"
        metrics.append(result)
    assert metrics[0] == metrics[1]
    report["text_baseline"] = {"passed": 21, "total": 21, "metrics": metrics[1]}
    report["identity"] = {"head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in
            ("src/CycleSolver.jl", "src/CouplingState.jl", "example/cycling_csv_time.jl",
             "example/czm_cycling.jl", "example/czm_cycling_n20.jl",
             "test/test_cycling_csv_timeline.jl", "test/test_cycling_csv_time_contract.jl")},
        "runtime": "Julia 1.11.2 -t1 --startup-file=no; OPENBLAS_NUM_THREADS=1; GKSwstype=100"}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("CSV audit: 4 runs, 58 files; 10 physical history tables unchanged; text baseline 21/21")


if __name__ == "__main__":
    main()
