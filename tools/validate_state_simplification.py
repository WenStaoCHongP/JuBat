"""Compare this batch's same-source 60 s text runs; never freeze a baseline."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path)
parser.add_argument("--reference-evidence")
parser.add_argument("--reference-worktree")
parser.add_argument("--prepare", action="store_true")
args = parser.parse_args()
if not args.prepare and (not args.reference_evidence or not args.reference_worktree):
    parser.error("comparison requires --reference-evidence and --reference-worktree")
OUT = args.output.resolve() if args.output else ROOT / "output" / Path(__file__).stem
OUT.mkdir(parents=True, exist_ok=True)
PATTERNS = {
    "theta_elements": r"周向单元数 n_theta: (\d+)",
    "thermal_elements": r"总单元数 ne: (\d+)",
    "thermal_nodes": r"总节点数 nT: (\d+)",
    "time_steps": r"总时间步数: (\d+)",
    "initial_voltage_v": r"初始电压: (\S+) V",
    "final_voltage_v": r"最终电压: (\S+) V",
    "voltage_drop_v": r"电压降: (\S+) V",
    "final_capacity_ah": r"最终容量: (\S+) Ah",
    "temperature_min_k": r"温度范围: (\S+) K ~ \S+ K",
    "temperature_max_k": r"温度范围: \S+ K ~ (\S+) K",
    "czm_d_max_percent": r"CZM 最终 D_max: (\S+)%",
    "czm_d_mean_percent": r"CZM 最终 D_mean: (\S+)%",
    "czm_max_effective_separation_m": r"CZM 最大有效分离位移: (\S+) m",
    "czm_fractured_elements": r"CZM 断裂单元数: (\d+)",
    "hoop_stress_min_mpa": r"最终环向应力范围: (\S+) ~ \S+ MPa",
    "hoop_stress_max_mpa": r"最终环向应力范围: \S+ ~ (\S+) MPa",
    "shear_stress_min_mpa": r"最终切向剪应力范围: (\S+) ~ \S+ MPa",
    "shear_stress_max_mpa": r"最终切向剪应力范围: \S+ ~ (\S+) MPa",
    "von_mises_max_mpa": r"全时域最大高斯点 von Mises 应力: (\S+) MPa",
    "equivalent_plastic_strain_max": r"全时域最大等效塑性应变: (\S+)",
}


def metrics(label):
    text = (OUT / f"{label}_60s.log").read_text(encoding="utf-8-sig")
    result = {}
    for name, pattern in PATTERNS.items():
        matches = re.findall(pattern, text)
        assert len(matches) == 1, (label, name, matches)
        result[name] = matches[0]
    result["exit_code"] = (OUT / f"{label}_60s.exit").read_text(encoding="utf-8-sig").strip()
    return result


if args.prepare:
    run={"source_sha256":{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/"src").glob("*.jl"))},
         "example_sha256":hashlib.sha256((ROOT/"example/testexample.jl").read_bytes()).hexdigest(),
         "dependency_sha256":{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ("Project.toml","Manifest.toml")},
         "head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()}
    (OUT/"candidate_run_identity.json").write_text(json.dumps(run,indent=2)+"\n",encoding="utf-8")
    print("Plain-state candidate source/environment identity captured")
    sys.exit(0)
reference, candidate = metrics("reference"), metrics("candidate")
checks = [dict(name=k, reference=v, candidate=candidate[k], passed=v == candidate[k])
          for k, v in reference.items()]
assert len(checks) == 21
diff = subprocess.check_output(["git", "diff", "--binary"], cwd=ROOT, stderr=subprocess.DEVNULL)
manifest = {
    "reference_head": "fd7f4f0ad35078296b5631d8a6b553d324417192",
    "reference_worktree": args.reference_worktree,
    "reference_evidence": args.reference_evidence,
    "candidate_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "candidate_tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
    "source_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in sorted((ROOT / "src").glob("*.jl"))},
    "example_sha256": hashlib.sha256((ROOT / "example/testexample.jl").read_bytes()).hexdigest(),
    "dependency_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ("Project.toml", "Manifest.toml") if (ROOT / name).exists()},
    "julia": "1.11.2", "threads": 1, "startup_file": False,
    "OPENBLAS_NUM_THREADS": "1", "GKSwstype": "100", "n_theta": 20, "duration_s": 60,
}
run=json.loads((OUT/"candidate_run_identity.json").read_text(encoding="utf-8"))
for key in ("source_sha256","example_sha256","dependency_sha256"):
    assert manifest[key]==run[key],("candidate changed during run",key)
assert manifest["candidate_head"]==run["head"]
manifest["reference_log_sha256"]=hashlib.sha256((OUT/"reference_60s.log").read_bytes()).hexdigest()
report = dict(official_baseline=False, refrozen=False, figures_generated=False,
              checks=checks, passed=sum(c["passed"] for c in checks), total=len(checks), identity=manifest)
(OUT / "comparison_60s.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
(OUT / "candidate_tracked.patch").write_bytes(diff)
print(f"Exact text comparison: {report['passed']}/{report['total']}; official_baseline=false; refrozen=false")
assert all(c["passed"] for c in checks) and candidate["exit_code"] == "0"
