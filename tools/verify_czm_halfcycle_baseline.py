"""Read-only exact comparison of saved Case2 discharge/rest half-cycle evidence.

No solver is called and no files are written. Each optional source root validates
that run's own manifest; changed source hashes between runs are reported, allowed.
"""
import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import sys
import tomllib

PHASES = ("discharge", "rest1")
SUMMARY = "phase,t_start_s,t_end_s,duration_s,V_start_V,V_end_V,capacity_Ah,termination_reason,T_max_K,T_mean_end_K,D_max_element,D_mean_element,D_max_gp,Dv_max_gp,n_fractured_gp,hoop_min_Pa,hoop_max_Pa,shear_min_Pa,shear_max_Pa,recovered_element_vm_max_history_Pa,n_recorded_times".split(",")
GP = "cohesive,gp,interface,D,D_visc,fractured,delta_max_n,delta_max_t,delta_max_eff,accumulated_damage".split(",")
STATE = {"t_global_s", "mechanics_initialized", "gp_history_present", "temperature_K_extrema", "element_D_extrema", "element_Dv_extrema", "element_fractured"}
IGNORED_METRICS = {"run_directory", "source_head"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value, where):
    try:
        result = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"{where}: invalid number {value!r}") from exc
    require(result.is_finite(), f"{where}: nonfinite number {value!r}")
    return result


def read_csv(path, columns, keys, strings=()):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames == columns, f"{path}: missing, extra, duplicate or reordered columns")
        result = {}
        for index, row in enumerate(reader, 2):
            require(None not in row and all(v is not None and v != "" for v in row.values()), f"{path}:{index}: incomplete row")
            parsed = {k: v if k in strings else number(v, f"{path}:{index}:{k}") for k, v in row.items()}
            key = tuple(parsed[k] for k in keys)
            require(key not in result, f"{path}:{index}: duplicate key {key}")
            result[key] = parsed
        require(result, f"{path}: no rows")
        return result


def compare(left, right, where):
    require(type(left) is type(right), f"{where}: value types differ")
    if isinstance(left, dict):
        require(left.keys() == right.keys(), f"{where}: keys differ")
        for key in left:
            compare(left[key], right[key], f"{where}/{key}")
    elif isinstance(left, (tuple, list)):
        require(len(left) == len(right), f"{where}: lengths differ")
        for index, (a, b) in enumerate(zip(left, right)):
            compare(a, b, f"{where}/{index}")
    else:
        require(left == right, f"{where}: {left!r} != {right!r}")


def finite_tree(value, where):
    if isinstance(value, dict):
        return {k: finite_tree(v, f"{where}/{k}") for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_tree(v, where) for v in value]
    if isinstance(value, (int, float, Decimal)) and not isinstance(value, bool):
        return number(value, where)
    return value


def read_state(path):
    result = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        key, value = line.split("=", 1)
        require(key not in result, f"{path}: duplicate state key {key}")
        if value in ("true", "false"):
            parsed = value == "true"
        elif value.startswith("("):
            require(value.endswith(")"), f"{path}: invalid extrema")
            parts = value[1:-1].split(",")
            require(len(parts) == 2, f"{path}: invalid extrema")
            parsed = [number(v.strip(), f"{path}/{key}") for v in parts]
        else:
            parsed = number(value, f"{path}/{key}")
        result[key] = finite_tree(parsed, f"{path}/{key}")
    require(result.keys() == STATE, f"{path}: state schema differs")
    return result


def read_metrics(path):
    with path.open("rb") as stream:
        record = tomllib.load(stream, parse_float=Decimal)
    for key in ("initial", "mesh", "phases", "interfaces", "exit_code", "first_recorded_nonzero_element_damage_snapshot"):
        require(key in record, f"{path}: missing metrics {key}")
    require(set(record["phases"]) == set(PHASES), f"{path}: phase schema differs")
    require(set(record["mesh"]) == {"n_theta", "cohesive_elements", "mechanical_bulk_elements", "mechanical_nodes", "thermal_elements", "thermal_nodes"}, f"{path}: mesh schema differs")
    require(record["exit_code"] == 0, f"{path}: metrics exit is not zero")
    return finite_tree({k: v for k, v in record.items() if k not in IGNORED_METRICS}, str(path))


def manifest(run, source_root):
    path = run / "source_manifest.csv"
    rows = read_csv(path, ["relative_path", "sha256"], ["relative_path"], ["relative_path", "sha256"])
    result = {}
    for (relative,), row in rows.items():
        p = Path(relative)
        require(not p.is_absolute() and ".." not in p.parts, f"{path}: unsafe relative source path")
        digest = row["sha256"].lower()
        require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest), f"{path}: invalid SHA256")
        result[relative] = digest
        if source_root is not None:
            actual = hashlib.sha256((source_root / p).read_bytes()).hexdigest()
            require(actual == digest, f"{path}: source hash differs for {relative}")
    return result


def read_run(run, source_root):
    require((run / "exit_code.txt").read_text(encoding="utf-8-sig").strip() == "0", f"{run}: exit code is not zero")
    require((run / "completion.txt").read_text(encoding="utf-8-sig").strip() == "SUCCESS: discharge and rest1 only; no charge/rest2 run.", f"{run}: missing successful completion")
    summary = read_csv(run / "phase_summary.csv", SUMMARY, ["phase"], ["phase", "termination_reason"])
    require(list(summary) == [(p,) for p in PHASES], f"{run}: missing, extra or reordered phases")
    metrics = read_metrics(run / "metrics.toml")
    gps = {}
    states = {name: read_state(run / f"{name}_state.txt") for name in ("initial", *PHASES)}
    for phase in PHASES:
        table = read_csv(run / f"{phase}_gp.csv", GP, ["cohesive", "gp"], ["interface", "fractured"])
        require(len(table) == 2 * metrics["mesh"]["cohesive_elements"], f"{run}/{phase}: GP count differs from mesh")
        expected = {(Decimal(e), Decimal(g)) for e in range(1, int(metrics["mesh"]["cohesive_elements"]) + 1) for g in (1, 2)}
        require(table.keys() == expected, f"{run}/{phase}: GP key coverage differs")
        require(all(v["fractured"] in ("true", "false") and v["interface"] in ("PE_PCC", "NE_NCC") for v in table.values()), f"{run}/{phase}: invalid GP label")
        gps[phase] = table
        values = summary[(phase,)]
        for key, value in values.items():
            if key != "phase":
                compare(value, metrics["phases"][phase][key], f"{run}/{phase}/metrics-summary/{key}")
        step = metrics["phases"][phase]
        require(step["n_time_intervals_excluding_initial"] == step["n_recorded_times"] - 1, f"{run}/{phase}: inconsistent step count")
        require(step["n_recorded_mechanical_snapshots"] == step["n_time_intervals_excluding_initial"], f"{run}/{phase}: inconsistent mechanical count")
    return dict(summary=summary, gps=gps, states=states, metrics=metrics), manifest(run, source_root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--reference-source-root", type=Path)
    parser.add_argument("--candidate-source-root", type=Path)
    args = parser.parse_args()
    left, left_sources = read_run(args.reference, args.reference_source_root)
    right, right_sources = read_run(args.candidate, args.candidate_source_root)
    compare(left, right, "halfcycle")
    changed = [p for p in sorted(left_sources.keys() | right_sources.keys()) if left_sources.get(p) != right_sources.get(p)]
    print(json.dumps(dict(passed=True, phases=list(PHASES), summary_fields=len(SUMMARY), gp_rows={p: len(left["gps"][p]) for p in PHASES}, gp_fields=len(GP), source_changes=changed, source_roots_verified=dict(reference=args.reference_source_root is not None, candidate=args.candidate_source_root is not None), ignored_metric_identity=sorted(IGNORED_METRICS)), ensure_ascii=False))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, SyntaxError, TypeError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
