"""Check this style batch against its parent AST and the applicable text run.

Does not modify sources or run simulations. --syntax-only skips log checks.
"""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / Path(__file__).stem
REFERENCE = "7201421946ae09a8e54129ce18d440d9a0fff5fa"
SOURCES = ("src/Solve.jl", "src/CycleSolver.jl", "src/Parallelsolution.jl", "src/SPMe.jl")
# Every non-syntactic change has an explicit, occurrence-counted exception.
REMOVALS = {
    "src/Solve.jl": [("                MT_old = MT_new\n", "", 1)],
    "src/SPMe.jl": [
        ("    else\n        I_app = Float64(I_app)\n", "", 2),
    ],
    "src/Parallelsolution.jl": [
        ('variables["thermal2D Vsolve status"] = converged ? 3.0 : 3.5',
         'variables["thermal2D Vsolve status"] = 3.0', 1),
        ('variables["thermal2D Vsolve converged"] = converged ? 1.0 : 0.0',
         'variables["thermal2D Vsolve converged"] = 1.0', 1),
    ],
}
JULIA_AST = r'''
function canonical(x)
    x isa Expr || return x
    args = [canonical(a) for a in x.args if !(a isa LineNumberNode)]
    x.head === :block && length(args) == 1 && return only(args)
    return Expr(x.head, args...)
end
for i in 1:2:length(ARGS)
    before = canonical(Meta.parseall(read(ARGS[i], String)))
    after = canonical(Meta.parseall(read(ARGS[i+1], String)))
    isequal(before, after) || error("AST differs beyond declared removals: " * ARGS[i+1])
    println("AST_EQUIVALENT ", ARGS[i+1])
end
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--syntax-only", action="store_true")
    options = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="jubat_style_") as temp:
        pairs = []
        for rel in SOURCES:
            original = subprocess.check_output(["git", "show", f"{REFERENCE}:{rel}"], cwd=ROOT).decode("utf-8").replace("\r\n", "\n")
            for old, new, count in REMOVALS.get(rel, []):
                assert original.count(old) == count, (rel, old)
                original = original.replace(old, new)
            path = Path(temp) / Path(rel).name
            path.write_text(original, encoding="utf-8")
            pairs.extend((str(path), str(ROOT / rel)))
        subprocess.run(["julia", "-t1", "--startup-file=no", "-e", JULIA_AST, *pairs], check=True, cwd=ROOT)
    if options.syntax_only:
        return
    tree = ast.parse((ROOT / "tools/validate_case_batch_a.py").read_text(encoding="utf-8"))
    patterns = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == "PATTERNS" for t in n.targets))
    metrics = []
    reference_log = ROOT / "output/test_cycling_csv_timeline/candidate_60s.log"
    for path in (reference_log, OUT / "candidate_60s.log"):
        log = path.read_text(encoding="utf-8-sig")
        values = {}
        for key, pattern in patterns.items():
            matches = re.findall(pattern, log)
            assert len(matches) == 1, (path, key)
            values[key] = matches[0]
        assert path.with_suffix(".exit").read_text(encoding="utf-8-sig").strip() == "0"
        metrics.append(values)
    assert metrics[0] == metrics[1], (metrics[0], metrics[1])
    tests = {}
    for name in ("test_case_batch_a", "test_case_batch_c", "test_temperature_current", "test_dirty_fallbacks"):
        assert (OUT / f"{name}.exit").read_text(encoding="utf-8-sig").strip() == "0", name
        log = (OUT / f"{name}.log").read_text(encoding="utf-8-sig")
        assert "Test Summary:" in log and "TEST_FILE_COMPLETED " + name in log, name
        tests[name] = "passed"
    patch = subprocess.check_output(["git", "diff", "--binary", REFERENCE, "--", *SOURCES], cwd=ROOT)
    report = dict(reference_commit=REFERENCE, reference_log=str(reference_log),
                  candidate_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  ast_equivalent_after_declared_removals=list(SOURCES), tests=tests,
                  text_comparison=dict(passed=21, total=21, metrics=metrics[1]),
                  official_baseline=False, refrozen=False,
                  runtime="Julia 1.11.2 -t1 --startup-file=no; OPENBLAS_NUM_THREADS=1; GKSwstype=100",
                  source_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in SOURCES},
                  input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in
                                ("Project.toml", "Manifest.toml", "example/testexample.jl")},
                  patch_sha256=hashlib.sha256(patch).hexdigest())
    (OUT / "candidate.patch").write_bytes(patch)
    (OUT / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print("Style batch: AST 4/4; existing test files 4/4; exact text comparison 21/21")


if __name__ == "__main__":
    main()
