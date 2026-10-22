"""Read-only checks for recorded source identity, scoped sync, and task documents.

This command prints JSON. Copying, report writing, deletion, and commits are separate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def identity(manifest):
    record = read_json(manifest)
    actual_sources = {p.relative_to(ROOT).as_posix(): sha(p) for p in (ROOT / "src").glob("*.jl")}
    assert actual_sources == record["source_sha256"], "source identity differs"
    assert sha(ROOT / "example/testexample.jl") == record["example_sha256"], "example identity differs"
    assert all(sha(ROOT / p) == value for p, value in record["dependency_sha256"].items()), "dependency identity differs"
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    assert head == record["head"], "HEAD differs"
    assert not subprocess.check_output(["git", "diff", "--cached", "--name-only"], cwd=ROOT), "index is not empty"
    return {"check": "identity", "head": head, "index_empty": True, "sources": len(actual_sources), "passed": True}


def sync(receipt, kb):
    record = read_json(receipt)
    checked = []
    for item in record["paths"]:
        local, remote = ROOT / item["path"], kb / item["kb_path"]
        expected = item["copied_sha256"]
        assert sha(local) == expected, f"local sync hash differs: {local}"
        assert sha(remote) == expected, f"KB sync hash differs: {remote}"
        checked.append({"path": item["path"], "sha256": expected})
    return {"check": "sync", "commit": record["commit"], "files": checked, "passed": True}


def documents(source, target, names, receipt):
    prior = None if receipt is None else read_json(receipt)
    if source is None:
        assert prior is not None, "docs requires --source or a verified archive --receipt"
        source = Path(prior["source_directory"])
    prior_rows = {} if prior is None else {x["name"]: x for x in prior["files"]}
    if source.exists():
        assert {p.name for p in source.iterdir()} == set(names), "source file set differs"
    else:
        assert prior is not None, "missing source requires its verified archive receipt"
    rows = []
    for name in names:
        a, b = source / name, target / name
        expected = prior_rows[name]["source_sha256"] if prior is not None else sha(a)
        if a.exists():
            assert sha(a) == expected, f"source differs from recorded archive identity: {a}"
        assert sha(b) == expected, f"archive hash differs: {b}"
        text = b.read_bytes().decode("utf-8", errors="strict")
        fence = None
        for line in text.splitlines():
            match = re.match(r"^\s*(`{3,}|~{3,})", line)
            if match:
                token = match.group(1)
                if fence is None:
                    fence = token
                elif token[0] == fence[0] and len(token) >= len(fence):
                    fence = None
        assert fence is None, f"unclosed code fence: {b}"
        links = []
        for raw in re.findall(r"(?<!!)\[[^\]]*\]\(([^)]+)\)", text):
            link = raw.strip().strip("<>")
            if link.startswith(("http://", "https://", "app://", "plugin://", "codex://", "#", "mailto:")):
                continue
            link = re.sub(r":\d+$", "", link.split("#", 1)[0])
            if not link:
                continue
            path = Path(link)
            resolved = path if path.is_absolute() else b.parent / path
            assert resolved.exists(), f"missing local document link: {b}: {link}"
            links.append({"link": link, "target": str(resolved.resolve())})
        rows.append({"name": name, "source": str(a), "target": str(b), "source_sha256": expected,
                     "target_sha256": sha(b), "utf8_valid": True, "fences_balanced": True, "links": links})
    return {"check": "documents", "source_directory": str(source), "target_directory": str(target),
            "source_exists": source.exists(), "files": rows, "passed": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    p = commands.add_parser("identity")
    p.add_argument("--manifest", required=True, type=Path)
    p = commands.add_parser("sync")
    p.add_argument("--receipt", required=True, type=Path)
    p.add_argument("--kb", required=True, type=Path)
    p = commands.add_parser("docs")
    p.add_argument("--source", type=Path)
    p.add_argument("--target", required=True, type=Path)
    p.add_argument("--receipt", type=Path)
    p.add_argument("--files", nargs="+", default=["task_plan.md", "findings.md", "progress.md", "completion_report.md"])
    args = parser.parse_args()
    if args.mode == "identity":
        result = identity(args.manifest)
    elif args.mode == "sync":
        result = sync(args.receipt, args.kb)
    else:
        result = documents(args.source, args.target, args.files, args.receipt)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
