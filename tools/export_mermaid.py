import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


MERMAID_BLOCK_RE = re.compile(r"```mermaid\s*\n(.*?)\n```", re.DOTALL | re.IGNORECASE)


def detect_name(index: int, content: str) -> str:
    head = content.strip().splitlines()
    first = next((ln.strip() for ln in head if ln.strip()), "").lower()
    base = f"figure_mermaid_{index+1}"
    if first.startswith("sequence"):
        return "figure_coupling_sequence"
    if first.startswith("flowchart"):
        # 粗略识别：含 Bus/Parallel 词汇则判作并联系统图
        low = content.lower()
        if "parallel" in low or "bus" in low or "branches" in low:
            return "figure_parallel_network"
        # 其余 flowchart 默认为耦合数据流示意
        return "figure_model_coupling_from_md"
    return base


def normalize_mermaid(content: str) -> str:
    """Normalize mermaid code to improve CLI parsing robustness.
    - Replace literal \n inside labels with <br/>
    - Strip BOM and trailing spaces on lines
    """
    # Replace literal backslash-n with mermaid-supported line break
    content = content.replace("\\n", "<br/>")
    # Normalize some Unicode that may break mermaid parser
    replacements = {
        "→": "->",
        "←": "<-",
        "↔": "<->",
        "（": "(",
        "）": ")",
        "Σ": "Sigma",
        "·": "*",
    }
    for k, v in replacements.items():
        content = content.replace(k, v)
    # Strip possible BOM
    if content and content[0] == "\ufeff":
        content = content[1:]
    # Trim trailing spaces per line
    content = "\n".join(ln.rstrip() for ln in content.splitlines()) + "\n"
    return content


def _sanitize_flowchart_basic(content: str) -> str:
    # Switch to graph keyword
    content = re.sub(r"^\s*flowchart\s+", "graph ", content, flags=re.IGNORECASE | re.MULTILINE)
    # Simplify special shape syntaxes: id((label)) -> id[label], id([label]) -> id[label], id[(label)] -> id[label]
    content = re.sub(r"\(\(([^)]+)\)\)", r"[\1]", content)
    content = re.sub(r"\(\[([^\]]+)\]\)", r"[\1]", content)
    content = re.sub(r"\[\(([^)]+)\)\]", r"[\1]", content)
    # Inside node labels, strip parentheses to avoid parser confusion
    def _strip_parens_in_labels(m: re.Match) -> str:
        label = m.group(1)
        label = label.replace("(", " ").replace(")", " ")
        label = label.replace("/", " / ")
        label = label.replace(",", ", ")
        label = re.sub(r"\s+", " ", label).strip()
        return f"[{label}]"
    content = re.sub(r"\[([^\]]+)\]", _strip_parens_in_labels, content)
    return content


def _resolve_npx_cmd() -> list[str]:
    """Return base npx invocation adapted for Windows (npx.cmd) or others (npx)."""
    candidates = [["npx"], ["npx.cmd"], ["cmd", "/c", "npx"]]
    for cand in candidates:
        try:
            subprocess.run(cand + ["-v"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return cand
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue
    raise FileNotFoundError("npx executable not found. Ensure Node.js is installed and npx is in PATH.")


def run_mmdc(mmd_path: Path, out_svg: Path, out_png: Path, config: Path | None, scale: float = 1.5):
    npx_cmd = _resolve_npx_cmd()
    cmd_base = npx_cmd + [
        "-y", "@mermaid-js/mermaid-cli",
        "-i", str(mmd_path),
        "--backgroundColor", "white",
    ]
    if config is not None:
        cmd_base += ["-C", str(config)]

    # SVG
    cmd_svg = cmd_base + ["-o", str(out_svg)]
    subprocess.run(cmd_svg, check=True)

    # PNG
    cmd_png = cmd_base + ["-o", str(out_png), "--scale", str(scale)]
    subprocess.run(cmd_png, check=True)


def main():
    ap = argparse.ArgumentParser(description="Export Mermaid code blocks in a Markdown file to PNG/SVG using mermaid-cli (mmdc)")
    ap.add_argument("md_file", help="Path to Markdown file containing ```mermaid blocks")
    ap.add_argument("--out", default=".", help="Output directory (default: current dir)")
    ap.add_argument("--config", default=None, help="Mermaid config JSON file (optional)")
    ap.add_argument("--only", default=None, help="Optional comma-separated indices (1-based) of mermaid blocks to export, e.g. 2,3")
    args = ap.parse_args()

    md_path = Path(args.md_file).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    config_path = Path(args.config).resolve() if args.config else None

    text = md_path.read_text(encoding="utf-8")
    blocks = MERMAID_BLOCK_RE.findall(text)
    if not blocks:
        print("No mermaid blocks found.")
        return 0

    idx_filter = None
    if args.only:
        try:
            idx_filter = {int(x.strip()) - 1 for x in args.only.split(',') if x.strip()}
        except Exception:
            print("Invalid --only format. Expected comma-separated integers, e.g. 2,3", file=sys.stderr)
            return 2

    exported = []
    for i, content in enumerate(blocks):
        if idx_filter is not None and i not in idx_filter:
            continue
        name = detect_name(i, content)
        out_svg = out_dir / f"{name}.svg"
        out_png = out_dir / f"{name}.png"
        # 写入临时 .mmd 文件
        with tempfile.TemporaryDirectory() as td:
            mmd_file = Path(td) / f"{name}.mmd"
            mmd_file.write_text(normalize_mermaid(content), encoding="utf-8")
            try:
                run_mmdc(mmd_file, out_svg, out_png, config_path)
                exported.append((name, out_svg, out_png))
            except FileNotFoundError as e:
                print("Error: npx not found. Please install Node.js and ensure 'npx' is available in PATH.", file=sys.stderr)
                return 127
            except subprocess.CalledProcessError as e:
                # Fallback: try to sanitize flowchart syntax and retry once
                raw = mmd_file.read_text(encoding="utf-8")
                if raw.lstrip().lower().startswith("flowchart") or raw.lstrip().lower().startswith("graph"):
                    sanitized = _sanitize_flowchart_basic(raw)
                    mmd_file.write_text(sanitized, encoding="utf-8")
                    try:
                        run_mmdc(mmd_file, out_svg, out_png, config_path)
                        exported.append((name + "_sanitized", out_svg, out_png))
                        continue
                    except subprocess.CalledProcessError as e2:
                        print(f"mmdc failed for block {i+1} ({name}) even after sanitize: {e2}", file=sys.stderr)
                        return e2.returncode
                print(f"mmdc failed for block {i+1} ({name}): {e}", file=sys.stderr)
                return e.returncode

    for name, svg, png in exported:
        print(f"Saved: {svg}")
        print(f"Saved: {png}")
    # 建议与 figure_model_coupling.* 保持同款尺寸/背景，可在 config 中调整主题变量
    return 0


if __name__ == "__main__":
    sys.exit(main())
