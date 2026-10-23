#!/usr/bin/env python3
"""Check the final handoff chain without rewriting scientific outputs."""

import hashlib
import json
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "docs/协作/主线验收核查.json"
GATES = [f"docs/协作/{name}-ready.json" for name in ("sequence", "methods", "report")]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    if OUTPUT.is_symlink() or OUTPUT.resolve() != OUTPUT or not OUTPUT.resolve().is_relative_to(ROOT):
        raise ValueError("Audit output must be an ordinary file in this repository")
    OUTPUT.write_text(json.dumps({"status": "running", "started_at_utc": datetime.now(timezone.utc).isoformat()}) + "\n")
    snapshots = {}
    original_paths = {}
    checks = {}
    gate_counts = {}
    gates = {}

    def read(relative):
        if Path(relative).is_absolute():
            raise ValueError(f"Expected a repository-relative input: {relative}")
        path = (ROOT / relative).resolve(strict=True)
        if not path.is_relative_to(ROOT):
            raise ValueError(f"Input escapes repository: {relative}")
        data = path.read_bytes()
        key = str(path.relative_to(ROOT))
        current = digest(data)
        if key in snapshots and current != snapshots[key]:
            raise ValueError(f"Input changed during audit: {key}")
        snapshots[key] = current
        if relative in original_paths and original_paths[relative] != path:
            raise ValueError(f"Input link changed during audit: {relative}")
        original_paths[relative] = path
        return data

    def record(name, passed):
        checks[name] = bool(passed)
        if not passed:
            raise ValueError(f"Check failed: {name}")

    own_path = str(Path(__file__).resolve().relative_to(ROOT))
    own_hash = digest(read(own_path))
    for relative in GATES:
        gate = json.loads(read(relative))
        gates[relative] = gate
        record(f"{relative}:ready", gate.get("status") == "ready")
        files = gate.get("files")
        record(f"{relative}:nonempty_bindings", isinstance(files, dict) and bool(files))
        for filename, expected in files.items():
            record(f"sha256:{filename}", digest(read(filename)) == expected)
        inputs = gate.get("inputs", {})
        record(f"{relative}:input_mapping", isinstance(inputs, dict))
        for filename, expected in inputs.items():
            record(f"input_sha256:{filename}", digest(read(filename)) == expected)
        gate_counts[relative] = len(files)

    for child, parents in ((GATES[1], GATES[:1]), (GATES[2], GATES[:2])):
        bound = {**gates[child]["files"], **gates[child].get("inputs", {})}
        for parent in parents:
            record(f"stage_link:{child}:{parent}", bound.get(parent) == digest(read(parent)))
        record(f"design_version:{child}", gates[child].get("design_version") == gates[GATES[0]].get("design_version") == "B-v1")

    primary = json.loads(read("设计/verification_summary.json"))
    independent = json.loads(read("设计/independent_verification.json"))
    record("primary_scientific_checks", primary.get("status") == "passed" and all(primary["checks"].values()))
    record("independent_scientific_checks", independent.get("status") == "pass")
    record("unchanged_ElegantPaper", digest(read("报告/elegantpaper.cls")) == "29370455a0cdf0af94039a6483ba1433f5b067666771b0f96c4cca0159830173")
    record("project_rules_equal", read("AGENTS.md") == read("CLAUDE.md"))

    reachable = set()
    active = set()

    def strip_comments(text):
        lines = []
        for line in text.splitlines():
            for index, char in enumerate(line):
                if char == "%":
                    prefix = line[:index]
                    slashes = len(prefix) - len(prefix.rstrip("\\"))
                    if slashes % 2 == 0:
                        line = line[:index]
                        break
            lines.append(line)
        return "\n".join(lines)

    def visit(relative):
        if relative in active:
            raise ValueError(f"Cyclic TeX input: {relative}")
        active.add(relative)
        reachable.add(relative)
        text = strip_comments(read(relative).decode())
        parts = [text]
        for name in re.findall(r"\\(?:input|include)\s*\{([^}]+)\}", text):
            path = Path("报告") / name
            if not path.suffix:
                path = path.with_suffix(".tex")
            parts.append(visit(str(path)))
        active.remove(relative)
        return "\n".join(parts)

    tex = visit("报告/main.tex")
    record("both_B_sections_reachable", {"报告/sections/methods-cloning.tex", "报告/sections/results-cloning.tex"}.issubset(reachable))
    bib = read("报告/references.bib").decode()
    keys = re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", bib)
    record("unique_bibliography_keys", all(n == 1 for n in Counter(keys).values()))
    cited = set()
    # Handle starred, spaced, optional-note and multiple-key citation forms;
    # reject unsupported syntax instead of silently omitting a citation.
    for match in re.finditer(r"\\([A-Za-z]*cite[A-Za-z]*)\*?\s*", tex, re.I):
        position = match.end()
        multiple = match[1].lower().endswith("cites")
        count = 0
        while True:
            group = re.match(r"\s*(?:\[[^\[\]]*\]\s*)*\{([^{}]+)\}", tex[position:])
            if group is None:
                if count == 0:
                    raise ValueError(f"Unsupported citation syntax: {tex[match.start():match.start()+90]}")
                break
            for key in group[1].split(","):
                key = key.strip()
                if key == "*" and match[1].lower() == "nocite":
                    cited.update(keys)
                elif re.fullmatch(r"[A-Za-z0-9_:./+\-]+", key):
                    cited.add(key)
                else:
                    raise ValueError(f"Unsupported citation key: {key}")
            count += 1
            position += group.end()
            if not multiple:
                break
    record("all_citation_keys_exist", bool(cited) and cited.issubset(set(keys)))
    graphics = re.findall(r"\\includegraphics\s*(?:\[[^\]]*\]\s*)?\{([^}]+)\}", tex)
    graphic_paths = set()
    for name in graphics:
        path = ROOT / "报告" / name
        if not path.suffix:
            matches = [path.with_suffix(ext) for ext in (".pdf", ".png", ".jpg") if path.with_suffix(ext).is_file()]
            if not matches:
                raise ValueError(f"Missing figure: {name}")
            path = matches[0]
        read(str(path.relative_to(ROOT)))
        graphic_paths.add(str(path.relative_to(ROOT)))
    record("data_figures_in_report", {"报告/images/B/plasmid-map.pdf", "报告/images/B/predicted-bands.pdf"}.issubset(graphic_paths))
    record("workflow_figure_in_report", any(Path(p).name.startswith("cloning-workflow") for p in graphic_paths))
    report_bound = {**gates[GATES[2]]["files"], **gates[GATES[2]].get("inputs", {})}
    required_report_files = reachable | graphic_paths | {"报告/main.pdf", "报告/references.bib", "报告/elegantpaper.cls"}
    record("actual_report_inputs_and_PDF_bound", required_report_files.issubset(report_bound))

    for name in ("methods-cloning.tex", "results-cloning.tex"):
        text = read(f"报告/sections/{name}").decode()
        record(f"no_substantive_placeholder:{name}", not re.search(r"\b(?:TBD|TODO|FIXME|placeholder)\b", text, re.I))
    record("working_draft_identity_retained", "Working draft" in tex)

    pdf = "报告/main.pdf"
    read(pdf)
    info = subprocess.run(["pdfinfo", str(ROOT / pdf)], check=True, capture_output=True, text=True).stdout
    match = re.search(r"^Pages:\s+(\d+)", info, re.M)
    record("nonempty_readable_PDF", match is not None and int(match[1]) > 0)
    extracted = subprocess.run(["pdftotext", "-layout", str(ROOT / pdf), "-"], check=True, capture_output=True, text=True).stdout
    record("PDF_has_extractable_content", len(extracted.strip()) > 1000)
    record("PDF_no_undefined_citation_marker", "[?]" not in extracted and "??" not in extracted)
    record("PDF_no_replacement_character", "\ufffd" not in extracted)

    # A second read rejects concurrent edits rather than binding a mixed snapshot.
    for original, resolved in original_paths.items():
        if (ROOT / original).resolve(strict=True) != resolved:
            raise ValueError(f"Input link changed before publication: {original}")
    for filename, expected in list(snapshots.items()):
        if digest((ROOT / filename).read_bytes()) != expected:
            raise ValueError(f"Input changed before audit publication: {filename}")
    payload = {
        "status": "passed",
        "scope": "Read-only final structural/hash/reference audit; manual science and page review are recorded separately in docs/B验收记录.md.",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "gate_file_counts": gate_counts,
        "checks": checks,
        "citation_keys_used": sorted(cited),
        "bibliography_entries": len(keys),
        "graphics_includes": graphics,
        "reachable_TeX_files": sorted(reachable),
        "pdf_pages": int(match[1]),
        "files": snapshots,
        "audit_script_sha256": own_hash,
        "limitations": ["No wet-lab data", "Does not validate independent A protocols", "Text extraction does not replace visual PDF inspection", "Does not assert Git/PR delivery"],
    }
    if OUTPUT.is_symlink() or OUTPUT.resolve() != OUTPUT:
        raise ValueError("Audit output path changed before publication")
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": payload["status"], "checks": len(checks), "bound_files": len(snapshots), "pdf_pages": payload["pdf_pages"]}))


if __name__ == "__main__":
    main()
