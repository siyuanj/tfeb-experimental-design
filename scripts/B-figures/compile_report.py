#!/usr/bin/env python3
"""Compile a stable report snapshot and bind its generated review artifacts.

This script performs no scientific calculation or visual inspection. Its passed
record only binds a completed build to its source files, stage gates and output
artifacts. Subsequent consistency and per-page visual reviews remain required.
"""

import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from uuid import uuid4
from zoneinfo import ZoneInfo

from handoff import ROOT, assert_unchanged, sha256, verify_gate


REPORT = ROOT / "报告"
BUILD = REPORT / "build"
RECORD = ROOT / "docs/协作/报告编译核验.json"
SCRIPT = Path(__file__).resolve()
IMAGE_SUFFIXES = {".pdf", ".svg", ".png", ".json"}


def now():
    """Return the user-facing local timestamp, including its UTC offset."""
    return datetime.now(ZoneInfo("America/New_York")).isoformat()


def relative(path):
    """Reject paths that resolve outside the project before recording them."""
    if not path.resolve(strict=True).is_relative_to(ROOT):
        raise ValueError(f"Path escapes this repository: {path}")
    return str(path.relative_to(ROOT))


def write_record(payload):
    """Atomically replace the record so an interrupted write cannot pass."""
    RECORD.parent.mkdir(parents=True, exist_ok=True)
    temporary = RECORD.with_name(f".{RECORD.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        temporary.replace(RECORD)
    finally:
        temporary.unlink(missing_ok=True)


def source_snapshot():
    """Hash the complete requested source set, excluding generated main.pdf."""
    paths = [
        REPORT / "main.tex",
        REPORT / "references.bib",
        REPORT / "elegantpaper.cls",
    ]
    sections = sorted((REPORT / "sections").glob("*.tex"))
    if not sections:
        raise ValueError("No report section sources were found")
    paths.extend(sections)
    image_directory = REPORT / "images/B"
    if not image_directory.is_dir():
        raise FileNotFoundError(image_directory)
    paths.extend(
        path for path in image_directory.rglob("*")
        if path.is_file() and path.suffix.lower() in IMAGE_SUFFIXES
    )
    result = {}
    for path in sorted(set(paths)):
        if not path.is_file():
            raise FileNotFoundError(path)
        result[relative(path)] = sha256(path)
    return result


def artifact(path):
    """Require a real, nonempty output before binding its bytes."""
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty build artifact: {path}")
    return {"path": relative(path), "sha256": sha256(path)}


def run_command(argv, commands, output_path=None):
    """Run an argument vector without a shell and record its actual outcome."""
    entry = {
        "argv": [str(item) for item in argv],
        "command": shlex.join([str(item) for item in argv]),
        "cwd": str(REPORT),
        "started_at": now(),
    }
    commands.append(entry)
    if output_path is None:
        completed = subprocess.run(
            entry["argv"], cwd=REPORT, shell=False, check=False,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
        )
        output = completed.stdout
    else:
        with output_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(
                entry["argv"], cwd=REPORT, shell=False, check=False,
                stdout=handle, stderr=subprocess.STDOUT,
                text=True, encoding="utf-8", errors="replace",
            )
        output = output_path.read_text(encoding="utf-8", errors="replace")
        entry["output_path"] = relative(output_path)
    entry["returncode"] = completed.returncode
    entry["completed_at"] = now()
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command exited with {completed.returncode}: {entry['command']}"
        )
    return output


def main():
    commands = []
    payload = {
        "status": "running",
        "started_at": now(),
        "scope": "Compilation and artifact binding only; visual review is separate",
        "actual_command": shlex.join([sys.executable, str(SCRIPT)]),
        "commands": commands,
    }
    # Invalidate a previous passed record before checking any prerequisite.
    write_record(payload)
    try:
        _, sequence = verify_gate("docs/协作/sequence-ready.json")
        _, methods = verify_gate("docs/协作/methods-ready.json")
        inputs = source_snapshot()
        script_hash = sha256(SCRIPT)
        payload.update(
            inputs=inputs, sequence_gate=sequence, methods_gate=methods,
            script_sha256=script_hash,
        )
        BUILD.mkdir(parents=True, exist_ok=True)
        tools = {}
        version_flags = {
            "tectonic": "--version", "pdftotext": "-v",
            "pdftoppm": "-v", "pdfinfo": "-v",
        }
        for name, flag in version_flags.items():
            executable = shutil.which(name)
            if executable is None:
                raise FileNotFoundError(f"Required executable not found: {name}")
            executable = str(Path(executable).resolve(strict=True))
            version = run_command([executable, flag], commands).strip()
            if not version:
                raise ValueError(f"No version output from {name}")
            tools[name] = {"path": executable, "version_output": version}
        payload["tools"] = tools

        compiler_log = BUILD / "tectonic-output.txt"
        run_command(
            [tools["tectonic"]["path"], "-X", "compile", "--keep-logs",
             "--keep-intermediates", "main.tex"],
            commands, output_path=compiler_log,
        )
        pdf = REPORT / "main.pdf"
        pdf_result = artifact(pdf)
        with pdf.open("rb") as handle:
            if not handle.read(5) == b"%PDF-":
                raise ValueError("Compiler output is not a PDF")

        text_path = BUILD / "main.txt"
        run_command(
            [tools["pdftotext"]["path"], "-layout", str(pdf), str(text_path)],
            commands,
        )
        # A fresh directory prevents old pages or different zero-padding from
        # being mistaken for artifacts of this run; old renderings are retained.
        render_directory = BUILD / ("compile-pages-" + uuid4().hex)
        render_directory.mkdir(exist_ok=False)
        run_command(
            [tools["pdftoppm"]["path"], "-png", "-r", "140", str(pdf),
             str(render_directory / "page")],
            commands,
        )
        info = run_command([tools["pdfinfo"]["path"], str(pdf)], commands)
        match = re.search(r"^Pages:\s+(\d+)\s*$", info, re.MULTILINE)
        if match is None or int(match.group(1)) <= 0:
            raise ValueError("pdfinfo did not report a positive page count")
        page_count = int(match.group(1))
        pdf_result["pages"] = page_count

        numbered_pages = {}
        for path in render_directory.glob("page-*.png"):
            match = re.fullmatch(r"page-(\d+)\.png", path.name)
            if match is None:
                continue
            number = int(match.group(1))
            if not 1 <= number <= page_count:
                continue
            if number in numbered_pages:
                raise ValueError(f"Duplicate raster for page {number}")
            with path.open("rb") as handle:
                if handle.read(8) != b"\x89PNG\r\n\x1a\n":
                    raise ValueError(f"Invalid PNG raster: {path}")
            numbered_pages[number] = {"page": number, **artifact(path)}
        if set(numbered_pages) != set(range(1, page_count + 1)):
            raise ValueError("Current render does not contain every PDF page")
        page_images = [numbered_pages[number] for number in range(1, page_count + 1)]
        text_result = artifact(text_path)
        log_result = artifact(compiler_log)

        # Compare the entire mapping to detect added/deleted source files too.
        if source_snapshot() != inputs:
            raise ValueError("Report source set changed during compilation")
        assert_unchanged(sequence)
        assert_unchanged(methods)
        if sha256(SCRIPT) != script_hash:
            raise ValueError("Compilation script changed during execution")
        for output in [pdf_result, text_result, log_result, *page_images]:
            if sha256(ROOT / output["path"]) != output["sha256"]:
                raise ValueError(f"Build artifact changed during execution: {output['path']}")

        payload.update(
            status="passed", completed_at=now(), pdf=pdf_result,
            pdf_text=text_result, compiler_log=log_result,
            page_images=page_images,
            limitations=[
                "Does not inspect scientific content or compiler warnings",
                "Does not perform per-page visual review",
                "External TeX packages and fonts are not vendored by this record",
                "Page rasters and build intermediates remain local project evidence",
            ],
        )
        write_record(payload)
        print(json.dumps({
            "status": "passed", "record": str(RECORD.relative_to(ROOT)),
            "pages": page_count, "pdf_sha256": pdf_result["sha256"],
        }, ensure_ascii=False))
        return 0
    except (Exception, KeyboardInterrupt) as error:
        payload.update(
            status="failed", completed_at=now(),
            error={"type": type(error).__name__, "message": str(error)},
        )
        write_record(payload)
        print(json.dumps({"status": "failed", "error": payload["error"]},
                         ensure_ascii=False), file=sys.stderr)
        return 130 if isinstance(error, KeyboardInterrupt) else 1


if __name__ == "__main__":
    raise SystemExit(main())
