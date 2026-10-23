#!/usr/bin/env python3
"""Render an immutable PDF snapshot for manual review, never the live build."""

import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def main():
    source = ROOT / "报告/main.pdf"
    data = source.read_bytes()
    pdf_hash = sha256(data)
    directory = ROOT / "docs/核查/B" / f"main-final-review-{pdf_hash[:12]}"
    directory.mkdir(parents=True, exist_ok=True)
    snapshot = directory / "snapshot.pdf"
    if snapshot.exists() and snapshot.read_bytes() != data:
        raise ValueError("Existing snapshot differs; do not overwrite review evidence")
    snapshot.write_bytes(data)
    if source.read_bytes() != data:
        raise ValueError("Source changed while taking the snapshot")
    info = subprocess.run(["pdfinfo", str(snapshot)], check=True, capture_output=True, text=True)
    page_match = re.search(r"^Pages:\s+(\d+)", info.stdout, re.M)
    if not page_match or int(page_match[1]) < 1:
        raise ValueError("PDF has no readable pages")
    render = subprocess.run(
        ["pdftoppm", "-r", "110", "-png", str(snapshot), str(directory / "page")],
        capture_output=True, text=True,
    )
    (directory / "render.stderr.txt").write_text(render.stderr)
    if render.returncode or render.stderr.strip():
        raise ValueError(f"Renderer returned {render.returncode}: {render.stderr}")
    pages = sorted(directory.glob("page-*.png"))
    numbers = [int(path.stem.split("-")[-1]) for path in pages]
    if numbers != list(range(1, int(page_match[1]) + 1)):
        raise ValueError("Page images do not cover the complete PDF")
    if source.read_bytes() != data or snapshot.read_bytes() != data:
        raise ValueError("PDF changed during rendering; snapshot is not final")
    payload = {
        "status": "rendered_not_visually_reviewed",
        "pdf_sha256": pdf_hash,
        "snapshot": str(snapshot.relative_to(ROOT)),
        "page_count": len(pages),
        "dpi": 110,
        "script_sha256": sha256(Path(__file__).read_bytes()),
        "pages": [{"page": n, "image": str(p.relative_to(ROOT)), "sha256": sha256(p.read_bytes())}
                  for n, p in zip(numbers, pages)],
    }
    (directory / "render-manifest.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "directory": str(directory),
                      "page_count": len(pages), "pdf_sha256": pdf_hash}, ensure_ascii=False))


if __name__ == "__main__":
    main()
