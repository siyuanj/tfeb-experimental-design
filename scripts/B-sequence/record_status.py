#!/usr/bin/env python3
"""Prepend a status entry with backup and concurrent-write protection."""
import argparse
import hashlib
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parents[2]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("entry", type=Path)
    args = parser.parse_args()
    path = ROOT / "docs/STATUS.md"
    before = path.read_bytes()
    heading, separator, body = before.partition(b"\n")
    if not separator or not heading.startswith(b"# "):
        raise ValueError("Missing status title")
    stamp = datetime.now(ZoneInfo("America/New_York")).strftime("%Y%m%d-%H%M%S-%f")
    backup = ROOT / "docs/备份" / f"STATUS-{stamp}.md"
    backup.parent.mkdir(parents=True, exist_ok=True)
    backup.write_bytes(before)
    addition = args.entry.read_bytes().strip() + b"\n\n"
    after = heading + b"\n\n" + addition + body
    if path.read_bytes() != before:
        raise RuntimeError("Concurrent status change; reread and merge")
    path.write_bytes(after)
    actual = path.read_bytes()
    if actual != after or not actual.endswith(body):
        raise RuntimeError("Status changed during write; inspect, do not restore blindly")
    print({"backup": str(backup.relative_to(ROOT)), "preserved_body_sha256": hashlib.sha256(body).hexdigest()})

if __name__ == "__main__":
    main()
