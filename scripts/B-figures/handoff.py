"""Read-only validation of stage gates used by report generation."""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_gate(relative_path):
    path = ROOT / relative_path
    gate = json.loads(path.read_text())
    if gate.get("status") != "ready":
        raise ValueError(f"Stage is not ready: {relative_path}")
    files = gate.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError(f"No bound input hashes: {relative_path}")
    files = dict(files)
    for filename, expected in gate.get("inputs", {}).items():
        if filename in files and files[filename] != expected:
            raise ValueError(f"Conflicting bound hashes: {filename}")
        files[filename] = expected
    for filename, expected in files.items():
        file_path = (ROOT / filename).resolve(strict=True)
        if not file_path.is_relative_to(ROOT):
            raise ValueError(f"Handoff path escapes this repository: {filename}")
        actual = sha256(file_path)
        if actual != expected:
            raise ValueError(f"Handoff input changed: {filename}")
    return gate, {"path": relative_path, "sha256": sha256(path), "files": files}


def assert_unchanged(snapshot):
    if sha256(ROOT / snapshot["path"]) != snapshot["sha256"]:
        raise ValueError(f"Handoff changed while rendering: {snapshot['path']}")
    for filename, expected in snapshot["files"].items():
        if sha256(ROOT / filename) != expected:
            raise ValueError(f"Input changed while rendering: {filename}")
