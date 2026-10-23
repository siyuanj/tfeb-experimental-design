#!/usr/bin/env python3
"""Repeat reviewed audits in an exclusive output directory without editing inputs.

This wrapper reuses the independent sequence and method implementations; it is
not a third independent sequence algorithm. Primer3 values are recalculated
directly, and PDF raster integrity is checked separately from visual review.
Run with the project's existing plasmid interpreter and -B.
"""
import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import primer3

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(relative):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    if len(sys.argv) != 1:
        raise ValueError("Run without arguments")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = ROOT / "docs/核查/B" / ("parent-review-" + stamp)
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "running", "checked_at_utc": stamp,
              "git_commit": subprocess.check_output(
                  ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "script_sha256": sha(Path(__file__)), "python": sys.version,
              "output_directory": str(output.relative_to(ROOT)),
              "limitations": ["Reuses reviewed independent audit implementations",
                              "Thermodynamics repeat the declared Primer3 model",
                              "Raster integrity is not visual inspection",
                              "No wet experiments or independent A protocol validation"]}

    def save(name, value):
        (output / name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")

    try:
        structural = load("scripts/B-audit/final_audit.py")
        structural.OUTPUT = output / "structural.json"
        with contextlib.redirect_stdout(io.StringIO()):
            structural.main()
        chain = json.loads(structural.OUTPUT.read_text())
        snapshot = dict(chain["files"])
        for path in ("scripts/B-sequence/independent_verify.py",
                     "scripts/B-methods/independent_review.py",
                     "scripts/B-supervision/review_handoffs.py",
                     "docs/协作/报告逐页核验.json",
                     "资料/B来源索引.tsv", "docs/10-23 B部分执行计划.md",
                     "docs/B执行进度.md"):
            snapshot[path] = sha(ROOT / path)
        sequence = {"inputs": [], "checks": [], "negative_controls": []}
        load("scripts/B-sequence/independent_verify.py").audit(sequence)
        sequence["status"] = "passed"
        save("sequence.json", sequence)
        process = subprocess.run(
            [sys.executable, "-B", str(ROOT / "scripts/B-methods/independent_review.py")],
            cwd=ROOT, text=True, capture_output=True, check=True)
        methods = json.loads(process.stdout)
        if methods["status"] != "passed":
            raise ValueError("Independent method audit failed")
        save("methods.json", methods)

        checks = json.loads((ROOT / "设计/primer_checks.json").read_text())
        model = checks["parameters"]
        args = dict(mv_conc=model["mv_conc_mM"], dv_conc=model["dv_conc_mM"],
                    dntp_conc=model["dntp_conc_mM_total"], dna_conc=model["dna_conc_nM"])
        structure_args = dict(args, temp_c=model["structure_temperature_C"])
        rows = list(csv.DictReader(io.StringIO((ROOT / "设计/primers.tsv").read_text()), delimiter="\t"))
        thermo = []

        def compare(name, observed, expected, tolerance=1e-6):
            passed = math.isfinite(observed) and abs(observed - expected) <= tolerance
            thermo.append(dict(name=name, observed=observed, expected=expected,
                               tolerance=tolerance, passed=passed))
            if not passed:
                raise ValueError("Thermodynamic mismatch: " + name)

        for row in rows:
            name, oligo = row["name"], row["sequence_5to3"]
            reference = next(item for item in checks["primers"] if item["primer"] == name)
            tm_args = dict(args, tm_method=model["tm_method"],
                           salt_corrections_method=model["salt_corrections_method"])
            compare(name + ":annealing_Tm", primer3.calc_tm(row["annealing_sequence_5to3"], **tm_args),
                    float(row["annealing_tm_C"]), 0.00051)
            compare(name + ":full_Tm", primer3.calc_tm(oligo, **tm_args),
                    float(row["full_oligo_tm_C_not_for_initial_annealing"]), 0.00051)
            for kind in ("hairpin", "homodimer"):
                actual = getattr(primer3, "calc_" + kind)(oligo, **structure_args)
                compare(name + ":" + kind, actual.tm, reference["full_oligo"][kind]["tm"])
        forward, reverse = [row["sequence_5to3"] for row in rows]
        compare("heterodimer", primer3.calc_heterodimer(forward, reverse, **structure_args).tm,
                checks["heterodimer"]["tm"])
        for index, (first, second) in enumerate(((forward, reverse), (reverse, forward))):
            compare("3prime_stability:" + str(index),
                    primer3.calc_end_stability(first, second, **structure_args).tm,
                    checks["directional_3prime_end_stability"][index]["tm"])
        parameters = json.loads((ROOT / "设计/methods_parameters.json").read_text())
        sanger = list(csv.DictReader(io.StringIO((ROOT / "设计/sanger_primers.tsv").read_text()), delimiter="\t"))
        model = parameters["sanger"]["Tm_model"]
        args = dict(mv_conc=model["mv_mM"], dv_conc=model["Mg_mM"],
                    dntp_conc=model["dNTP_mM"], dna_conc=model["primer_nM"])
        for row in sanger:
            name, oligo = row["primer_id"], row["sequence_5to3"]
            compare(name + ":Sanger_Tm", primer3.calc_tm(oligo, **args,
                    tm_method="santalucia", salt_corrections_method="santalucia"), float(row["tm_C"]))
            for kind in ("hairpin", "homodimer"):
                compare(name + ":Sanger_" + kind,
                        getattr(primer3, "calc_" + kind)(oligo, temp_c=37, **args).tm,
                        float(row[kind + "_tm_C"]))
        model = parameters["colony_PCR"]["Tm_model"]
        args = dict(mv_conc=model["mv_mM"], dv_conc=model["Mg_mM"],
                    dntp_conc=model["total_dNTP_mM"], dna_conc=model["primer_nM"])
        oligos = {row["primer_id"]: row["sequence_5to3"] for row in sanger}
        for name, pair in parameters["colony_PCR"]["pairs"].items():
            sequences = [oligos[identifier] for identifier in pair["primers"]]
            for index, oligo in enumerate(sequences):
                compare(name + ":Tm:" + str(index), primer3.calc_tm(oligo, **args,
                        tm_method="santalucia", salt_corrections_method="santalucia"),
                        pair["model_tm_C"][index])
            compare(name + ":heterodimer", primer3.calc_heterodimer(*sequences, temp_c=37, **args).tm,
                    pair["heterodimer_tm_C"])
        save("thermodynamics.json", thermo)

        sources = list(csv.DictReader(io.StringIO((ROOT / "资料/B来源索引.tsv").read_text()), delimiter="\t"))
        source_hashes = {}
        for row in sources:
            if re.fullmatch(r"[0-9a-f]{64}", row["raw_file_sha256"]):
                relative = row["raw_file_path"]
                current = sha(ROOT / relative)
                if current != row["raw_file_sha256"]:
                    raise ValueError("Source ledger hash mismatch: " + row["source_id"])
                snapshot[relative] = current
                source_hashes[row["source_id"]] = dict(path=relative, sha256=current)
        plan = (ROOT / "docs/10-23 B部分执行计划.md").read_text()
        progress = (ROOT / "docs/B执行进度.md").read_text()
        planned = re.findall(r"^- \[[ x]\] (.*)$", plan, re.M)
        tracked = re.findall(r"^- \[[ x]\] (.*)$", progress, re.M)
        if planned != tracked or "- [ ]" in progress:
            raise ValueError("Incomplete or changed P0-P9 checklist")

        def table(relative):
            return list(csv.DictReader(io.StringIO((ROOT / relative).read_text()), delimiter="\t"))

        map_data = json.loads((ROOT / "报告/images/B/plasmid-map.json").read_text())
        pipeline = json.loads((ROOT / "报告/images/B/pipeline-figures.json").read_text())
        features = table("设计/construct_features.tsv")
        expected_features = {(row["name"], int(row["start_1based"]), int(row["end_1based"]),
                              int(row["strand"])) for row in features}
        shown = {(row["label"], row["start_1based"], row["end_inclusive_1based"], row["strand"])
                 for row in map_data["features_displayed"]}
        if shown != expected_features or len(shown) != len(features):
            raise ValueError("Plasmid feature coordinates or directions differ")
        fragments = table("设计/expected_fragments.tsv")
        lane_count = 0
        for panel in pipeline["gel"]["panels"]:
            for lane in panel["lanes"]:
                expected = [int(row["size_bp"]) for row in fragments if
                            row["assay"] == lane["assay"] and row["sample"] == lane["sample"]
                            and int(row["size_bp"]) > 0]
                if lane["bands_bp"] != expected:
                    raise ValueError("Plotted bands differ from fragment table")
                lane_count += 1
        if pipeline["gel"]["marker_bp"] != [int(row["size_bp"]) for row in fragments if row["assay"] == "marker"]:
            raise ValueError("Plotted marker differs from reference table")
        build = json.loads((ROOT / "docs/协作/报告编译核验.json").read_text())
        for relative, expected in build["inputs"].items():
            if sha(ROOT / relative) != expected:
                raise ValueError("Compilation input mismatch: " + relative)
        for entry in (build["pdf"], build["pdf_text"], build["compiler_log"], *build["page_images"]):
            if sha(ROOT / entry["path"]) != entry["sha256"]:
                raise ValueError("Compilation output mismatch: " + entry["path"])

        pages = json.loads((ROOT / "docs/协作/报告逐页核验.json").read_text())
        if pages["reviewed_pdf_sha256"] != sha(ROOT / "报告/main.pdf"):
            raise ValueError("Raster review targets a different PDF")
        if [page["page"] for page in pages["pages"]] != list(range(1, chain["pdf_pages"] + 1)):
            raise ValueError("Missing PDF raster page")
        for page in pages["pages"]:
            if sha(ROOT / page["image"]) != page["sha256"]:
                raise ValueError("Raster hash mismatch")
        for path, expected in snapshot.items():
            if sha(ROOT / path) != expected:
                raise ValueError("Input changed during review: " + path)
        result.update(status="passed", structural_checks=len(chain["checks"]),
                      bound_files=len(chain["files"]), sequence_groups=len(sequence["checks"]),
                      negative_controls=sequence["negative_controls"],
                      method_checks=methods["checks"], thermodynamic_checks=len(thermo),
                      checklist_items=len(planned), source_ledger_entries=len(sources),
                      source_hashes=source_hashes,
                      plotted_features=len(shown), plotted_lanes=lane_count,
                      compilation_inputs=len(build["inputs"]),
                      pdf_pages=chain["pdf_pages"], pdf_sha256=pages["reviewed_pdf_sha256"],
                      page_rasters=pages["pages"], input_hashes=snapshot,
                      evidence_files={name: sha(output / name) for name in
                                      ("structural.json", "sequence.json", "methods.json", "thermodynamics.json")})
    except Exception as error:
        result.update(status="failed", error=str(error), error_type=type(error).__name__)
        save("summary.json", result)
        raise
    save("summary.json", result)
    print(json.dumps({key: value for key, value in result.items() if key not in
                      {"input_hashes", "page_rasters", "negative_controls", "python"}}, ensure_ascii=False))


if __name__ == "__main__":
    main()
