#!/usr/bin/env python3
"""Refresh the generated protocol/reference block from local evidence tables."""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
START = "% BEGIN GENERATED SOURCE RECORDS\n"
END = "% END GENERATED SOURCE RECORDS\n"


def text_escape(value):
    return value.replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")


def online(key, author, title, url, accessed="2026-10-23", note=None, year=None):
    fields = [("author", "{" + text_escape(author) + "}"),
              ("title", "{" + text_escape(title) + "}"),
              ("url", url), ("urldate", accessed)]
    if note:
        fields.append(("note", text_escape(note)))
    if year:
        fields.append(("year", str(year)))
    return "@online{" + key + ",\n" + ",\n".join(
        "  " + name + " = {" + value + "}" for name, value in fields) + "\n}\n"


def main():
    with (ROOT / "docs/协作/方法来源.tsv").open() as handle:
        methods = list(csv.DictReader(handle, delimiter="\t"))
    with (ROOT / "资料/B来源索引.tsv").open() as handle:
        sources = {row["source_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    entries = []
    for key, source_id, author, title in [
        ("TFEBRef", "S01", "NCBI RefSeq", "Homo sapiens TFEB transcript variant 1, NM_007162.3"),
        ("VectorRef", "S04", "GenBank", "Cloning vector pEGFP-N1, complete sequence, U55762.1"),
        ("VectorManual", "S05", "Clontech Laboratories", "pEGFP-N1 Vector Information, PT3027-5 / PR93634"),
        ("IDTTemplate", "S08", "Integrated DNA Technologies", "gBlocks and gBlocks HiFi Gene Fragments"),
        ("Primer3API", "S09", "Primer3-py developers", "Primer3-py 2.4.0 thermodynamic bindings"),
    ]:
        row = sources[source_id]
        entries.append(online(key, author, title, row["url_or_doi"],
            note="Manufacturer document dated 19 March 1999, accessed through a third-party mirror" if source_id == "S05" else None,
            year=1999 if source_id == "S05" else None))
    for row in methods:
        if row["source_id"] == "M17":
            continue  # The original paper already has its article citation.
        author = "Eurofins Genomics" if row["source_id"] in {"M13", "M14", "M15"} else "New England Biolabs"
        entries.append(online("Method" + row["source_id"], author, row["title"], row["url"], row["accessed_date"]))
    path = ROOT / "报告/references.bib"
    original = path.read_text()
    if START in original:
        before, remainder = original.split(START, 1)
        _, after = remainder.split(END, 1)
    else:
        before, after = original.rstrip() + "\n\n", ""
    path.write_text(before + START + "\n".join(entries) + END + after)
    print(f"Updated {len(entries)} generated reference records; preserved the manual bibliography block.")


if __name__ == "__main__":
    main()
