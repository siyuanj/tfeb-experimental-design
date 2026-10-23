"""Independent sequence primitives for method design (not an assembly library)."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import primer3


def sha256(path: Path) -> str:
    """Return a file-content digest without parsing or normalizing bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rc(sequence: str) -> str:
    """Reverse-complement an unambiguous DNA sequence."""
    if not sequence or set(sequence.upper()) - set("ACGT"):
        raise ValueError("Expected nonempty unambiguous DNA")
    return sequence.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def occurrences(sequence: str, motif: str, circular: bool = False) -> list[int]:
    """Return zero-based motif starts, including origin-spanning matches."""
    if not motif or len(motif) > len(sequence):
        raise ValueError("Invalid motif length")
    target = sequence + sequence[:len(motif) - 1] if circular else sequence
    return [match.start() for match in re.finditer(f"(?={motif})", target)
            if match.start() < len(sequence)]


def binding_sites(sequence: str, primer: str, circular: bool = True) -> list[dict]:
    """Find exact full-length primer sites on both strands."""
    sites = []
    for strand, query in ((1, primer), (-1, rc(primer))):
        for start in occurrences(sequence, query, circular):
            sites.append({"start0": start, "end0": start + len(primer),
                          "strand": strand})
    return sites


def pcr_products(sequence: str, forward: str, reverse: str,
                 circular: bool = True, maximum: int = 4000) -> list[dict]:
    """Enumerate inward-facing exact-site amplicons below a length limit.

    Both primers are considered on both strands; this detects unexpected products
    when an insert reverses or a primer has an additional full-length match.
    No mismatches, secondary structures or amplification efficiencies are modeled.
    """
    n = len(sequence)
    output = []
    for fprimer, rprimer in ((forward, reverse), (reverse, forward)):
        for fsite in occurrences(sequence, fprimer, circular):
            for rsite in occurrences(sequence, rc(rprimer), circular):
                distance = (rsite - fsite) % n if circular else rsite - fsite
                size = distance + len(rprimer)
                if distance >= len(fprimer) and size <= min(maximum, n):
                    output.append({"forward_primer": fprimer,
                                   "reverse_primer": rprimer,
                                   "forward_start_1based": fsite + 1,
                                   "reverse_binding_start_1based": rsite + 1,
                                   "length_bp": size})
    return output


def digest_literal(sequence: str, enzymes: list[tuple[str, str, int]],
                   circular: bool = True) -> dict:
    """Digest from literal motifs and top-strand cut offsets.

    Enzyme tuples contain name, recognition sequence, and bases retained before
    the top-strand cut. Lengths are the top-strand cut-to-cut convention and are
    not pydna's end-to-end duplex span or average of strand lengths.
    """
    n = len(sequence)
    cut_events = []
    for name, motif, offset in enzymes:
        if motif != rc(motif):
            raise ValueError("This independent checker expects palindromic sites")
        for start in occurrences(sequence, motif, circular):
            cut = start + offset
            cut_events.append({"enzyme": name, "site_start_1based": start + 1,
                               "cut_top_0based_boundary": cut % n if circular else cut})
    cuts = sorted(set(e["cut_top_0based_boundary"] for e in cut_events))
    if circular:
        fragments = ([n] if len(cuts) == 1 else
                     [((cuts[(i + 1) % len(cuts)] - c) % n) for i, c in enumerate(cuts)])
    else:
        cuts = [0] + cuts + [n]
        fragments = [b - a for a, b in zip(cuts, cuts[1:])]
    return {"cuts": cut_events, "fragments_bp": sorted(fragments, reverse=True),
            "uncut": not cut_events, "length_definition": "top-strand cut-to-cut"}


def sanger_metrics(primer: str) -> dict:
    """Compute descriptive Sanger-primer metrics under a declared salt model.

    Salt input is a reproducible comparison model, not Eurofins chemistry.
    Structure screens are calculated at 37 C; Tm cutoffs are design choices.
    """
    conditions = dict(mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=50)
    structure = dict(**conditions, temp_c=37)
    return {"length_nt": len(primer),
            "gc_percent": 100 * (primer.count("G") + primer.count("C")) / len(primer),
            "tm_C": primer3.calc_tm(primer, **conditions, tm_method="santalucia",
                                    salt_corrections_method="santalucia"),
            "hairpin_tm_C": primer3.calc_hairpin(primer, **structure).tm,
            "homodimer_tm_C": primer3.calc_homodimer(primer, **structure).tm}


def choose_primer(sequence: str, window: tuple[int, int], strand: int,
                  preferred_start: int, name: str, must_cover_0: int | None = None) -> dict:
    """Choose one sequence-exact primer from an explicit start window.

    Ranking targets 58 C and the requested location. Hard constraints follow the
    selected Sanger service ranges, plus conservative structure design screens.
    This screens uniqueness within the final plasmid only, not a transcriptome.
    """
    candidates = []
    for start in range(max(0, window[0]), min(len(sequence) - 18, window[1]) + 1):
        for length in range(18, 26):
            if start + length > len(sequence):
                continue
            if strand == 1:
                first, last = start + length + 50, start + length + 649
            else:
                first, last = start - 650, start - 51
            if first < 0 or last >= len(sequence):
                continue
            if must_cover_0 is not None and not first <= must_cover_0 <= last:
                continue
            primer = sequence[start:start + length]
            if strand == -1:
                primer = rc(primer)
            if re.search(r"(.)\1{3}", primer):
                continue
            if primer[-1] not in "GC" or re.search(r"[GC]{3}$", primer):
                continue
            metrics = sanger_metrics(primer)
            if not (50 <= metrics["tm_C"] <= 62 and 35 <= metrics["gc_percent"] <= 60):
                continue
            if max(metrics["hairpin_tm_C"], metrics["homodimer_tm_C"]) >= 40:
                continue
            sites = binding_sites(sequence, primer)
            if len(sites) != 1 or sites[0]["strand"] != strand:
                continue
            score = abs(metrics["tm_C"] - 58) + abs(start - preferred_start) / 30
            candidates.append((score, start, length, primer, metrics))
    if not candidates:
        raise ValueError(f"No primer meets declared constraints: {name} {window}")
    _, start, length, primer, metrics = min(candidates)
    if strand == 1:
        first, last = start + length + 50, start + length + 649
    else:
        first, last = start - 650, start - 51
    assert 0 <= first <= last < len(sequence)
    return {"primer_id": name, "sequence_5to3": primer,
            "binding_start_1based": start + 1, "binding_end_1based": start + length,
            "direction": "+" if strand == 1 else "-",
            "usable_start_1based": first + 1, "usable_end_1based": last + 1,
            **metrics}
