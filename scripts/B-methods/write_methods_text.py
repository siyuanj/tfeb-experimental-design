#!/usr/bin/env python3
"""Render the English P6 handoff from the verified parameter and primer tables."""
import csv
import json
from pathlib import Path
from sequence_checks import sha256

R = Path(__file__).resolve().parents[2]
P = json.loads((R / "设计/methods_parameters.json").read_text())
V = json.loads((R / "docs/协作/方法独立核验.json").read_text())
# This prose contains deliberate B-v1 constants and must fail on changed inputs.
# A later design must review the prose and update this signature explicitly.
if sha256(R / "设计/methods_parameters.json") != "9eacf2114c3993d6511795852ee54c4e652ebf560b0f6159eea1b844b071e32c":
    raise RuntimeError("B-v1 parameters changed; review the narrative before rendering")
if V["status"] != "passed" or not all(V["checks"].values()):
    raise RuntimeError("P6 numerical verification is not passed")
for name, expected in V["inputs"].items():
    if sha256(R / name) != expected:
        raise RuntimeError(f"Changed input: {name}")
for name, expected in V["outputs"].items():
    if sha256(R / name) != expected:
        raise RuntimeError(f"Changed P6 table: {name}")
with (R / "设计/sanger_primers.tsv").open() as stream:
    primers = list(csv.DictReader(stream, delimiter="\t"))
L = P["ligation"]
text = f'''# English P6 methods and expected screening results

Design version: **B-v1, pEGFP-N1-TFEB-Bv1**. This is the method-line handoff for integration into the existing report, not a separate report. All steps describe proposed experiments. Sequence checks are computational; no amplification, transformation, Sanger sequencing or cellular testing has been performed. Source keys M01–M20 map to `docs/协作/方法来源.tsv`; candidate KpnI source M03 is historical and is not part of the final method.

## Materials and sequence-defined construct

Use a planned synthetic linear dsDNA template containing the unmodified 1,431 nt human TFEB CDS from NM_007162.3, including its native TGA. No codon optimization is planned. The supplier QC and delivered sequence must be checked before use; the material has not been acquired. The reverse cloning primer excludes TGA to express **TFEB–EGFP**, with TFEB at the N terminus and EGFP at the C terminus. The empty pEGFP-N1 backbone is U55762.1 (4,733 bp). The proposed construct is 6,134 bp and encodes TFEB 476 aa–GDPPVAT linker 7 aa–EGFP 239 aa, followed by the vector TAA. Its 722-aa sequence has a calculated average molecular mass of 80.425 kDa; this does not predict exact gel migration. The historical TFEB–EGFP construct supports the tag orientation but is not claimed to be sequence-identical to this design [M17].

The selected workflow is Q5 High-Fidelity 2X Master Mix (NEB M0492) → HindIII-HF R3104/BamHI-HF R3136 in rCutSmart → Monarch cleanup/gel extraction → T4 DNA Ligase M0202 → NEB 5-alpha C2987 (DH5α derivative) → colony screening → Monarch T1110 miniprep → Sanger primer walking. OneTaq M0482 is used only for screening and its products never enter the construct.

## Cloning primers and high-fidelity PCR

The primers below are written 5′→3′; spaces separate protection, recognition site, deliberate addition, and template-annealing sequence.

| Primer | Complete sequence partitioned 5′→3′ | Annealing Tm |
|---|---|---|
| TFEB-HindIII-F | ACTGAC AAGCTT ACC ATGGCGTCACGCATAGGGTTGC | 69.498°C |
| TFEB-BamHI-R | TATGAC GGATCC CC CAGCACATCGCCCTCCTCCATG | 69.050°C |

The forward ACC supplies the initiation context; reverse-primer CC introduces coding-strand GG to maintain the fusion frame. Six protective bases are deliberately used for each enzyme. NEB's end-cleavage table assigns both enzymes +++ at three flanking base pairs; +++ means 50–100% under its oligonucleotide test conditions, so six bases are a conservative choice, not proof of complete digestion [M02, M04, M19].

Both annealing segments are 22 nt. Tm values use Primer3-py 2.4.0/SantaLucia with 50 mM monovalent ions, 2.0 mM Mg2+, 0.8 mM total dNTP and 500 nM primer. The monovalent-ion value is an explicit modeling assumption, not the proprietary Q5 composition. Full tailed-oligo Tm is not used to select the initial annealing temperature. Neither transcriptome-wide specificity nor an optimal wet-lab annealing temperature is claimed.

| Component | Volume in 50 µL | Final/target concentration |
|---|---:|---|
| Q5 M0492 2X Master Mix | 25.0 µL | 1X; Mg2+ 2.0 mM; each dNTP 200 µM |
| Each cloning primer, 10 µM | 2.5 µL each | 0.5 µM each |
| Synthetic template dilution | 1.0 µL | Target 1 ng total from 1 ng/µL stock |
| Nuclease-free water | 19.0 µL | Total 50.0 µL |

The template concentration is a proposed dilution target, not a measured value. Calculate dilution from measured DNA concentration. A no-template control replaces the template with water. Assemble on ice and use 98°C 30 s; 30 cycles of 98°C 10 s, 65°C 20 s, 72°C 45 s; then 72°C 2 min. The 45 s extension is the 1,457 bp expected product × 30 s/kb rounded upward. Initial 65°C annealing and a 62/65/68°C troubleshooting gradient are design choices rather than Q5-calculator outputs or measured optima [M01].

## Electrophoresis and purification

Analyze 5 µL PCR product on 1.0% agarose in 1X TAE with 1X loading dye. A starting field of 5 V/cm electrode separation is a design choice; stop when relevant bands resolve before the tracking front leaves the gel. Use 0.5 µg NEB N3232 ladder, a compatible DNA stain at its recommended concentration, and minimal illumination of DNA to be cloned. The marker contains 10,002, 8,001, 6,001, 5,001, 4,001, 3,001, 2,000, 1,500, 1,000, 517 and 500 bp fragments; the 517/500 pair normally co-migrates. Uncut plasmids are not sized against a linear-fragment calibration [M16].

For a single 1,457 bp band, use Monarch Spin PCR & DNA Cleanup T1130: sample 20–100 µL, five volumes prepared BZ, S1A column, two 200 µL WZ washes, and 20 µL EY elution after a 1 min wait; centrifuge according to the current kit protocol (16,000 ×g, room temperature). Cleanup removes primers, dNTPs, salts and polymerase. If extra PCR bands occur, excise the target and use T1120: four volumes BY per gel-slice volume, dissolve near 50°C, load S1A, wash twice with 200 µL WZ and elute in 20 µL EY. Prepare alcohol-containing buffers using the current bottle label. Measure recovered DNA; do not assume a yield or recovery fraction [M09, M10, M18].

## Double digestion and ligation

Digest PCR product and empty vector separately. Each 50 µL reaction contains 500 ng DNA (target 10.0 µL at 50 ng/µL), 5.0 µL 10X rCutSmart, 1.0 µL HindIII-HF at 20 U/µL, 1.0 µL BamHI-HF at 20 U/µL and 33.0 µL water. Recalculate DNA/water volumes from measured concentration; if insufficient material is available, repeat or scale the preparation without inventing yield. Incubate at 37°C for 60 min. This keeps enzyme-stock glycerol at 2%. Include uncut and single-digest analytical controls if troubleshooting. BamHI-HF cannot be heat-inactivated; purify before ligation. For direct preparative or diagnostic gel loading of unpurified BamHI-HF digests, use SDS-containing Purple Gel Loading Dye B7024 (6X) at 1X: add one volume dye to five volumes sample. This dissociates residual enzyme-bound DNA; do not substitute the no-SDS B7025 product. Add loading dye only to gel samples, not to ligation. If SYBR Safe or GelRed is used, stain after electrophoresis as directed for B7024; retain the planned gel recovery or column cleanup before ligation [M02, M05, M19, M20].

Gel-purify the 4,695 bp backbone, separating residual circular vector and the excised 38 nt MCS. Purify the digested insert with T1130; its nominal top-strand length is 1,439 nt, while discarded terminal pieces are 7 and 11 nt. The 38/7/11 nt species are not expected to be usefully visible on this 1% agarose gel. Both retained fragments have equal-length strands: pydna reports full duplex spans of 4,699 and 1,443 bases because of staggered ends. Report nominal 4,695/1,439 bp sizes for gels and use mean strand length for molarity; the two definitions coincide for this pair. Distinct 5′ AGCT and GATC ends support directional ligation. Routine dephosphorylation is omitted; vector-only background is used to detect incomplete digestion/circular-vector carryover.

Set vector to 20 fmol and insert to 60 fmol (3:1). For mean strand length L, mass (ng) = fmol × L × 0.00066, giving **{L['vector_ng']:.3f} ng backbone and {L['insert_ng']:.4f} ng insert**. A practical 20 µL recipe uses the target stocks below; concentrations remain measurements to be made [M06].

| Component | Volume (rounded to 0.01 µL) |
|---|---:|
| Backbone at target 20 ng/µL | {L['vector_uL']:.2f} µL |
| Insert at target 10 ng/µL | {L['insert_uL']:.2f} µL |
| 10X T4 ligase buffer | 2.00 µL |
| T4 DNA ligase, 400 U/µL | 1.00 µL |
| Water | {L['water_uL']:.2f} µL |
| Total | 20.00 µL |

The unrounded calculation gives {L['total_DNA_ng_per_uL']:.3f} ng/µL total DNA, within the manufacturer's 1–10 ng/µL starting range. Mix fully thawed ATP-containing buffer; incubate at 20–25°C for 10 min, then 65°C for 10 min. Prepare a matched vector-only ligation by replacing insert volume with water. No PNK is needed after the planned restriction digestion.

## Transformation and antibiotic selection

Add 2 µL ligation to 50 µL NEB 5-alpha C2987 competent cells. Use ice 30 min → 42°C 30 s → ice 5 min → 950 µL SOC → 37°C, 250 rpm, 60 min. Plate 100 µL undiluted recovery and 1:10/1:100 SOC dilutions on LB-kanamycin, then incubate at 37°C overnight. Use **kanamycin 30 µg/mL** consistently for bacterial selection under the chosen C2987 guidance; G418 is reserved for mammalian selection. The 2 µL input contains approximately 11.90 ng DNA at the design concentration [M07, M08].

Process the vector-only ligation, no-DNA control and pUC19 transformation control in parallel. pUC19 is plated on ampicillin 100 µg/mL, not kanamycin. These controls distinguish cloning background, contamination/selection failure and competent-cell handling; no efficiency or successful-clone fraction is assumed.

## Colony PCR, miniprep and diagnostic digest

Screen eight isolated colonies initially and retain numbered patches. Use trace colony material in 25 µL OneTaq M0482: 12.5 µL 2X mix, 0.5 µL of each 10 µM screening primer and 11.5 µL water. Final primers are 0.2 µM. Program: 94°C 5 min for colony lysis; 30 cycles of 94°C 30 s, 58°C 30 s, 68°C 160 s; then 68°C 5 min. The 160 s extension covers the longest modeled 2,608 bp product at ≥1 min/kb. Initial 58°C and a 55/58/61°C optimization gradient are design choices. Calculated pair Tm values under an explicit 50 mM monovalent/1.8 mM Mg/0.8 mM dNTP/200 nM primer model are 60.04/62.55°C for S-F1/S-R1 and 66.09/62.55°C for S-F2/S-R1; this is not the Sanger Tm model [M12].

| Assay | Correct construct | Empty vector | Reversed-TFEB model | Retained-stop model |
|---|---:|---:|---:|---:|
| S-F1 + S-R1 flanking PCR | 2,605 bp | 1,204 bp | 2,605 bp | 2,608 bp |
| S-F2 + S-R1 orientation PCR | 2,093 bp | No specific product | No specific product | 2,096 bp |
| HindIII-HF + BamHI-HF digest | 4,695 + 1,439 bp | 4,695 + 38 bp | 4,695 + 1,439 bp | 4,695 + 1,442 bp |

The reversed-TFEB sequence is a hypothetical failure challenge, not a predicted product of normal ligation with incompatible ends. Flanking PCR and the cloning-enzyme digest cannot establish orientation; three-base stop differences are not resolved by the proposed gel. Include empty-vector and no-template PCR controls. A negative orientation assay is provisional until an orientation-positive control is available: in the first construction round, analyze flanking-positive minipreps by sequencing if all orientation reactions are negative. An accepted clone can serve as the positive control in later rounds.

Grow up to three screening-positive colonies separately in 3 mL LB-kanamycin at 37°C for 12–16 h with aeration; if fewer are recovered, analyze those available. For routine minipreps use 1.5 mL culture per T1110 prep: resuspend in 200 µL B1, add 200 µL B2 for 1 min with gentle inversion, neutralize with 400 µL B3 for 2 min, clarify, bind S2D, wash with 200 µL BZ then 400 µL WZ, and elute in ≥30 µL EY. B1 contains RNase A. For DNA intended for transfection, add the manufacturer's 5 min BZ incubation; this is not an endotoxin-free certification. Digest 500 ng miniprep using the defined 50 µL HindIII/BamHI conditions and analyze an aliquot. Size agreement cannot exclude substitutions, short indels or every rearrangement [M11, M18].

## Sanger primer walking and success criteria

Sequence each candidate with the ten primers below in separate reactions. For the selected Eurofins Standard Read service, the 6.134 kb plasmid falls in the 150–250 ng/µL stock range. A design target is 200 ng/µL: for ten reads from one sample, the separate-template format calls for 10 + 9×5 = 55 µL DNA. If one miniprep is insufficient, pool independently prepared aliquots from the same single-colony culture and measure/concentrate them; do not pool different clones. Supply each primer separately at 10 µM, at least 10 µL for its first read, with 5 µL for each additional use. These are submission targets, not claimed yields [M13].

Sanger-primer Tm values below use Primer3/SantaLucia, 50 mM monovalent ions, zero added Mg/dNTP and 50 nM oligo; they are a reproducible design comparison, not proprietary sequencing chemistry. Selection required 18–25 nt, 35–60% GC, Tm 50–62°C, no run of four identical bases, a terminal G/C without three consecutive terminal G/C bases, and modeled hairpin/homodimer Tm below 40°C. Each primer has one exact full-length site in the construct; genome-wide specificity was not tested.

The conservative usable read window is extended bases **51–650 after the primer 3′ end** (600 bases), consistent with excluding the initial low-quality portion and staying within a standard-read length. Walks were selected with at least 50 nt overlap between consecutive usable windows. The complete target **602–2799 (2,198 nt)** includes 20 upstream vector bases, the 5′ junction/initiation context, the entire 2,169 nt fusion ORF, both fusion interfaces and final stop. Every target base has at least one forward and one reverse planned read. Binding bases inside the target are covered by other primers. The remainder of the vector is not fully sequenced by this plan [M13–M15].

| ID | Sequence 5′→3′ | Binding coordinates/strand | Planned usable coordinates | Model Tm |
|---|---|---|---|---|
'''
for p in primers:
    text += f"| {p['primer_id']} | {p['sequence_5to3']} | {p['binding_start_1based']}–{p['binding_end_1based']} ({p['direction']}) | {p['usable_start_1based']}–{p['usable_end_1based']} | {float(p['tm_C']):.2f}°C |\n"
text += '''
Accept a clone only after actual chromatograms provide unambiguous, concordant bidirectional Q20-or-better calls across the target, no unexpected nucleotide changes, and the predicted restriction pattern. Inspect all mismatches and junctions manually after aligning quality-trimmed reads to the exact B-v1 sequence. Q20 describes a base-call error probability; it is not a guarantee that a molecule is error-free. Low-quality regions, single-strand gaps and ambiguous peaks require additional reads/primers. Reproducible mixed peaks require re-isolation of a single colony and fresh prep. An unambiguous unwanted coding/junction mutation rejects that clone; a shared error in independent clones prompts investigation of the template and design. Predicted coverage never substitutes for actual trace quality.

## Troubleshooting and limits

No cloning PCR band prompts checks of template identity, dilution, primer matching and annealing conditions; extra bands prompt a gradient and target-band purification. High vector-only background prompts checks of digestion and circular-vector carryover. No colonies in both test and positive transformation control points to handling/selection rather than proving failed ligation. A correct gel pattern with sequence disagreement fails acceptance. Sequence-accepted DNA is the input to the proposed HeLa experiment; its cellular expression, localization and biological behavior remain to be tested.

The supporting artifacts are `设计/expected_fragments.tsv`, `设计/sanger_primers.tsv`, `设计/sanger_coverage.tsv`, `设计/methods_parameters.json` and `docs/协作/方法独立核验.json`. Predicted bands and all schematic figures must be labeled as expected, not measured.
'''
output = R / "docs/协作/方法正文.md"
# The file is exclusively owned by this method line; read current bytes first.
if output.exists():
    output.read_bytes()
output.write_text(text)
print(f"Rendered {output.relative_to(R)} ({len(primers)} primers)")
