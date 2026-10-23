# TFEB report figures

These scripts use the existing local plasmid environment; they do not install dependencies or launch a web application. All drawings are predictions or schematics, never experimental images.

`render_plasmid_map.py` is adapted from the course renderer with original SHA-256 `e79564a42038eb9440236b2f09672de885f2ded1ba280d0b8a2ae3a31895b6c7`. It draws every located, non-omitted annotation from its actual coordinates. A multipart `join` is split into its real parts; only the last part in biological extraction order carries the terminal arrow. Mixed-strand joins, fuzzy positions, and unsupported location operators are rejected. Enzyme labels are the 1-based nucleotide immediately after reference-strand cleavage, not recognition-site starts.

`check_renderer.py` checks positive/reverse origin-spanning joins, an internal gap, rejection of mixed-strand joins, and EcoRI cleavage-label conventions using small synthetic inputs. It does not validate TFEB sequence design. Run from the repository root:

An explicit annotation note `no directional arrow assigned` suppresses an arrow even if a GenBank reader returns `strand=+1`. This handles the known loss of unstranded (`0`) intervals during GenBank round-tripping. The manifest retains both the parsed strand and the displayed strand.

```bash
'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' scripts/B-figures/check_renderer.py
```

Final figure generation is gated by the sequence and methods handoff manifests. `handoff.py` checks both `files` and upstream `inputs`, rejects non-ready stages, and checks that the snapshot did not change during generation. A retracted P6 gate was observed to stop rendering before any new output was saved.

From the repository root, use the existing interpreter (the quoted path above) for these scripts in order:

1. `write_construct_values.py`: export the accepted construct values to LaTeX macros.
2. `sync_references.py`: retain the manual bibliography and replace its generated source-record block.
3. `render_accepted_map.py`: render every accepted annotation and cut label; retain input, script and output hashes.
4. `render_pipeline_figures.py`: render the method workflow and predicted gel lanes from accepted tables and parameters.
5. `compile_report.py`: bind all source/image inputs, compile, extract text and render pages; verify that the inputs and accepted gates remained unchanged throughout the build.
6. `check_report_consistency.py`: check figure hashes/coordinates/strands, lane sizes, full oligos in PDF text and the captured compiler output.
7. After actual inspection of every final page, record `docs/协作/报告逐页核验.json` and the review narrative, then run `finalize_report.py` to bind P7/P8. This script checks the record against the current PDF; it does not perform visual inspection.

Compile from the repository root:

```bash
'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' -B scripts/B-figures/compile_report.py
```

The script records the actual Tectonic, pdftotext, pdftoppm and pdfinfo commands in `docs/协作/报告编译核验.json`, along with source, PDF, text, log and page-image hashes. Read the entire compiler output and each page raster. Build logs, extracted text and page PNGs are local ignored intermediates; their final hashes and inspection notes are retained in the compact report records. The three figure families are retained as PDF/SVG/PNG with JSON input/output manifests. The consistency audit and finalizer verify the same build snapshot. A new source or PDF invalidates the old inspection record and requires appropriate revalidation before publishing a new ready file.

The report remains a working draft because independent part-A assay parameters and author confirmation are outstanding. P7/P8 report readiness does not certify the main line's independent P9 or completion of the entire assignment.

The selected conceptual workflow is `报告/images/B/cloning-workflow-imagegen-v2.png`, generated with the built-in ImageGen tool under the user's additional illustration preference. Its original and revision prompts are `docs/协作/报告ImageGen提示词.txt` and `报告ImageGen修订提示词.txt`; the source/output review is `报告/images/B/imagegen-workflow.json`. The first output is retained for provenance and was rejected for unrequested peak traces. The second removes all peak traces and uses two full-target arrows. Pixel generation is not deterministic, so the exact outputs and hashes are retained. Symbolic icons never replace the sequence-derived map or predicted gel data. The code-generated workflow remains a reproducible detailed companion.
