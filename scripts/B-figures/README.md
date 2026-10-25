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

The report remains a working draft: exact siRNA/qPCR oligos, pilot-validated reagent settings and author confirmation are outstanding. P7/P8 report readiness does not certify the main line's independent P9 or completion of the entire assignment.

The historical B-v1 conceptual workflow was `报告/images/B/cloning-workflow-imagegen-v2.png`, generated with the built-in ImageGen tool under the user's additional illustration preference. Its original and revision prompts are `docs/协作/报告ImageGen提示词.txt` and `报告ImageGen修订提示词.txt`; the source/output review is `报告/images/B/imagegen-workflow.json`. The first output is retained for provenance and was rejected for unrequested peak traces. The second removes all peak traces and uses two full-target arrows. Pixel generation is not deterministic, so the exact outputs and hashes are retained. Symbolic icons never replace the sequence-derived map or predicted gel data. The code-generated workflow remains a reproducible detailed companion.


## 2026-10-25 publication revision

The current report uses ten native SnapGene Viewer 8.2.3 maps and ten main-text figure families. Explanations belong in LaTeX captions. All chromatic figure colors are defined by `publication_style.py` (ggsci NPG/nrc, fixed upstream commit); neutral text/axes remain black or gray. The final user-directed workflow selects `cloning-workflow-imagegen-v3.pdf` and `evidence-design-imagegen-v3.pdf`. Deterministic code versions and older ImageGen drafts remain available as provenance.

`build_control_panel.py` derives WT/S211A/S142A/double sequences and donors without modifying B-v1. `scripts/B-audit/check_control_panel.py` independently checks nucleotide/protein differences and donor reassembly. `build_research_panel.py` adds five C-terminal/untagged constructs; `scripts/B-audit/check_research_panel.py` independently verifies exact DNA edits, translation, tag removal, restriction changes and backbone features. Run both builders/audits before `prepare_snapgene_maps.py`, which prepares sequence-identical simplified GenBank–SnapGene display copies. Open those ten files in free Viewer, hide all enzymes, choose largest font, and export each as `<stem>-publication.svg` to `报告/质粒图-SnapGene/`. This GUI export was tested without a paid license; the CLI requires activation.

`render_publication_figures.py` preserves native geometry/vendor credit, translates only the size unit in derived SVG, composes the circular panels and draws fusion-cassette, schedule and evidence schematics. Its extra dependencies are installed under ignored `tools/figure-libs/` using `pip install --target tools/figure-libs -r scripts/B-figures/requirements-publication.txt`; the existing course environment is not changed. Exact original SVG exports remain available. Native GUI appearance is version-specific; reproduce from the retained SVG if no compatible Viewer is available.

Run the control builder/audit and prepare display copies before native export; then run the publication, accepted-map and pipeline renderers. Compile, check consistency, visually inspect every final page, and finalize the new report record. New figure manifests and mutant audits are checked by the report consistency script. Historical B-v1 P9 is not independent acceptance of the expanded cell design.

## ImageGen to editable SVG

The user subsequently requested ImageGen artwork converted to SVG, with GitHub scientific-illustration references. Bioicons and D2 informed the visual approach; no third-party icon was copied and D2 was not installed as a renderer. The built-in ImageGen tool generated the original strips. The first cloning strip's incomplete ladder motif was corrected with the same tool before tracing.

Immutable originals are in `报告/images/B/imagegen-v3/`; exact initial/edit/evidence prompts are in `docs/协作/流程图ImageGen提示词-v3.json`. Run `vectorize_imagegen_flows.py` after installing the pinned extra requirements. VTracer 0.6.15 fits spline paths; the script maps fills to the NPG palette and adds editable SVG text and connectors. Finished SVGs contain no embedded raster image. The manifest binds sources, parameters, path/text counts and outputs. Tracing approximates contours and requires visual inspection; it does not validate biological measurements.

The report explicitly requests PDF 1.7 for the vector maps. Tectonic's TeX stage still emits inclusion-version warnings before xdvipdfmx applies the setting; the final PDF version and page appearance are checked separately. Do not describe this as a warning-free compilation.
