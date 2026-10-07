# Ancestry Atlas

A local, auditable workflow from a four-column genotype export to an offline HTML research report. The report combines input QC, GRCh37 assembly evidence, global reference PCA, and a separately fitted East Asian reference PCA. Sample genotypes do not participate in reference fitting or feature selection.

This repository contains source code, templates, public coordinate snapshots, licensed font subsets, and checks. Personal inputs, derived genotypes, computed reports, and numerical model files belong outside the Git publication scope.

## Windows / D drive

The Windows migration package installs a source project and a separate private research workspace under `D:\AncestryAtlas`. See [Windows setup and public GitHub publication](docs/WINDOWS.md). Run `Publish-GitHub.cmd` inside the source project after signing in with GitHub CLI.

## Run locally

Python 3.12 is the tested runtime. Install the pinned dependencies in a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
mkdir -p input
```

Place the TSV export in `input/sample.tsv`, with the exact header `gxid\tchromosome\tposition\tgenotype`. Preserve the original file. Calls are two-character A/C/G/T, I/D, or `--` encodings.

```bash
python src/qc.py input/sample.tsv --out results
python src/download_anchors.py
python src/build_evidence.py input/sample.tsv
python src/download_reference.py
python src/decompress_reference.py
python src/extract_reference.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python src/expanded_pca.py
python src/figures.py
python src/build_report.py
python -m unittest discover -s tests -p 'test_*.py' -v
python tools/check_publish.py
```

Open `report/Sample01_reference_report_v05.html` directly in a browser. Charts, fonts, controls, printing, and JSON export work offline. The report does not make background network requests.

`build_evidence.py` uses seven fixed, archived public anchors from the original study; it checks that the supplied sample contains them at the expected coordinates. It is not a general-purpose build detector. Adapt and independently verify the anchors before analyzing a different export. No automatic lift-over or ambiguous strand inference is performed.

Full public references total approximately 3.06 GB compressed and 9.02 GB after decompression. Allow at least 16 GB of free disk space and 8 GB of memory. Download sizes and SHA-256 values are pinned; decompression validates full output lengths. `extract_reference.py` saves the alignment audit and selected variant indices before reading the PGEN matrix. No large references or personal analysis results are checked into Git.

## Checks and reproduction

Source-only CI runs QC edge cases and publication guard checks. Six numerical checks are skipped when private analysis artifacts are absent. In the private research package these checks verify source calls, reference labels, frequencies, retained-window LD, saved PCA coordinates, eigenvector residuals, distances, and diagnostic accounting.

For independent genotype decoding checks, obtain PLINK 2 from its official download page, then run:

```bash
python tools/plink_crosscheck.py --binary /absolute/path/to/plink2
```

Browser checks need Playwright and Chromium. Set `ATLAS_CHROMIUM` to a compatible executable, generate the report first, then run `node tests/browser_check.cjs`. They cover 360 / 390 / 768 / 1440 px, independent chart controls, keyboard interaction, JSON export, reduced motion, print behavior, no external requests, and no runtime errors.

The bundled CJK font subsets are renamed Atlas Serif SC and covered by the included SIL OFL. `src/subset_fonts.py` can regenerate them from Fontsource Noto Serif SC 5.3.0 assets placed in `references/fontsource/package`. Keep the OFL notice when redistributing the subsets.

## Interpretation

PCA is a relative reference-space comparison. Distances and internal holdout agreement are not ancestry proportions, identity labels, or personal accuracy probabilities. Nearby population centers can change order under solver or training-reference changes. The exported report keeps limitations next to the relevant figures.

No ADMIXTURE, NNLS ancestry proportions, ancient-DNA affinity, Y haplogroup, mitochondrial haplogroup, or migration chronology is inferred by this version.

See [methods](docs/METHODS.md), [design](docs/DESIGN.md), and [sources](docs/SOURCES.md). Original code is MIT-licensed; reference data and fonts retain their respective terms.
