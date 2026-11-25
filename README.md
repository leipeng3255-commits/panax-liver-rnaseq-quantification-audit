# Reproducing the mouse liver reanalysis

Run these commands from this directory with Rscript and Python 3 on PATH. Use the versions in ENVIRONMENT_LOCK.json. Set R_LIBS if your installed R packages are outside the default search path. The analysis runners use the Python standard library. Select a fresh writable output directory.

```sh
python3 run_models.py ../reanalysis_run
python3 run_closure.py ../reanalysis_run
```

The first command refits the root and flower treatment models. Results, comparison reports and R session information are written into the root and flower subdirectories of the selected output. Each model retains the original gene domain, sample order, contrasts, random seed and statistical core. The flower model also uses a frozen root reference for its shared-gene diagnostic. Its own reference outputs are used only during verification.

The second command runs the methodological sensitivity analyses, including all four count-scale procedures, disease models and CAMERA, cross-study correlations, 18 balanced partitions per flower contrast, and the Nat8l locus and transcript abundance shares. The treatment offset comparisons report edgeR QL and TREAT. For the flower disease offset contrast, voom uses library-size and gene-by-sample length terms.

The full-precision flower input contains three gene-level matrices from tximport 1.40.0 across all 12 libraries: TPM, unscaled estimated counts and sample-specific average effective lengths. The count-scale analysis uses this object to preserve transcript-import precision. The disease and calibration calculations use their original decimal exports. Root disease inputs contain the six relevant libraries. The Nat8l inputs include the GENCODE M37 locus and the 24 original transcript–library quantification records. Coordinates are one-based and inclusive. Raw-read processing and read-level alignment diagnostics are outside this reproduction tier.

The closure verifier compares all nine CSVs by row, column, identifier and numeric field (relative tolerance 1e−8, absolute tolerance 1e−12). The eight methodological reference tables matched the recalculated values exactly. Across the 24-row transcript-share table the largest difference was 2.22e−16. Both treatment models passed their original comparisons without numerical tolerance.

Main-figure source scripts and plotting inputs are under main_figures/source_code. Follow main_figures/source_code/README.md to generate vector graphics and 300 dpi previews. The Hallmark gene sets and attribution are under each model’s resources directory. ENVIRONMENT_LOCK.json records the portable inputs and upstream source identities. The package-level manifest lists the size and SHA-256 digest of each payload file.

Genome-wide differential expression results for the root and flower studies are provided in archived_full_tables/. This directory also contains the full leave-one-library-out fits and individual subsampling draws.
