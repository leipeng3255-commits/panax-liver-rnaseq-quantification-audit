# Main figure reproduction

Run the following commands from this directory with the Python plotting environment specified in ENVIRONMENT_LOCK.json. All inputs are included. The scripts write each vector PDF, editable SVG and 300 dpi PNG preview to the adjacent exports directory.

```sh
python draw_figure1.py
python draw_figure2.py
python draw_figure3.py
python draw_figure4.py
```

Run `python3 draw_figure1.py --svg-only` for the editable diagram, or `python3 draw_figure1.py` to also export a vector PDF and a 300 dpi PNG. Use Python 3.9 or newer. PDF and PNG export requires Inkscape and DejaVu Sans regular and bold. Use `--inkscape` for a nonstandard executable location and `--font-dir` to supply the TTF directory directly. The default font directory comes from Matplotlib. Use `--out DIRECTORY` to choose an output directory.

The layout template was reconstructed from the approved editable diagram. Library sizes, gene-domain sizes, partition and transcript counts, count schemes and Hallmark coverage are read from the accompanying analysis tables. The diagram icons are schematic. With these inputs, the SVG matches the submission SVG byte for byte. PDF metadata and font encoding depend on the exporter.

The design_counts table uses the root and flower treatment-model gene results and distinct transcript names in Supplementary Table S8. The additional disease CAMERA input comes from closure/reference/disease_positive_control_camera.csv. Data for the signed significance panel use the root A and flower A disease models under fixed zero correlation.
