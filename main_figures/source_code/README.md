# Main figure reproduction

Run the following commands from this directory with the Python plotting environment specified in ENVIRONMENT_LOCK.json. All inputs are included. The scripts write each vector PDF, editable SVG and 300 dpi PNG preview to the adjacent exports directory.

```sh
python draw_figure1.py
python draw_figure2.py
python draw_figure3.py
python draw_figure4.py
```

The design_counts table uses the root and flower treatment-model gene results and distinct transcript names in Supplementary Table S18. The additional disease CAMERA input comes from closure/reference/disease_positive_control_camera.csv. Data for the signed significance panel use the root A and flower A disease models under fixed zero correlation.
