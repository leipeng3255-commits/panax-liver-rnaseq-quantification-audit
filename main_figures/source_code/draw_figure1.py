"""Render the submission design diagram from its editable vector layout and data."""
import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from string import Template
from xml.sax.saxutils import escape

HERE = Path(__file__).resolve().parent


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def diagram_labels(inputs):
    fig12 = inputs / "fig12"
    groups = Counter(row["group"] for row in read_csv(fig12 / "library_groups.csv"))
    design = {row["metric"]: int(row["value"]) for row in read_csv(fig12 / "design_counts.csv")}
    domains = {row["statistic"]: int(row["n_genes"])
               for row in read_csv(fig12 / "cross_context_snr_calibration.csv")
               if row["stratum"] == "all" and row["contrast"] == "treatment"}
    splits = read_csv(fig12 / "flower_split_half_reliability.csv")
    pairs = {row["pair"] for row in splits if row["row_type"] == "pair" and row["contrast"] == "treatment"}
    schemes = {row["scheme"] for row in read_csv(inputs / "tximport_scale_sensitivity_13genes.csv")}
    with gzip.open(inputs / "flower_treatment_CAMERA.tsv.gz", "rt") as stream:
        sets = {row["gene_set"] for row in csv.DictReader(stream, delimiter="\t")}
    return {**{group: groups[group] for group in ("Control", "Model", "PNS", "NC", "DM", "PNF")},
            "root_domain": format(design["root_treatment_genes"], ","),
            "flower_domain": format(design["flower_treatment_genes"], ","),
            "transcripts": design["flower_transcripts"],
            "tpm_domain": format(domains["log2(1+TPM) mean difference"], ","),
            "model_domain": format(domains["edgeR_logFC"], ","),
            "splits": len(pairs), "schemes": len(schemes), "scheme_names": ", ".join(sorted(schemes)),
            "hallmarks": len(sets)}


def render(svg, output, inkscape, font_dir):
    if not all((font_dir / name).is_file() for name in ("DejaVuSans.ttf", "DejaVuSans-Bold.ttf", "DejaVuSans-Oblique.ttf")):
        raise ValueError("Provide the DejaVu Sans regular, bold and oblique TTF files using --font-dir.")
    # The approved SVG uses CSS shorthand. Specify the PDF typeface explicitly
    # during export to avoid the system's default sans-serif substitution.
    def expand_font(match):
        weight = "bold" if match[1] else "normal"
        return f"font-family:DejaVu Sans; font-weight:{weight}; font-size:{match[2]}px"

    styled = re.sub(r"font: (700 )?([\d.]+)px 'sans-serif'", expand_font, svg)
    with tempfile.TemporaryDirectory() as temp:
        work = Path(temp)
        source = work / "figure.svg"
        source.write_text(styled, encoding="utf-8")
        conf = work / "fonts.conf"
        conf.write_text(f"<fontconfig><dir>{escape(str(font_dir.resolve()))}</dir>"
                        f"<cachedir>{escape(str(work))}</cachedir></fontconfig>", encoding="utf-8")
        env = dict(os.environ, FONTCONFIG_FILE=str(conf))
        subprocess.run([inkscape, str(source), f"--export-filename={output / 'Figure1.pdf'}"],
                       env=env, check=True)
        previews = output / "preview_png"
        previews.mkdir(exist_ok=True)
        subprocess.run([inkscape, str(source), "--export-dpi=300",
                        f"--export-filename={previews / 'Figure1.png'}"], env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE.parent / "exports")
    parser.add_argument("--svg-only", action="store_true")
    parser.add_argument("--inkscape", default=shutil.which("inkscape"))
    parser.add_argument("--font-dir", type=Path, help="Directory containing DejaVu Sans TTF files")
    args = parser.parse_args()
    layout = Template((HERE / "figure1_layout.svg.in").read_text(encoding="utf-8"))
    svg = layout.substitute(diagram_labels(HERE / "inputs"))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "Figure1.svg").write_bytes(svg.encode("utf-8"))
    if args.svg_only:
        return
    if not args.inkscape:
        parser.error("Install Inkscape or supply its executable with --inkscape.")
    if args.font_dir is None:
        import matplotlib
        args.font_dir = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    render(svg, args.out.resolve(), args.inkscape, args.font_dir)


if __name__ == "__main__":
    main()
