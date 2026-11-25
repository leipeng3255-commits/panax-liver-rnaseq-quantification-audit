"""Recompute scale, disease, split-half and annotation results."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys

BASE = Path(__file__).resolve().parent
SCRIPTS = BASE / 'closure/scripts'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='Directory containing regenerated root/ and flower/ models')
    args = parser.parse_args()
    output = args.output.resolve()
    for name in ['root', 'flower']:
        if not (output / name / 'gene_results_all.tsv').is_file():
            parser.error('Run run_models.py first')
    work = output / 'round71b'
    for step in ['task1', 'task2', 'root_qa', 'structure']:
        (work / step).mkdir(parents=True, exist_ok=False)
    def r(name):
        subprocess.run(['Rscript','--vanilla',str(SCRIPTS / name), str(output)], cwd=BASE, check=True)
    r('scale.R')
    r('calibrate.R')
    for name, target in [('cross_context_snr_calibration','cross'),('disease_positive_control_camera','camera')]:
        (output / (name+'.csv')).rename(work / 'task2' / ('pre71b_'+target+'.csv'))
    r('disease.R')
    r('disease_split.R')
    path = output / 'flower_split_half_reliability.csv'
    with (work / 'root_qa/flower_disease_A_split_half.csv').open() as infile, path.open('a') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writerows(row for row in reader if row['statistic'] == 'voom_t')
    subprocess.run([sys.executable,str(SCRIPTS/'structure.py'),str(output)],cwd=BASE,check=True)
    subprocess.run([sys.executable,str(SCRIPTS/'transcript_shares.py'),str(output)],cwd=BASE,check=True)
    subprocess.run([sys.executable,str(SCRIPTS/'verify.py'),str(output)],cwd=BASE,check=True)


if __name__ == '__main__':
    main()
