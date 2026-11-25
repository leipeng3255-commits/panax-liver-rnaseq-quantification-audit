"""Run and verify one or both baseline models from the frozen filtered counts."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run_model(tag, output):
    module = ROOT / 'models' / tag
    target = output.resolve() / tag
    target.mkdir(parents=True, exist_ok=False)
    (target / 'RUN_CLAIM.json').write_text(json.dumps({
        'bundle': str(module), 'tier': 'filtered model-input tier'
    }))
    env = dict(os.environ, D03_REPRO_BUNDLE=str(module))
    subprocess.run(['Rscript', '--vanilla', 'scripts/run_model.R', '--output', str(target)],
                   env=env, cwd=module, check=True)
    subprocess.run([sys.executable, 'scripts/compare_outputs.py', str(target)],
                   cwd=module, check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--model', choices=['root', 'flower', 'both'], default='both')
    args = parser.parse_args()
    for tag in (['root', 'flower'] if args.model == 'both' else [args.model]):
        run_model(tag, args.output)
