"""Compare complete CSV tables in source order, including every identifier."""
import csv
import math
from pathlib import Path
import sys

expected = Path(__file__).resolve().parents[1] / 'reference'
run = Path(sys.argv[1])
for ref in sorted(expected.glob('*.csv')):
    matches = [p for p in run.rglob(ref.name) if 'pre71b' not in p.name]
    assert len(matches) == 1, (ref.name, matches)
    with ref.open() as f, matches[0].open() as g:
        rows = list(csv.reader(f))
        computed = list(csv.reader(g))
    assert len(rows) == len(computed), (ref.name,len(rows),len(computed))
    maxdiff = 0.
    for index, (old, new) in enumerate(zip(rows, computed)):
        assert len(old) == len(new), (ref.name,index)
        for field, x, y in zip(rows[0], old, new):
            if x == y:
                continue
            try:
                a, b = float(x), float(y)
            except ValueError:
                raise ValueError((ref.name,index,field,x,y)) from None
            assert math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-12), (ref.name,index,field,x,y)
            maxdiff = max(maxdiff,abs(a-b))
    print(f'{ref.name}: {len(rows)-1} data rows, all fields passed, largest absolute difference {maxdiff:.3g}')
