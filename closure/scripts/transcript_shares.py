"""Calculate within-gene abundance and estimated fragment fractions from Salmon records."""
import csv
from collections import defaultdict
from pathlib import Path
import sys

inputs = Path(__file__).resolve().parents[1] / 'inputs'
with (inputs / 'Nat8l_transcript_records.csv').open() as file:
    reader = csv.DictReader(file)
    fields = reader.fieldnames + ["gene_TPM","unscaled_estimated_fragments","fraction_gene_TPM","fraction_gene_estimated_fragments"]
    records = list(reader)

by_run = defaultdict(list)
for record in records:
    by_run[record['run']].append(record)

assert len(by_run) == 12
for rows in by_run.values():
    assert {x['transcript_name'] for x in rows} == {'Nat8l-201','Nat8l-202'}
    for source, total, share in [('TPM','gene_TPM','fraction_gene_TPM'),
                                 ('NumReads_estimated_fragments','unscaled_estimated_fragments','fraction_gene_estimated_fragments')]:
        denominator = sum(float(x[source]) for x in rows)
        for row in rows:
            row[total] = denominator
            row[share] = float(row[source]) / denominator if denominator else ''
with (Path(sys.argv[1]) / 'Nat8l_transcript_shares.csv').open('w') as file:
    writer = csv.DictWriter(file, fieldnames=fields)
    writer.writeheader()
    writer.writerows(records)
