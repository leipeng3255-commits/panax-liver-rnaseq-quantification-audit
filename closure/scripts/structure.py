from pathlib import Path
import csv, gzip, hashlib, json, re, statistics
from collections import defaultdict
import sys
o=Path(sys.argv[1])
p=Path("closure/inputs/Nat8l_M37.gtf.gz")
sha="75c3cf3c940b3b305ac1bf1ee48aa50c8c1535134e1f92117b8e9f864dc8bbe4"  # SHA-256 of the complete source M37 GTF, before locus selection.
gene='ENSMUSG00000048142.10'
data=defaultdict(list)
with gzip.open(p,'rt') as f:
 for line in f:
  if gene not in line:continue
  chrom,source,typ,start,end,score,strand,phase,att=line.rstrip().split('\t');attrs=dict(re.findall(r'(\w+) "([^"]+)";',att))
  if typ in ('transcript','exon','CDS'):data[attrs['transcript_name']].append((typ,chrom,int(start),int(end),strand,attrs))
assert sorted(data)==['Nat8l-201','Nat8l-202']
def save(name,rows):
 dest=o/name
 with dest.open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
rows=[]
exons={t:sorted((e[2],e[3]) for e in entries if e[0]=='exon') for t,entries in data.items()}
for t in sorted(data):
 g,=[i for i in data[t] if i[0]=='transcript'];ex=exons[t];assert all(x[1]<y[0] for x,y in zip(ex,ex[1:]))
 overlap=sum(max(0,min(b,d1)-max(a,c)+1) for a,b in ex for c,d1 in exons['Nat8l-201'])
 rows.append(dict(transcript=t,transcript_id=g[5]['transcript_id'],gene_id=gene,full_length_nt=sum(b-a+1 for a,b in ex),biotype=g[5]['transcript_type'],chromosome=g[1],strand=g[4],genomic_start=g[2],genomic_end=g[3],exon_count=len(ex),exon_intervals='|'.join(f'{a}-{b}' for a,b in ex),exonic_overlap_with_201_nt=overlap,CDS_segments=sum(i[0]=='CDS' for i in data[t]),evidence=f'GRCm39 GENCODE M37 one-based closed coordinates; SHA256 {sha}'))
assert [r['full_length_nt'] for r in rows]==[6529,205]
save('nat8l_annotation_structure.csv',rows);print(rows)
with gzip.open("closure/inputs/flower_samples.csv.gz","rt") as f:
 labels={x['run']:x['isolate'] for x in csv.DictReader(f)}
out=[]
for x in csv.DictReader((o/'flower_split_half_reliability.csv').open()):
 if x['contrast']!='treatment' or x['row_type']!='pair':continue
 z={k:x[k] for k in ('pair','statistic','rho')};z['rho']=float(z['rho'])
 for k in ('half1_hi','half1_lo','half2_hi','half2_lo'):z[k]='|'.join(f'{labels[i]} ({i})' for i in x[k].split('|'))
 hs=next(k for k in ('half1_hi','half2_hi') if 'PNF20_1' in z[k]);z['PNF20_partner']=next(labels[i] for i in x[hs].split('|') if labels[i]!='PNF20_1')
 out.append(z)
for z in out:z['rank_of_18']=1+sum(t['rho']<z['rho'] for t in out if t['statistic']==z['statistic'])
assert len(out)==36
save('round71b/structure/split_treatment_labels.csv',out)
for z in out:
 if z['rank_of_18']==1:
  print(z);print('same PNF partner range',min(y['rho'] for y in out if y['statistic']==z['statistic'] and y['PNF20_partner']==z['PNF20_partner']),max(y['rho'] for y in out if y['statistic']==z['statistic'] and y['PNF20_partner']==z['PNF20_partner']))
