#!/usr/bin/env python3
"""Compare every tabular field and summary value with frozen round40 references.

No tolerance is used to declare numerical equality. Nonzero numerical changes
require investigation; 1e-10 bounds below are reported only to characterize them.
"""
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
TEXT = {
    'gene_results_all.tsv': {'gene_id','symbol','gene_type','unique_symbol_in_domain'},
    'voom_expression.tsv': {'gene_id'},
    'voom_weights.tsv': {'gene_id'},
    'pathway_coverage.tsv': {'gene_set','eligible'},
    'pathway_members.tsv': {'gene_set','gene_id','symbol'},
    'pathway_CAMERA_all.tsv': {'scenario','correlation_setting','gene_set','Direction'},
    'gene_leave_one_all.tsv': {'omitted','gene_id'},
    'pathway_CAMERA_leave_one.tsv': {'scenario','correlation_setting','gene_set','Direction'},
    'gene_leave_one_summary.tsv': {'gene_id'},
    'design_matrix.tsv': {'run'},
}

def rows(blob):
    return list(csv.reader(io.StringIO(blob.decode()),delimiter='\t'))

def numeric(a,b):
    if a == 'NA' or b == 'NA': return (a==b,0.0,0.0,a==b)
    x,y = float(a),float(b)
    if not math.isfinite(x) or not math.isfinite(y):
        return (x==y and not math.isnan(x),0.0,0.0,x==y and not math.isnan(x))
    delta = abs(x-y)
    relative = delta/max(abs(x),abs(y),1e-300)
    return x==y,delta,relative,delta <= 1e-10+1e-10*max(abs(x),abs(y))

def compare_table(name, actual):
    expected_bytes = gzip.decompress((ROOT/'reference'/(name+'.gz')).read_bytes())
    actual_bytes = actual.read_bytes()
    expected, observed = rows(expected_bytes),rows(actual_bytes)
    report = {'file':name,'expected_rows':len(expected)-1,'actual_rows':len(observed)-1,
        'header_and_column_order_equal':expected[0]==observed[0],
        'bytes_identical':expected_bytes==actual_bytes,
        'expected_uncompressed_sha256':hashlib.sha256(expected_bytes).hexdigest(),
        'actual_sha256':hashlib.sha256(actual_bytes).hexdigest(),
        'numeric_cells_checked':0,'text_ID_boolean_cells_checked':0,
        'numeric_different_cells':0,'text_ID_boolean_different_cells':0,
        'max_abs_difference':0.0,'max_relative_difference':0.0,
        'outside_descriptive_1e_10_bounds':0,'first_differences':[]}
    if expected[0]!=observed[0] or len(expected)!=len(observed):
        report['status']='FAIL_SCHEMA_OR_ROW_COUNT';return report
    for i,(erow,arow) in enumerate(zip(expected[1:],observed[1:]),start=2):
        if len(erow)!=len(arow):
            report['status']='FAIL_ROW_WIDTH';report['bad_line']=i;return report
        for column,e,a in zip(expected[0],erow,arow):
            if column in TEXT[name]:
                report['text_ID_boolean_cells_checked']+=1;equal=e==a
                report['text_ID_boolean_different_cells']+=int(not equal)
            else:
                report['numeric_cells_checked']+=1
                try: equal,absolute,relative,within=numeric(e,a)
                except ValueError:
                    equal,absolute,relative,within=False,0.0,0.0,False
                report['numeric_different_cells']+=int(not equal)
                report['max_abs_difference']=max(report['max_abs_difference'],absolute)
                report['max_relative_difference']=max(report['max_relative_difference'],relative)
                report['outside_descriptive_1e_10_bounds']+=int(not within)
            if not equal and len(report['first_differences'])<20:
                report['first_differences'].append({'line':i,'column':column,'expected':e,'actual':a})
    report['status']='PASS_ALL_FIELDS_EQUAL' if not (report['numeric_different_cells']+report['text_ID_boolean_different_cells']) else 'DIFFERENCE_REQUIRES_REVIEW'
    return report

def compare_json(expected,actual,path='$',details=None):
    if details is None: details=[]
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or list(expected)!=list(actual):
            details.append({'path':path,'equal':False,'kind':'object_keys_or_order'});return details
        for key in expected: compare_json(expected[key],actual[key],path+'.'+key,details)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(expected)!=len(actual):
            details.append({'path':path,'equal':False,'kind':'list_length'});return details
        for i,(e,a) in enumerate(zip(expected,actual)):compare_json(e,a,f'{path}[{i}]',details)
    else:
        equal = (type(expected) is type(actual) or (type(expected) in (int,float) and type(actual) in (int,float))) and expected==actual
        details.append({'path':path,'expected':expected,'actual':actual,'equal':equal})
    return details

def compare_y_samples(output):
    qc=list(csv.DictReader(io.StringIO(gzip.decompress((ROOT/'reference/sample_QC.tsv.gz').read_bytes()).decode()),delimiter='\t'))
    current=list(csv.DictReader((output/'reconstructed_y_samples.tsv').open(),delimiter='\t'))
    values=[]
    if len(qc)!=len(current):return {'status':'FAIL_ROW_COUNT'}
    for e,a in zip(qc,current):
        for actual_key,expected_key in [('run','run'),('group','group'),('lib.size','filtered_total'),('norm.factors','TMM')]:
            eq=e[expected_key]==a[actual_key] if actual_key in ('run','group') else numeric(e[expected_key],a[actual_key])[0]
            values.append({'run':e['run'],'field':actual_key,'expected':e[expected_key],'actual':a[actual_key],'equal':eq})
    return {'status':'PASS_ALL_FIELDS_EQUAL' if all(v['equal'] for v in values) else 'DIFFERENCE_REQUIRES_REVIEW','values':values}

def main():
    if len(sys.argv)!=2:raise SystemExit('Usage: compare_outputs.py completed_output_directory')
    output=Path(sys.argv[1]).resolve()
    tables=[compare_table(name,output/name) for name in TEXT]
    expected_summary=json.loads((ROOT/'reference/RESULTS_summary.json').read_text())
    actual_summary=json.loads((output/'RESULTS_summary.json').read_text())
    summary=compare_json(expected_summary,actual_summary)
    normalization=compare_y_samples(output)
    passed=all(x['status']=='PASS_ALL_FIELDS_EQUAL' for x in tables) and all(x['equal'] for x in summary) and normalization['status']=='PASS_ALL_FIELDS_EQUAL'
    report={'status':'PASS_EXACT_NUMERIC_AND_ID_ORDER_REPRODUCTION' if passed else 'DIFFERENCE_REQUIRES_REVIEW',
        'tier':'filtered model-input tier','acceptance':'Every compared numeric value and every ID/text/boolean field must be equal; row and column order must match. No tolerance-based pass.',
        'reference_values_are_not_used_by_statistical_core':True,
        'tables':tables,'summary_all_values':summary,'reconstructed_normalization_and_groups':normalization,
        'numeric_table_cells_checked':sum(t['numeric_cells_checked'] for t in tables),
        'text_ID_boolean_table_cells_checked':sum(t['text_ID_boolean_cells_checked'] for t in tables),
        'summary_leaf_values_checked':len(summary),'new_model_objects_sha256':hashlib.sha256((output/'model_objects.rds').read_bytes()).hexdigest(),
        'model_object_and_PDF_bytes_not_used_as_acceptance_criterion':True,
        'does_not_revalidate_original_filtering_tximport_FASTQ_or_experimental_units':True}
    (output/'COMPARISON_QA.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('tables','summary_all_values','reconstructed_normalization_and_groups')},indent=2))
    return 0 if passed else 1

if __name__=='__main__':raise SystemExit(main())
