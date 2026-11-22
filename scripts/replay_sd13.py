#!/usr/bin/env python3
"""Verify a completed SD13 REVIEW export and recount archived pileups read-only.

The original byte-identical independent parser writes a report. It therefore runs
on a temporary copy, not on the frozen supplied input. No production engine import.
"""
import argparse, csv, hashlib, json, os, shutil, subprocess, sys, tempfile
from pathlib import Path
EXCLUDED = ['FSP348','FSP350','FSP351','FSP353','FSP357','FSP360','FSP363','FSP369']
PROFILES = ['MQ20_BQ20_noBAQ','MQ40_BQ20_noBAQ','MQ40_BQ20_redoBAQ']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def check_manifest(root, name, exclude_subtrees=()):
    root=Path(root).resolve();manifest=root/name
    if not manifest.is_file():raise FileNotFoundError('Required manifest not found: '+str(manifest))
    seen=set()
    for line in manifest.read_text().splitlines():
        expected,rel=line.split('  ',1);parts=Path(rel).parts;p=(root/rel).resolve()
        if not parts or root not in p.parents or rel in seen:raise ValueError('Unsafe or repeated manifest member: '+rel)
        if '__pycache__' in parts or p.suffix in ('.pyc','.pyo'):raise ValueError('Runtime cache in source manifest: '+rel)
        if not p.is_file() or sha(p)!=expected:raise ValueError('Checksum mismatch: '+rel)
        seen.add(rel)
    if not seen:raise ValueError('Empty manifest')
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p!=manifest and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo') and not any(p.relative_to(root).parts[0]==d for d in exclude_subtrees)}
    if actual!=seen:raise ValueError('Missing or unlisted files in '+str(root))
    return len(seen)

def rows(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def table(p,data):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)

def replay(run, out):
    if not __debug__:raise RuntimeError('Run without -O: independent parser uses assertions')
    run=Path(run).resolve();out=Path(out).resolve()
    if not run.is_dir():raise FileNotFoundError('SD13 run directory not found: '+str(run))
    if out==run or run in out.parents:raise ValueError('Do not place replay output within frozen SD13 input')
    if out.exists():raise FileExistsError('Choose a new SD13 output directory: '+str(out))
    counts={'review':check_manifest(run,'REVIEW_SHA256SUMS.txt',exclude_subtrees=('results',)),'code':check_manifest(run/'code','SHA256SUMS.txt'),'results':check_manifest(run/'results','SHA256SUMS.txt')}
    status=json.loads((run/'runner_status.json').read_text())
    if status['exit_code']!=0 or status['status'] not in ('COMPLETE_READ_SUPPORT_MASK_PENDING','COMPLETE_READ_SUPPORT_REVIEW_REQUIRED'):
        raise ValueError('Input run did not complete the read audit: '+str(status))
    # Binding the parser to the executed code prevents silent substitutions.
    parser=Path(__file__).with_name('verify_sd13.py')
    if sha(parser)!=sha(run/'code/verify_support.py'):raise ValueError('Independent parser differs from executed SD13 parser')
    result=run/'results';samples=rows(result/'samples.tsv');targets=rows(result/'targets.tsv')
    ids=[r['fsp'] for r in samples]
    if len(ids)!=142 or len(set(ids))!=142 or not set(EXCLUDED)<=set(ids) or len(targets)!=425:raise ValueError('Unexpected fixed target/sample dimensions')
    if sum(r['polymorphic_in_focus']=='1' for r in targets)!=421:raise ValueError('Expected 421 focus-polymorphic targets')
    out.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix='sd13_recount_') as tmp:
        temp=Path(tmp)
        for n in ('targets.tsv','samples.tsv','target_genotypes.tsv.gz','support_by_group_and_baseline.tsv','site_support_by_group.tsv'):
            shutil.copy2(result/n,temp/n)
        for sid in ids:
            dest=temp/'samples'/sid;dest.mkdir(parents=True)
            for name in ['support.tsv.gz']+[p+'.pileup.gz' for p in PROFILES]:shutil.copy2(result/'samples'/sid/name,dest/name)
        cmd=[sys.executable,str(parser),str(temp),'--exclude',*EXCLUDED]
        with (out/'independent_parser.log').open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        report=json.loads((temp/'independent_verification.json').read_text())
        if report['status']!='PASS' or report['support_cells']!=181050:raise ValueError('Unexpected independent-parser result')
        save(out/'independent_verification.json',report)
    # Rich manuscript summaries read only independently verified support labels.
    import summarize_targeted_support
    manuscript_summary = summarize_targeted_support.summarize(result,out/'paper_summary')
    for name in ('support_by_group_and_baseline.tsv','site_support_by_group.tsv'):
        shutil.copy2(result/name,out/name)
    mask=json.loads((result/'common_mask_scope.json').read_text());save(out/'common_mask_scope.json',mask)
    # Recheck integrity after the external verifier to demonstrate read-only inputs.
    post_counts={'review':check_manifest(run,'REVIEW_SHA256SUMS.txt',exclude_subtrees=('results',)),'code':check_manifest(run/'code','SHA256SUMS.txt'),'results':check_manifest(run/'results','SHA256SUMS.txt')}
    if post_counts!=counts:raise ValueError('Input changed during replay')
    summary={'status':'PASS','manifest_files':counts,'independent_verification':report,'sample_count':142,'target_sites':425,'polymorphic_in_focus':421,'manuscript_summary':manuscript_summary,'common_mask_status':mask['status'],'mask_recomputed_locally':False,'source_review_manifest_sha256':sha(run/'REVIEW_SHA256SUMS.txt'),'scope':'Archived-pileup replay and recorded mask status; no original BAM extraction, no new variant calls, no new mask calculation.'}
    save(out/'summary.json',summary);return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();print(json.dumps(replay(a.run,a.out),indent=2))
