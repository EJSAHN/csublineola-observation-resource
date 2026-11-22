#!/usr/bin/env python3
"""Standalone private-candidate check: export identity and raw pileup recount.
No import of sd12.py, common.py or the SD11 production engine.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path
import verify_sd11

def rows(p):
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt',newline='') as f:yield from csv.DictReader(f,delimiter='\t')

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def check(root,out):
    if sys.flags.optimize:raise RuntimeError('Assertions must be enabled; do not use Python -O')
    if out.exists():raise RuntimeError('Use a new verification output directory')
    out.mkdir(parents=True)
    hashes=0;listed=set()
    for line in (root/'SHA256SUMS.txt').read_text().splitlines():
        h,rel=line.split('  ',1);p=root/rel
        assert p.resolve().is_relative_to(root.resolve()) and rel not in listed
        assert digest(p)==h,rel
        hashes+=1;listed.add(rel)
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and '__pycache__' not in p.parts and p.suffix!='.pyc'}
    assert listed==actual,'Unlisted or missing candidate file'
    master=list(rows(root/'01_samples/all245_sample_master.tsv'))
    matrix=list(rows(root/'01_samples/matrix_samples.tsv'))
    ids=[r['fsp'] for r in matrix];assert len(ids)==len(set(ids))==142
    expected={(r['fsp'],h):r[h] for r in master for h in ['SC326-6','SC328-C','SC414-12E','BTX623','SC112-12E']}
    reactions=list(rows(root/'01_samples/host_reactions.tsv'))
    assert len(master)==245 and len(reactions)==len(expected)==1225
    assert {(r['fsp'],r['host']):r['reaction'] for r in reactions}==expected
    genotype_dir=root/'02_variants';selected=iter(rows(genotype_dir/'selected_panel_genotypes.tsv.gz'))
    count=0;missing=0;catalogue_count=0
    with (genotype_dir/'selected_panel_GT_only.vcf').open() as vf:
        while True:
            line=next(vf)
            if line.startswith('#CHROM'):break
        assert line.rstrip('\n').split('\t')[9:]==ids
        for original in rows(genotype_dir/'emitted_biallelic_catalogue.tsv.gz'):
            catalogue_count+=1
            if original['partition']!='SELECTED':continue
            r=next(selected);assert r==original,'Selected TSV differs from catalogue'
            fields=next(vf).rstrip('\n').split('\t')
            assert fields[:6]==[r['chrom'],r['pos1'],'SD10_'+r['site_index0'],r['ref'],r['alt'],r['QUAL']]
            assert fields[6:9]==['.','SRCIDX='+r['site_index0']+(';INCDS' if r['in_CDS']=='1' else ''),'GT']
            assert fields[9:]==[r[sid] for sid in ids],'Exported genotypes differ'
            missing+=fields[9:].count('.');count+=1
        assert next(selected,None) is None and not vf.read().strip()
    assert catalogue_count==493398 and count==43814
    # Reconstruct the original logical SD11 layout from the portable tables.
    # Per-sample status comes from the shipped extraction table, not an invented success flag.
    audit=root/'03_read_audit';ex={r['fsp']:r for r in rows(audit/'per_sample_extraction.tsv')}
    assert set(ex)==set(ids)
    with tempfile.TemporaryDirectory(prefix='sd12_read_check_') as tmp:
        r=Path(tmp)/'results';r.mkdir()
        for name in ['targets.tsv','samples.tsv','target_genotypes.tsv.gz','primary_strata.tsv','primary_transitions.tsv','site_support.tsv','witness_requests.tsv','witness_read_support.tsv.gz','witness_summary.tsv']:
            shutil.copyfile(audit/name,r/name)
        current=None;handle=None;writer=None;counts=collections.Counter();closed=set()
        try:
            for entry in rows(audit/'read_support_long.tsv.gz'):
                sid=entry.pop('fsp');group=entry.pop('source_group')
                assert sid in ex and group==ex[sid]['source_group']
                if sid!=current:
                    if handle:handle.close();closed.add(current)
                    assert sid not in closed,'Non-contiguous duplicate sample block'
                    d=r/'samples'/sid;d.mkdir(parents=True)
                    handle=gzip.open(d/'support.tsv.gz','wt',newline='');writer=csv.DictWriter(handle,list(entry),delimiter='\t',lineterminator='\n');writer.writeheader();current=sid
                writer.writerow(entry);counts[sid]+=1
        finally:
            if handle:handle.close()
        assert set(counts)==set(ids) and sum(counts.values())==281160
        for sid in ids:
            d=r/'samples'/sid
            assert counts[sid]==int(ex[sid]['cells'])==1980 and ex[sid]['status']=='COMPLETE'
            (d/'status.json').write_text(json.dumps({'status':ex[sid]['status']}))
            for profile in verify_sd11.PROFILES:shutil.copyfile(audit/'pileups'/sid/(profile+'.pileup.gz'),d/(profile+'.pileup.gz'))
        verify_sd11.run(Path(tmp),out/'read_recount')
    summary=dict(status='PASS',candidate_files_verified=hashes,phenotype_cells=1225,catalogue_rows=catalogue_count,selected_TSV_VCF_rows=count,selected_genotype_cells=count*len(ids),selected_missing_cells_preserved=missing,read_support_cells=281160,imports_production_engines=False,original_BAM_extraction_rerun=False)
    (out/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate_root',type=Path);p.add_argument('output_root',type=Path);a=p.parse_args();check(a.candidate_root.resolve(),a.output_root.resolve())
