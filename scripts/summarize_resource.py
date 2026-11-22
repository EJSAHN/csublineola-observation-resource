#!/usr/bin/env python3
"""Recompute manuscript source tables from an intact, verified resource export.
All calculations are descriptive; no imputation, association or new variant calls.
"""
import collections
import csv
import gzip
import hashlib
import json
import statistics
from pathlib import Path

HOSTS=('SC326-6','SC328-C','SC414-12E','BTX623','SC112-12E')
PARTS=('SELECTED','OUTSIDE_MASK_ONLY','LOW_CALL_COUNT_ONLY','BOTH',
       'MONOMORPHIC_REF','MONOMORPHIC_ALT','NO_RETAINED')
PROFILES=('MQ20_BQ20_noBAQ','MQ40_BQ20_noBAQ','MQ40_BQ20_redoBAQ')

def rows(path):
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',newline='') as f: yield from csv.DictReader(f,delimiter='\t')

def table(path,data,fields=None):
    data=list(data)
    fields=fields or list(data[0])
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(data)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def classify(nref,nalt,in_mask,required=128):
    if nref and nalt:
        return ('SELECTED' if in_mask else 'OUTSIDE_MASK_ONLY') if nref+nalt>=required else ('LOW_CALL_COUNT_ONLY' if in_mask else 'BOTH')
    return 'MONOMORPHIC_REF' if nref else 'MONOMORPHIC_ALT' if nalt else 'NO_RETAINED'

def reaction_counts(data,ids):
    counts=collections.Counter()
    for r in data:
        if r['fsp'] in ids:
            if r['reaction'] not in ('R','S','V'):raise ValueError('Invalid reaction')
            counts[r['host'],r['reaction']]+=1
    return [dict(host=h,n=len(ids),R=counts[h,'R'],S=counts[h,'S'],V=counts[h,'V']) for h in HOSTS]

def summarize(data,out):
    if out.exists():raise ValueError('Use a new summary directory')
    out.mkdir(parents=True)
    sm=list(rows(data/'01_samples/sample_manifest.tsv'))
    ids={r['fsp'] for r in sm};assert len(ids)==len(sm)==245
    primary={r['fsp'] for r in sm if r['in_primary_genotype_panel']=='1'}
    sequenced={r['fsp'] for r in sm if r['recorded_BAM_available']=='1'}
    assert primary<=sequenced<=ids
    cohort=[dict(cohort=k,isolates=len(v)) for k,v in [('all_phenotyped',ids),('sequenced',sequenced),('primary_genotype',primary),('sparse_sequenced',sequenced-primary),('phenotype_only',ids-sequenced)]]
    table(out/'cohort_counts.tsv',cohort)
    reactions=list(rows(data/'01_samples/host_reactions.tsv'))
    assert len(reactions)==len({(r['fsp'],r['host']) for r in reactions})==1225
    host=[]
    for name,group in [('all245',ids),('primary142',primary)]:
        host.extend(dict(cohort=name,**r) for r in reaction_counts(reactions,group))
    table(out/'host_reaction_counts.tsv',host)
    master=list(rows(data/'01_samples/all245_sample_master.tsv'))
    observed=[]
    for r in master:
        if r['fsp'] not in sequenced:continue
        observed.append(dict(fsp=r['fsp'],source_group=r['source_group'],primary=int(r['fsp'] in primary),mean_aligned_depth=float(r['MQ20_BQ20_genome_mean_depth']),fraction_reference_DP5=float(r['MQ20_BQ20_genome_fraction_DP5']),fraction_CDS_DP5=float(r['MQ20_BQ20_CDS_fraction_DP5'])))
    table(out/'observation_by_isolate.tsv',observed)
    stats=[]
    for group in sorted({r['source_group'] for r in observed}):
        g=[r for r in observed if r['source_group']==group]
        for metric in ['mean_aligned_depth','fraction_reference_DP5','fraction_CDS_DP5']:
            v=[r[metric] for r in g]
            stats.append(dict(source_group=group,isolates=len(g),metric=metric,minimum=min(v),median=statistics.median(v),maximum=max(v)))
    table(out/'observation_by_source_group.tsv',stats)
    count=collections.Counter();missing=0;total=0
    for r in rows(data/'02_variants/emitted_biallelic_catalogue.tsv.gz'):
        a,b=int(r['n_REF']),int(r['n_ALT'])
        values=[r[s] for s in primary]
        assert values.count('0')==a and values.count('1')==b
        assert set(values)<= {'0','1','.'}
        assert classify(a,b,r['mask_selected']=='1')==r['partition']
        count['REFERENCE',r['partition']]+=1
        if r['in_CDS']=='1':count['CDS',r['partition']]+=1
        if r['partition']=='SELECTED':missing+=values.count('.')
        total+=1
    parts=[dict(panel=p,partition=k,sites=count[p,k]) for p in ('REFERENCE','CDS') for k in PARTS]
    table(out/'catalogue_partitions.tsv',parts)
    assert total==493398 and count['REFERENCE','SELECTED']==43814 and count['CDS','SELECTED']==9580
    audit=collections.defaultdict(collections.Counter)
    for r in rows(data/'03_read_audit/read_support_long.tsv.gz'):
        if r['primary_sample']!='1' or r['baseline']=='.':continue
        key=(r['partition'],r['profile'],'REF' if r['baseline']=='0' else 'ALT')
        c=audit[key];c['baseline_retained']+=1;c[r['baseline_relation']]+=1
    evidence=[]
    for part in PARTS[:4]:
        for profile in PROFILES:
            for allele in ('REF','ALT'):
                c=audit[part,profile,allele];n=c['baseline_retained']
                assert n==c['SAME_ALLELE']+c['OPPOSITE_ALLELE']+c['UNRESOLVED_READ_SUPPORT']
                evidence.append(dict(partition=part,profile=profile,baseline_allele=allele,baseline_retained=n,same=c['SAME_ALLELE'],unresolved=c['UNRESOLVED_READ_SUPPORT'],opposite=c['OPPOSITE_ALLELE'],same_fraction=c['SAME_ALLELE']/n if n else ''))
    table(out/'read_support_by_allele.tsv',evidence)
    existing={(r['partition'],r['profile'],r['baseline_allele']):r for r in rows(data/'04_validation/independent_recount/primary_allele_specific.tsv')}
    for r in evidence:
        old=existing[r['partition'],r['profile'],r['baseline_allele']]
        for k in ('baseline_retained','same','unresolved','opposite'):assert r[k]==int(old[k])
    sites=collections.defaultdict(collections.Counter)
    for r in rows(data/'03_read_audit/site_support.tsv'):
        if r['primary_sample']=='1':
            c=sites[r['partition'],r['profile']];c['sampled_sites']+=1;c['both_alleles_supported']+=int(r['both_alleles_supported'])
    table(out/'read_support_by_site.tsv',[dict(partition=p,profile=q,**sites[p,q]) for p in PARTS[:4] for q in PROFILES])
    report=dict(status='PASS',cohorts={r['cohort']:r['isolates'] for r in cohort},catalogue_records=total,selected_sites=count['REFERENCE','SELECTED'],selected_CDS_sites=count['CDS','SELECTED'],selected_missing_genotypes=missing,source_release_manifest_sha256=sha(data/'SHA256SUMS.txt'),phenotype_model_fitted=False,new_variant_calls=0)
    (out/'metrics.json').write_text(json.dumps(report,indent=2)+'\n')
    used=['01_samples/sample_manifest.tsv','01_samples/host_reactions.tsv','01_samples/all245_sample_master.tsv','02_variants/emitted_biallelic_catalogue.tsv.gz','03_read_audit/read_support_long.tsv.gz','03_read_audit/site_support.tsv','04_validation/independent_recount/primary_allele_specific.tsv']
    table(out/'input_provenance.tsv',[dict(file=f,sha256=sha(data/f)) for f in used])
    return report
