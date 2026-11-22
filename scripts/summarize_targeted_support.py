#!/usr/bin/env python3
"""Recount SD13 exported support by within-focus polymorphism, without production imports.

This summary reads exported support labels. Run the independent raw-pileup verifier
first; this program alone does not validate the underlying BAM extraction.
"""
import argparse,collections,csv,gzip,hashlib,json
from pathlib import Path
PROFILES=('MQ20_BQ20_noBAQ','MQ40_BQ20_noBAQ','MQ40_BQ20_redoBAQ')
EXCLUDED={'FSP348','FSP350','FSP351','FSP353','FSP357','FSP360','FSP363','FSP369'}
CLASSES=('POLYMORPHIC421','MONOMORPHIC_ALT4')

def rows(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path,encoding='utf-8')) as f:
        yield from csv.DictReader(f,delimiter='\t')

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def require(value,message):
    if not value:raise ValueError(message)

def write_tsv(path, records):
    require(bool(records),'No summary rows')
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)

def summarize(sd13,out):
    sd13=Path(sd13).resolve();out=Path(out).resolve()
    result=sd13/'results' if (sd13/'results').is_dir() else sd13
    require(not out.exists(),'Choose a new output directory')
    require(out!=result and result not in out.parents,'Do not write into source results')
    targets={r['site_index0']:r for r in rows(result/'targets.tsv')}
    samples=list(rows(result/'samples.tsv'));ids=[r['fsp'] for r in samples]
    require(len(ids)==len(set(ids))==142 and EXCLUDED<=set(ids),'Expected frozen142 sample set')
    require(len(targets)==425,'Expected425 targets')
    nclass=collections.Counter('POLYMORPHIC421' if r['polymorphic_in_focus']=='1' else 'MONOMORPHIC_ALT4' for r in targets.values())
    require(nclass=={'POLYMORPHIC421':421,'MONOMORPHIC_ALT4':4},'Expected421/4 split')
    counts=collections.defaultdict(collections.Counter);byiso=collections.defaultdict(collections.Counter)
    # For every site, preserve supported alleles regardless of baseline separately
    # from supported alleles matching the original retained calls.
    site=collections.defaultdict(collections.Counter);base=collections.defaultdict(lambda:collections.defaultdict(set))
    altiso=collections.defaultdict(set);unresolvedALT=collections.defaultdict(set)
    total=0
    for sid in ids:
        group='COMPARISON8' if sid in EXCLUDED else 'FOCUS134';seen=set()
        for r in rows(result/'samples'/sid/'support.tsv.gz'):
            key=(r['site_index0'],r['profile']);require(key not in seen,'Duplicate support row');seen.add(key)
            require(r['site_index0'] in targets and r['profile'] in PROFILES,'Unexpected target/profile')
            require(r['baseline'] in ('.','0','1') and r['supported_allele'] in ('.','0','1'),'Invalid allele state')
            t=targets[r['site_index0']];cat='POLYMORPHIC421' if t['polymorphic_in_focus']=='1' else 'MONOMORPHIC_ALT4'
            total+=1
            for g in (group,'ALL142'):
                k=(g,r['profile'],cat,r['baseline']);c=counts[k];c['cells']+=1;c[r['baseline_relation']]+=1;c['one_strand_only']+=int(r['one_strand_only']);c['zero_depth']+=int(r['depth']=='0')
                sk=(g,r['profile'],cat);site[(g,r['profile'],r['site_index0'])][r['supported_allele']]+=1
                if r['baseline_relation']=='SAME_ALLELE':
                    base[sk][r['baseline']].add(r['site_index0'])
                    if r['baseline']=='1':altiso[sk].add(sid)
                if r['baseline']=='1' and r['baseline_relation']=='UNRESOLVED_READ_SUPPORT':unresolvedALT[sk].add(r['site_index0'])
            if r['baseline']=='1':
                c=byiso[(sid,group,r['profile'],cat)];c['baseline_ALT']+=1;c[r['baseline_relation']]+=1
        require(len(seen)==425*3,'Incomplete sample support matrix: '+sid)
    require(total==181050,'Unexpected support cell count')
    site_rows=[]
    for g in ('FOCUS134','COMPARISON8','ALL142'):
        for profile in PROFILES:
            for cat in CLASSES:
                target_ids=[i for i,t in targets.items() if ('POLYMORPHIC421' if t['polymorphic_in_focus']=='1' else 'MONOMORPHIC_ALT4')==cat]
                counters=[site[(g,profile,i)] for i in target_ids];k=(g,profile,cat)
                refset=base[k]['0'];altset=base[k]['1']
                site_rows.append(dict(group=g,profile=profile,site_class=cat,sites=len(target_ids),sites_any_supported_REF=sum(c['0']>0 for c in counters),sites_any_supported_ALT=sum(c['1']>0 for c in counters),sites_both_supported_any_baseline=sum(c['0']>0 and c['1']>0 for c in counters),sites_with_baseline_ALT_supported=len(altset),sites_both_baseline_REF_ALT_supported=len(refset&altset),isolates_with_baseline_ALT_supported=len(altiso[k]),sites_any_unresolved_baseline_ALT=len(unresolvedALT[k]),sites_no_baseline_ALT_supported=len(target_ids)-len(altset)))
    fields=['cells','SAME_ALLELE','OPPOSITE_ALLELE','UNRESOLVED_READ_SUPPORT','BASELINE_MISSING','one_strand_only','zero_depth']
    support_rows=[dict(zip(['group','profile','site_class','baseline'],k),**{f:c[f] for f in fields}) for k,c in sorted(counts.items())]
    iso_rows=[dict(zip(['fsp','group','profile','site_class'],k),**{f:c[f] for f in ['baseline_ALT','SAME_ALLELE','OPPOSITE_ALLELE','UNRESOLVED_READ_SUPPORT']}) for k,c in sorted(byiso.items())]
    table=[]
    for profile in PROFILES:
        p=counts[('FOCUS134',profile,'POLYMORPHIC421','1')];m=counts[('FOCUS134',profile,'MONOMORPHIC_ALT4','1')]
        ps=next(r for r in site_rows if r['group']=='FOCUS134' and r['profile']==profile and r['site_class']=='POLYMORPHIC421')
        table.append(dict(profile=profile,focus_isolates=134,target_sites=425,all_ALT_supported=p['SAME_ALLELE']+m['SAME_ALLELE'],all_ALT_calls=p['cells']+m['cells'],polymorphic_ALT_supported=p['SAME_ALLELE'],polymorphic_ALT_calls=p['cells'],monomorphic_ALT_supported=m['SAME_ALLELE'],monomorphic_ALT_calls=m['cells'],polymorphic_sites_both_alleles_supported=ps['sites_both_baseline_REF_ALT_supported'],polymorphic_sites=421,polymorphic_sites_any_unresolved_ALT=ps['sites_any_unresolved_baseline_ALT'],polymorphic_sites_without_supported_ALT=ps['sites_no_baseline_ALT_supported'],polymorphic_unresolved_ALT_calls=p['UNRESOLVED_READ_SUPPORT'],monomorphic_unresolved_ALT_calls=m['UNRESOLVED_READ_SUPPORT'],all_opposite_allele_calls=p['OPPOSITE_ALLELE']+m['OPPOSITE_ALLELE']))
    out.mkdir(parents=True)
    for name,recs in [('table3_focus134_read_support.tsv',table),('support_by_polymorphism.tsv',support_rows),('site_support_by_polymorphism.tsv',site_rows),('isolate_ALT_support_by_polymorphism.tsv',iso_rows)]:write_tsv(out/name,recs)
    summary=dict(status='PASS_SUMMARY_RECOUNT',support_cells=total,target_sites=425,focus_isolates=134,comparison_isolates=8,focus_polymorphic_sites=421,focus_monomorphic_ALT_sites=4,profiles=list(PROFILES),table3=table,mask=json.load(open(result/'common_mask_scope.json')),input_results_manifest_sha256=sha(result/'SHA256SUMS.txt'),input_raw_parser_report=json.load(open(result/'independent_verification.json')),scope='Recount of exported labels; independent raw-pileup verification must precede this step. Same BAMs, not independent sequencing, association, or contamination diagnosis. No independence-based confidence interval calculated.')
    (out/'targeted_support_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sd13',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();print(json.dumps(summarize(a.sd13,a.out),indent=2))
