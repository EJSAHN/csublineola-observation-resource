#!/usr/bin/env python3
"""Independent recount of targeted ALT pileups and group summaries.
Parser functions reused from SD12 verify_sd11.py; aggregation adapted for SD13.
No imports from the production package; standard-library Python only.
"""
import argparse,collections,csv,gzip,hashlib,json,math,re,statistics
from pathlib import Path

def load_tsv(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def table(p,rs):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def check(a,b,label):
    if str(a)!=str(b):raise AssertionError(f'{label}: {a!r} != {b!r}')
COUNT=['REF_F','REF_R','ALT_F','ALT_R','OTHER_F','OTHER_R','N','DEL','SKIP']
OTHER=['depth','base_quality_sum','base_quality_n','map_quality_sum','pileup_row_present']
PROFILES=['MQ20_BQ20_noBAQ','MQ40_BQ20_noBAQ','MQ40_BQ20_redoBAQ']

def unpack(line,t):
    a=line.rstrip('\n').split('\t');assert len(a)==7
    check(a[2].upper(),t['ref'],'pileup reference');n=int(a[3]);d={k:0 for k in COUNT+OTHER};d.update(depth=n,pileup_row_present=1)
    if n==0:
        assert a[4:]==['*','*','*'];return d
    text=re.sub(r'\^.', '',a[4]).replace('$','')
    while True:
        m=re.search(r'[+-](\d+)',text)
        if not m:break
        end=m.end()+int(m.group(1));assert end<=len(text);text=text[:m.start()]+text[end:]
    assert len(text)==len(a[5])==len(a[6])==n
    for c,q,m in zip(text,a[5],a[6]):
        d['map_quality_sum']+=ord(m)-33
        if c in '*#':d['DEL']+=1
        elif c in '<>':d['SKIP']+=1
        elif c in 'Nn':d['N']+=1
        else:
            assert c in '.,ACGTacgt',c
            strand='F' if c=='.' or c.isupper() else 'R'
            base=t['ref'] if c in '.,' else c.upper()
            key='REF' if base==t['ref'] else 'ALT' if base==t['alt'] else 'OTHER'
            d[key+'_'+strand]+=1;d['base_quality_sum']+=ord(q)-33;d['base_quality_n']+=1
    assert sum(d[k] for k in COUNT)==n
    return d

def label(d,baseline):
    n=d['depth'];a='.'
    if n==0:s='NO_RETAINED_OBSERVATION'
    elif n<5:s='LOW_DEPTH'
    elif 10*(d['REF_F']+d['REF_R'])>=9*n:s='SUPPORTED_REF';a='0'
    elif 10*(d['ALT_F']+d['ALT_R'])>=9*n:s='SUPPORTED_ALT';a='1'
    elif 10*(d['OTHER_F']+d['OTHER_R'])>=9*n:s='OTHER_BASES_DOMINATE'
    else:s='MIXED_OR_NONBASE'
    relation='BASELINE_MISSING' if baseline=='.' else 'UNRESOLVED_READ_SUPPORT' if a=='.' else 'SAME_ALLELE' if a==baseline else 'OPPOSITE_ALLELE'
    side='REF' if a=='0' else 'ALT'
    return dict(baseline=baseline,support_status=s,supported_allele=a,baseline_relation=relation,one_strand_only=int(a!='.' and min(d[side+'_F'],d[side+'_R'])==0))

def run(result,excluded):
    if not __debug__:raise RuntimeError('Independent verifier requires Python without -O')
    result=Path(result);excluded=set(excluded)
    targets=load_tsv(result/'targets.tsv');samples=load_tsv(result/'samples.tsv')
    target={int(t['site_index0']):t for t in targets};coords={(t['chrom'],int(t['pos1'])):int(t['site_index0']) for t in targets}
    baseline={int(t['site_index0']):t for t in load_tsv(result/'target_genotypes.tsv.gz')}
    aggregate=collections.defaultdict(collections.Counter);site=collections.defaultdict(collections.Counter);n=fields=0
    for sm in samples:
        sid=sm['fsp'];group='COMPARISON8' if sid in excluded else 'FOCUS134';folder=result/'samples'/sid
        reported=load_tsv(folder/'support.tsv.gz');rep={(int(r['site_index0']),r['profile']):r for r in reported}
        assert len(rep)==len(reported)==len(targets)*3
        for profile in PROFILES:
            obs={}
            with gzip.open(folder/(profile+'.pileup.gz'),'rt') as f:
                for line in f:
                    a=line.split('\t',2);i=coords[(a[0],int(a[1]))];assert i not in obs;obs[i]=unpack(line,target[i])
            for i,t in target.items():
                d=obs.get(i,{k:0 for k in COUNT+OTHER});d=dict(d);d.update(label(d,baseline[i][sid]));r=rep[(i,profile)];n+=1
                for k,v in d.items():check(v,r[k],(sid,i,profile,k));fields+=1
                for k in ('chrom','pos1','ref','alt','partition','in_CDS','primary_sample','witness'):check(t[k],r[k],(sid,i,k));fields+=1
                for g in (group,'ALL142'):
                    c=aggregate[(g,profile,d['baseline'])];c['cells']+=1;c[d['baseline_relation']]+=1;c['one_strand_only']+=d['one_strand_only'];c['zero_depth']+=int(d['depth']==0)
                    for allele,name in [('0','supported_REF'),('1','supported_ALT'),('.','unresolved')]:c[name]+=int(d['supported_allele']==allele)
                    site[(g,i,profile)][d['supported_allele']]+=1
    summary=load_tsv(result/'support_by_group_and_baseline.tsv');assert len(summary)==len(aggregate)
    for r in summary:
        key=tuple(r[k] for k in ('group','profile','baseline'));c=aggregate[key]
        for k in list(r)[3:]:
            if k=='same_fraction':
                if r['baseline']=='.':check(r[k],'','missing fraction')
                else:assert math.isclose(float(r[k]),c['SAME_ALLELE']/c['cells'],abs_tol=1e-14)
            else:check(r[k],c[k],(key,k))
    ss=load_tsv(result/'site_support_by_group.tsv');assert len(ss)==len(site)
    for r in ss:
        c=site[(r['group'],int(r['site_index0']),r['profile'])]
        for k,v in [('supported_REF',c['0']),('supported_ALT',c['1']),('unresolved',c['.']),('both_alleles_supported',int(c['0']>0 and c['1']>0))]:check(v,r[k],('site summary',k))
    report=dict(status='PASS',support_cells=n,support_and_metadata_fields_compared=fields,group_baseline_summary_rows=len(summary),site_group_profile_rows=len(ss),imports_production_code=False,mismatches=0,scope='Independent parsing of archived pileups and group summaries; same BAMs, not independent sequencing validation.')
    save(result/'independent_verification.json',report);return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('result_directory',type=Path);p.add_argument('--exclude',nargs='+',required=True);a=p.parse_args();print(json.dumps(run(a.result_directory,a.exclude),indent=2))
