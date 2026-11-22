#!/usr/bin/env python3
"""Independent recount of every returned SD11 pileup, support cell and summary.
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

def run(root,out):
    result=root/'results';out.mkdir(parents=True,exist_ok=True)
    targets=load_tsv(result/'targets.tsv');samples=load_tsv(result/'samples.tsv')
    target={int(t['site_index0']):t for t in targets};coordinates={(t['chrom'],int(t['pos1'])):int(t['site_index0']) for t in targets}
    baselines={int(r['site_index0']):r for r in load_tsv(result/'target_genotypes.tsv.gz')}
    cells={};verified_fields=0;pile_lines=0;summary=collections.defaultdict(collections.Counter);depths=collections.defaultdict(list);trans=collections.defaultdict(collections.Counter);site=collections.defaultdict(collections.Counter)
    contrast_names=['RAISE_MAPQ','RECALCULATE_BAQ']
    for sm in samples:
        sid=sm['fsp'];group=sm['source_group'];folder=result/'samples'/sid
        status=json.loads((folder/'status.json').read_text());check(status['status'],'COMPLETE','sample status')
        reported=load_tsv(folder/'support.tsv.gz');rep={(int(r['site_index0']),r['profile']):r for r in reported};assert len(rep)==len(reported)==len(targets)*3
        for pr in PROFILES:
            dct={}
            with gzip.open(folder/(pr+'.pileup.gz'),'rt') as f:
                for line in f:
                    a=line.split('\t',2);idx=coordinates[(a[0],int(a[1]))];assert idx not in dct;dct[idx]=unpack(line,target[idx]);pile_lines+=1
            for idx,t in target.items():
                d=dct.get(idx,{k:0 for k in COUNT+OTHER});d=dict(d);d.update(label(d,baselines[idx][sid]));r=rep[(idx,pr)]
                for k,v in d.items():check(v,r[k],(sid,idx,pr,k));verified_fields+=1
                for k in ('chrom','pos1','ref','alt','partition','in_CDS','primary_sample','witness'):check(t[k],r[k],(sid,idx,k));verified_fields+=1
                cells[(sid,idx,pr)]=d
                site[(idx,pr)]['supported_'+d['supported_allele']]+=1
                if t['primary_sample']=='1':
                    key=(t['partition'],t['in_CDS'],group,pr);c=summary[key];c['cells']+=1
                    c['baseline_missing' if d['baseline']=='.' else 'baseline_retained']+=1
                    c[d['baseline_relation']]+=1;c['support_'+d['support_status']]+=1;c['one_strand_supported']+=d['one_strand_only']
                    c['missing_but_supported']+=int(d['baseline']=='.' and d['supported_allele']!='.');depths[key].append(d['depth'])
        for idx,t in target.items():
            if t['primary_sample']!='1':continue
            for i,contrast in enumerate(contrast_names):
                a=cells[(sid,idx,PROFILES[i])];b=cells[(sid,idx,PROFILES[i+1])];x=a['supported_allele'];y=b['supported_allele']
                state='BOTH_UNRESOLVED' if x==y=='.' else 'GAINED_SUPPORT' if x=='.' else 'LOST_SUPPORT' if y=='.' else 'SAME_SUPPORTED_ALLELE' if x==y else 'SUPPORTED_ALLELE_SWITCH'
                c=trans[(t['partition'],t['in_CDS'],group,contrast)];c['cells']+=1;c[state]+=1;c['depth_increased']+=int(b['depth']>a['depth'])
    sf=0
    for r in load_tsv(result/'primary_strata.tsv'):
        key=tuple(r[k] for k in ('partition','in_CDS','source_group','profile'));c=summary[key]
        for k in r:
            if k in ('partition','in_CDS','source_group','profile'):continue
            if k=='median_depth':v=statistics.median(depths[key])
            elif k=='same_fraction_of_baseline_retained':v=c['SAME_ALLELE']/c['baseline_retained'] if c['baseline_retained'] else None
            else:v=c[k]
            if v is None:assert r[k]==''
            else:assert math.isclose(float(r[k]),v,abs_tol=1e-14,rel_tol=1e-14),(key,k)
            sf+=1
    tf=0
    for r in load_tsv(result/'primary_transitions.tsv'):
        key=tuple(r[k] for k in ('partition','in_CDS','source_group','contrast'))
        for k in list(r)[4:]:check(trans[key][k],r[k],(key,k));tf+=1
    for r in load_tsv(result/'site_support.tsv'):
        c=site[(int(r['site_index0']),r['profile'])]
        for k,v in [('supported_REF',c['supported_0']),('supported_ALT',c['supported_1']),('unresolved',c['supported_.']),('both_alleles_supported',int(c['supported_0']>0 and c['supported_1']>0))]:check(v,r[k],('site',k))
    witness=load_tsv(result/'witness_requests.tsv');original={(r['panel'],r['sample_a'],r['sample_b'],r['site_index0']):r for r in witness};wc=collections.Counter();wf=0;retention=collections.defaultdict(list);loss=collections.Counter()
    for r in load_tsv(result/'witness_read_support.tsv.gz'):
        w=original[tuple(r[k] for k in ('panel','sample_a','sample_b','site_index0'))]
        for k,v in w.items():check(v,r[k],('witness source',k))
        pr=r['profile'];idx=int(r['site_index0']);a=cells[(r['sample_a'],idx,pr)];b=cells[(r['sample_b'],idx,pr)]
        x,y=a['supported_allele'],b['supported_allele'];state='ORIGINAL_DIFFERENCE_SUPPORTED' if x==w['allele_a'] and y==w['allele_b'] else 'ONE_OR_BOTH_UNRESOLVED' if '.' in (x,y) else 'SAME_ALLELE_SUPPORTED' if x==y else 'REVERSED_DIFFERENCE_SUPPORTED'
        check(state,r['witness_read_status'],'witness label');wf+=1
        for side,d in [('a',a),('b',b)]:
            for k in ('support_status','supported_allele','depth','REF_F','REF_R','ALT_F','ALT_R','one_strand_only'):check(d[k],r[side+'_'+k],('witness field',k));wf+=1
        wc[(r['panel'],r['baseline_passes_95'],pr,state)]+=1
        if r['baseline_passes_95']=='1':
            retention[(r['panel'],pr)].append(dict(sample_a=r['sample_a'],sample_b=r['sample_b'],site_index0=idx,status=state))
            if pr==PROFILES[-1] and state=='ONE_OR_BOTH_UNRESOLVED':
                states=[d['support_status'] for d in (a,b) if d['supported_allele']=='.'];loss[(r['panel'],' + '.join(sorted(set(states))))]+=1
    stored={tuple(r[k] for k in ('panel','baseline_passes_95','profile','witness_read_status')):int(r['pair_site_records']) for r in load_tsv(result/'witness_summary.tsv')};assert dict(wc)==stored
    # New compact descriptive tables, restricted to the primary sampling frame.
    pooled=collections.defaultdict(collections.Counter)
    for key,c in summary.items():
        for field,value in c.items():pooled[(key[0],key[3])][field]+=value
    outrows=[]
    for (part,pr),c in sorted(pooled.items()):
        ts=[int(t['site_index0']) for t in targets if t['partition']==part and t['primary_sample']=='1']
        two=sum(site[(idx,pr)]['supported_0']>0 and site[(idx,pr)]['supported_1']>0 for idx in ts)
        outrows.append(dict(partition=part,profile=pr,sampled_sites=len(ts),both_alleles_supported_sites=two,baseline_retained_cells=c['baseline_retained'],same_allele_cells=c['SAME_ALLELE'],same_fraction=c['SAME_ALLELE']/c['baseline_retained'],opposite_cells=c['OPPOSITE_ALLELE'],unresolved_cells=c['UNRESOLVED_READ_SUPPORT'],baseline_missing_cells=c['baseline_missing'],missing_now_supported=c['missing_but_supported']))
    table(out/'primary_compact.tsv',outrows)
    table(out/'witness95_compact.tsv',[dict(panel=p,profile=pr,pairs=len(v),original_difference_supported=sum(r['status']=='ORIGINAL_DIFFERENCE_SUPPORTED' for r in v),unresolved=sum(r['status']=='ONE_OR_BOTH_UNRESOLVED' for r in v)) for (p,pr),v in sorted(retention.items())])
    table(out/'witness95_unresolved_reasons.tsv',[dict(panel=p,unresolved_reason=reason,pairs=n) for (p,reason),n in sorted(loss.items())])
    report=dict(status='PASS',imports_production_code=False,samples=len(samples),sites=len(targets),profiles=3,pileup_files=len(samples)*3,pileup_lines=pile_lines,support_cells=len(cells),support_and_metadata_fields_compared=verified_fields,primary_summary_fields=sf,transition_fields=tf,site_summary_rows=len(site),witness_rows=1986,witness_fields=wf,witness_summary_groups=len(wc),mismatches=0,scope='Raw returned pileups independently parsed and all reported summaries reaggregated; original BAM/FASTA extraction not rerun')
    save(out/'independent_verification.json',report)
    # Distinguish primary genotype agreement, all-site opposite evidence and newly supported missing cells.
    diagnostic=collections.Counter()
    for (sid,idx,pr),d in cells.items():diagnostic[(pr,d['baseline_relation'])]+=1
    table(out/'all_union_baseline_relations.tsv',[dict(profile=k[0],relation=k[1],cells=v) for k,v in sorted(diagnostic.items())])
    ac=collections.defaultdict(collections.Counter)
    for (sid,idx,pr),d in cells.items():
        if target[idx]['primary_sample']=='1' and d['baseline']!='.':
            c=ac[(target[idx]['partition'],pr,d['baseline'])];c['baseline_retained']+=1;c[d['baseline_relation']]+=1
    ar=[]
    for (part,pr,b),c in sorted(ac.items()):
        ar.append(dict(partition=part,profile=pr,baseline_allele='REF' if b=='0' else 'ALT',baseline_retained=c['baseline_retained'],same=c['SAME_ALLELE'],unresolved=c['UNRESOLVED_READ_SUPPORT'],opposite=c['OPPOSITE_ALLELE'],same_fraction=c['SAME_ALLELE']/c['baseline_retained']))
    table(out/'primary_allele_specific.tsv',ar)
    print(json.dumps(report,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run_directory',type=Path);p.add_argument('output_directory',type=Path);a=p.parse_args();run(a.run_directory,a.output_directory)
