#!/usr/bin/env python3
"""Independent descriptive recount from frozen SD12 exports; no production imports.

Pairwise counts use jointly retained calls only and record denominators. Read support
is recomputed from the exported base-count columns, not independently extracted BAMs.
"""
import argparse,csv,gzip,hashlib,itertools,json,statistics
from collections import Counter,defaultdict
from pathlib import Path

EIGHT=('FSP348','FSP350','FSP351','FSP353','FSP357','FSP360','FSP363','FSP369')
FOUR_A=set(EIGHT[:4]);FOUR_B=set(EIGHT[4:]);EIGHT_SET=set(EIGHT)
HOSTS=('SC326-6','SC328-C','SC414-12E','BTX623','SC112-12E')

def rows(path):
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',newline='') as f:yield from csv.DictReader(f,delimiter='\t')

def write_tsv(path,entries,columns):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(entries)

def distribution(values):
    return {'n':len(values),'min':min(values),'median':statistics.median(values),'max':max(values)}

def bit_count(value):
    # Python 3.9 remains available on the recorded Ceres environment.
    return value.bit_count() if hasattr(int,'bit_count') else bin(value).count('1')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();data=args.data;out=args.out;out.mkdir(parents=True,exist_ok=False)
    ids=[r['fsp'] for r in rows(data/'01_samples/matrix_samples.tsv')]
    if len(ids)!=142 or not EIGHT_SET.issubset(ids):raise RuntimeError('Unexpected sample list')
    other=[s for s in ids if s not in EIGHT_SET]
    count=Counter();mac=Counter();call_per_sample={s:Counter() for s in ids}
    alt_bits={s:bytearray() for s in ids};called_bits={s:bytearray() for s in ids}
    targets=[];all_alt_outside=[];by_cds=Counter();titv=defaultdict(Counter)
    for index,r in enumerate(rows(data/'02_variants/selected_panel_genotypes.tsv.gz')):
        if index%8==0:
            for s in ids:alt_bits[s].append(0);called_bits[s].append(0)
        mask=1<<(index%8);byte=index//8
        for s in ids:
            g=r[s]
            if g not in ('0','1','.'):raise RuntimeError('Unexpected genotype')
            call_per_sample[s][g]+=1
            if g!='.':called_bits[s][byte]|=mask
            if g=='1':alt_bits[s][byte]|=mask
        nalt8=sum(r[s]=='1' for s in EIGHT);nref8=sum(r[s]=='0' for s in EIGHT)
        nalt134=sum(r[s]=='1' for s in other);nref134=sum(r[s]=='0' for s in other)
        total_alt=nalt8+nalt134;total_ref=nref8+nref134
        if total_alt!=int(r['n_ALT']) or total_ref!=int(r['n_REF']):
            raise RuntimeError('Stored site allele counts disagree')
        if min(total_alt,total_ref)==0:raise RuntimeError('Selected site is not polymorphic')
        mutation='transition' if frozenset((r['ref'],r['alt'])) in (frozenset('AG'),frozenset('CT')) else 'transversion'
        titv['ALL'][mutation]+=1;titv['CDS' if r['in_CDS']=='1' else 'NON_CDS'][mutation]+=1
        mac[min(total_alt,total_ref)]+=1
        count['selected_sites']+=1
        if nalt134==0:
            key='ALT_only_in_eight'
        elif nalt8==0:key='ALT_only_in_other134'
        else:key='ALT_in_both_groups'
        count[key]+=1;by_cds[(key,r['in_CDS'])]+=1
        if nalt8==8:count['ALT_in_all_eight']+=1
        if nalt134:
            entry={k:r[k] for k in ['site_index0','chrom','pos1','ref','alt','in_CDS']}
            entry.update(n_ALT_eight=nalt8,n_REF_eight=nref8,n_ALT_other134=nalt134,
                         n_REF_other134=nref134,n_missing_other134=134-nalt134-nref134,
                         polymorphic_other134=int(nref134>0),MAC_other134=min(nalt134,nref134))
            targets.append(entry)
            if nref134>0:
                count['polymorphic_within_other134']+=1
                if min(nref134,nalt134)>=2:count['polymorphic_other134_MAC_at_least_2']+=1
            else:all_alt_outside.append(entry)
    ai={s:int.from_bytes(alt_bits[s],'little') for s in ids}
    ci={s:int.from_bytes(called_bits[s],'little') for s in ids}
    pairs=[];group_values=defaultdict(list);group_comparable=defaultdict(list)
    for left,right in itertools.combinations(ids,2):
        joint=ci[left]&ci[right];n=bit_count(joint);d=bit_count((ai[left]^ai[right])&joint)
        if left in other and right in other:kind='within_other134'
        elif left in EIGHT_SET and right in EIGHT_SET:
            kind='within_four' if (left in FOUR_A)==(right in FOUR_A) else 'between_fours'
        else:kind='eight_vs_other134'
        pairs.append(dict(left=left,right=right,pair_group=kind,differences=d,jointly_retained_sites=n,difference_fraction=d/n if n else ''))
        group_values[kind].append(d);group_comparable[kind].append(n)
    isolate_rows=[]
    for s in ids:
        c=call_per_sample[s]
        isolate_rows.append(dict(fsp=s,descriptive_group='eight_A' if s in FOUR_A else 'eight_B' if s in FOUR_B else 'other134',n_REF=c['0'],n_ALT=c['1'],n_missing=c['.'],missing_fraction=c['.']/count['selected_sites'],ALT_fraction_called=c['1']/(c['0']+c['1'])))
    supports=defaultdict(Counter);support_mismatches=0;strand_mismatches=0;all_opposite=0;support_quality=defaultdict(Counter)
    for r in rows(data/'03_read_audit/read_support_long.tsv.gz'):
        n=sum(int(r[k]) for k in ['REF_F','REF_R','ALT_F','ALT_R','OTHER_F','OTHER_R','N','DEL','SKIP'])
        ref=int(r['REF_F'])+int(r['REF_R']);alt=int(r['ALT_F'])+int(r['ALT_R'])
        called='0' if n>=5 and 10*ref>=9*n else '1' if n>=5 and 10*alt>=9*n else '.'
        expected={'SUPPORTED_REF':'0','SUPPORTED_ALT':'1'}.get(r['support_status'],'.')
        if called!=expected:support_mismatches+=1
        allele_base='REF' if called=='0' else 'ALT'
        strand=int(called in ('0','1') and (int(r[allele_base+'_F'])==0 or int(r[allele_base+'_R'])==0))
        if strand!=int(r['one_strand_only']):strand_mismatches+=1
        for scope in ('ALL_AUDIT_CELLS','PRIMARY_SELECTED' if r['primary_sample']=='1' and r['partition']=='SELECTED' else ''):
            if not scope:continue
            q=support_quality[(scope,r['profile'])]
            q['cells']+=1;q[r['support_status']]+=1;q['supported_cells']+=called in ('0','1');q['one_strand_supported_cells']+=strand
        if r['baseline'] in ('0','1') and called not in (r['baseline'],'.'):all_opposite+=1
        if r['primary_sample']=='1' and r['baseline'] in ('0','1'):
            group='eight' if r['fsp'] in EIGHT_SET else 'other134'
            c=supports[(r['partition'],r['profile'],group,r['baseline'])]
            c['baseline']+=1;c['same']+=called==r['baseline'];c['unresolved']+=called=='.';c['opposite']+=called not in (r['baseline'],'.')
    support_rows=[dict(partition=k[0],profile=k[1],isolate_group=k[2],baseline_allele=k[3],**v) for k,v in sorted(supports.items())]
    master=list(rows(data/'01_samples/all245_sample_master.tsv'))
    profiles=Counter(tuple(r[h] for h in HOSTS[:3]) for r in master)
    v_rows=[r['fsp'] for r in master if any(r[h]=='V' for h in HOSTS)]
    filt=Counter(r['FILTER'] for r in rows(data/'02_variants/emitted_biallelic_catalogue.tsv.gz'))
    for quality in support_quality.values():
        for label in ['SUPPORTED_REF','SUPPORTED_ALT','LOW_DEPTH','NO_RETAINED_OBSERVATION','MIXED_OR_NONBASE','OTHER_BASES_DOMINATE']:
            quality.setdefault(label,0)
    result={
        'status':'COMPLETE','eight_ids':EIGHT,'counts':dict(count),'fraction_ALT_only_eight':count['ALT_only_in_eight']/count['selected_sites'],
        'counts_by_CDS':[dict(classification=k[0],in_CDS=k[1],sites=v) for k,v in sorted(by_cds.items())],
        'MAC_all142':dict(sorted(mac.items())),
        'TiTv_selected_sites':{k:dict(**v,ratio=v['transition']/v['transversion']) for k,v in titv.items()},
        'MAC_bins':{'1':mac[1],'2':mac[2],'3_5':sum(mac[k] for k in range(3,6)),'6_14':sum(mac[k] for k in range(6,15)),'15plus':sum(v for k,v in mac.items() if k>=15)},
        'other134_ALT_site_count':len(targets),'other134_monomorphic_ALT_sites':all_alt_outside,
        'ALT_per_isolate':{'eight':distribution([call_per_sample[s]['1'] for s in EIGHT]),'other134':distribution([call_per_sample[s]['1'] for s in other])},
        'isolates_with_ALT_at_most_24':[s for s in ids if call_per_sample[s]['1']<=24],
        'other134_ALT_outlier':[{k:r[k] for k in ['fsp','n_ALT','ALT_fraction_called']} for r in isolate_rows if r['fsp'] in other and r['n_ALT']>24],
        'selected_missingness':dict(distribution([r['missing_fraction'] for r in isolate_rows]),min_isolates=[r['fsp'] for r in isolate_rows if r['n_missing']==min(x['n_missing'] for x in isolate_rows)],max_isolates=[r['fsp'] for r in isolate_rows if r['n_missing']==max(x['n_missing'] for x in isolate_rows)],missing_cells=sum(r['n_missing'] for r in isolate_rows),matrix_cells=count['selected_sites']*len(ids)),
        'pairwise':{k:{'differences':distribution(v),'jointly_retained_sites':distribution(group_comparable[k]),'zero_difference_pairs':sum(x==0 for x in v)} for k,v in group_values.items()},
        'support_recount_from_exported_base_counts':{'mismatches':support_mismatches,'strand_mismatches':strand_mismatches,'opposite_all_audit':all_opposite},
        'support_quality_by_profile':[dict(scope=k[0],profile=k[1],**v,mixed_fraction_all_cells=v['MIXED_OR_NONBASE']/v['cells'],one_strand_fraction_supported_cells=v['one_strand_supported_cells']/v['supported_cells']) for k,v in sorted(support_quality.items())],
        'primary_support_by_carrier_group':support_rows,
        'host_profiles_three_variable_hosts':len(profiles),'isolates_with_any_V':len(v_rows),
        'host_profile_counts':[dict(profile='/'.join(k),isolates=v) for k,v in sorted(profiles.items())],
        'FILTER_values_full_catalogue':dict(filt),
        'interpretation':'Descriptive retained-call groupings only; no taxonomy/clonality assignment. 425 sites with ALT outside the eight includes four monomorphic-ALT sites in other134, so these are not all rare polymorphisms within that subset.',
        'scope':'Independent recount of frozen exported tables; no production imports, no new BAM extraction or raw-pileup token parsing.'}
    (out/'structure_recount.json').write_text(json.dumps(result,indent=2)+'\n')
    write_tsv(out/'other134_ALT_targets.tsv',targets,list(targets[0]))
    write_tsv(out/'primary142_isolate_call_counts.tsv',isolate_rows,list(isolate_rows[0]))
    write_tsv(out/'pairwise_differences_and_denominators.tsv',pairs,list(pairs[0]))
    write_tsv(out/'primary_support_by_carrier_group.tsv',support_rows,list(support_rows[0]))
    quality_columns=['scope','profile','cells','SUPPORTED_REF','SUPPORTED_ALT','LOW_DEPTH','NO_RETAINED_OBSERVATION','MIXED_OR_NONBASE','OTHER_BASES_DOMINATE','supported_cells','one_strand_supported_cells','mixed_fraction_all_cells','one_strand_fraction_supported_cells']
    write_tsv(out/'support_quality_by_profile.tsv',result['support_quality_by_profile'],quality_columns)
    hashes=[]
    for rel in ['01_samples/matrix_samples.tsv','01_samples/all245_sample_master.tsv','02_variants/selected_panel_genotypes.tsv.gz','02_variants/emitted_biallelic_catalogue.tsv.gz','03_read_audit/read_support_long.tsv.gz']:
        hashes.append(dict(file=rel,sha256=hashlib.sha256((data/rel).read_bytes()).hexdigest()))
    (out/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['counts','MAC_bins','host_profiles_three_variable_hosts','isolates_with_any_V','FILTER_values_full_catalogue']},indent=2))

if __name__=='__main__':main()
