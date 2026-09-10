#!/usr/bin/env python3
import argparse, bisect, csv, gzip, math, os
from collections import defaultdict

def norm_chr(x):
    s=str(x).strip()
    if s.lower().startswith('chr'): s=s[3:]
    return s.upper()

def num(x, kind=float):
    try: return kind(float(x)) if kind is int else float(x)
    except: return None if kind is int else math.nan

def load_dmrs(path):
    out=[]
    with open(path, newline='', encoding='utf-8-sig') as f:
        r=csv.DictReader(f)
        req={'DMR_region_id','DMR_chr','DMR_region_start','DMR_region_end'}
        miss=req-set(r.fieldnames or [])
        if miss: raise ValueError(f'Missing DMR columns: {sorted(miss)}')
        for i,row in enumerate(r):
            s=num(row['DMR_region_start'],int); e=num(row['DMR_region_end'],int)
            if s is None or e is None: raise ValueError(f'Bad coordinates on row {i+2}')
            if e<s: s,e=e,s
            out.append({'idx':i,'id':row['DMR_region_id'],'chr':norm_chr(row['DMR_chr']),'start':s,'end':e})
    return out

def build_index(dmrs):
    by=defaultdict(list)
    for d in dmrs: by[d['chr']].append(d)
    idx={}
    for ch,a in by.items():
        a.sort(key=lambda z:(z['start'],z['end']))
        starts=[z['start'] for z in a]
        pref=[]; m=-1
        for z in a:
            m=max(m,z['end']); pref.append(m)
        idx[ch]=(a,starts,pref)
    return idx

def point_hits(idx,ch,pos):
    if ch not in idx: return []
    a,starts,pref=idx[ch]
    i=bisect.bisect_right(starts,pos)-1
    out=[]
    while i>=0 and pref[i]>=pos:
        d=a[i]
        if d['start']<=pos<=d['end']: out.append(d)
        i-=1
    return out

def main():
    ap=argparse.ArgumentParser(description='Stream EPIGEN meQTLs and overlap target CpGs with GAM-DMRs')
    ap.add_argument('--dmr',required=True)
    ap.add_argument('--meqtl',required=True)
    ap.add_argument('--out',default='EPIGEN_DMR_overlap')
    ap.add_argument('--progress-every',type=int,default=1000000)
    args=ap.parse_args()
    os.makedirs(args.out,exist_ok=True)
    dmrs=load_dmrs(args.dmr); idx=build_index(dmrs)
    print(f'Loaded {len(dmrs):,} DMRs')

    st={d['idx']:{'all_cpg':set(),'cis_cpg':set(),'trans_cpg':set(),'all_snp':set(),'cis_snp':set(),'trans_snp':set(),'ld':set(),
                  'n':0,'cis_n':0,'trans_n':0,'min_p':math.inf,'min_fdr':math.inf,
                  'top_p':math.inf,'top_cpg':'','top_snp':'','top_type':'','top_beta':''} for d in dmrs}

    assoc_path=os.path.join(args.out,'EPIGEN_DMR_meQTL_associations.tsv.gz')
    total=0; matched=0
    with gzip.open(args.meqtl,'rt',encoding='utf-8',newline='') as fi, gzip.open(assoc_path,'wt',encoding='utf-8',newline='') as fo:
        r=csv.DictReader(fi,delimiter='\t')
        req={'SNP','CpG','type','beta','se','p-value','FDR','chr_cpg','pos_cpg','chr_snp','pos_snp','n_studies','n_samples','effects','MAF','LD_clump'}
        miss=req-set(r.fieldnames or [])
        if miss: raise ValueError(f'Missing EPIGEN columns: {sorted(miss)}')
        fields=['DMR_region_id','DMR_chr','DMR_region_start','DMR_region_end']+list(r.fieldnames)
        w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t'); w.writeheader()
        for row in r:
            total+=1
            pos=num(row['pos_cpg'],int)
            if pos is None: continue
            hs=point_hits(idx,norm_chr(row['chr_cpg']),pos)
            if not hs:
                if args.progress_every and total%args.progress_every==0: print(f'Scanned {total:,}; matched {matched:,}')
                continue
            typ=row['type'].strip().lower(); p=num(row['p-value']); fdr=num(row['FDR'])
            for d in hs:
                matched+=1
                out={'DMR_region_id':d['id'],'DMR_chr':d['chr'],'DMR_region_start':d['start'],'DMR_region_end':d['end']}; out.update(row); w.writerow(out)
                s=st[d['idx']]; cpg=row['CpG']; snp=row['SNP']; ld=row.get('LD_clump','')
                s['all_cpg'].add(cpg); s['all_snp'].add(snp); s['n']+=1
                if ld: s['ld'].add(ld)
                if typ=='cis': s['cis_cpg'].add(cpg); s['cis_snp'].add(snp); s['cis_n']+=1
                elif typ=='trans': s['trans_cpg'].add(cpg); s['trans_snp'].add(snp); s['trans_n']+=1
                if math.isfinite(p): s['min_p']=min(s['min_p'],p)
                if math.isfinite(fdr): s['min_fdr']=min(s['min_fdr'],fdr)
                if math.isfinite(p) and p<s['top_p']:
                    s['top_p']=p; s['top_cpg']=cpg; s['top_snp']=snp; s['top_type']=typ; s['top_beta']=row['beta']
            if args.progress_every and total%args.progress_every==0: print(f'Scanned {total:,}; matched {matched:,}')

    rows=[]
    for d in dmrs:
        s=st[d['idx']]
        rows.append({
            'DMR_region_id':d['id'],'DMR_chr':d['chr'],'DMR_region_start':d['start'],'DMR_region_end':d['end'],
            'EPIGEN_meQTL_overlap':'Yes' if s['all_cpg'] else 'No',
            'EPIGEN_cis_meQTL_overlap':'Yes' if s['cis_cpg'] else 'No',
            'EPIGEN_trans_meQTL_overlap':'Yes' if s['trans_cpg'] else 'No',
            'n_unique_meQTL_CpGs':len(s['all_cpg']),'n_cis_meQTL_CpGs':len(s['cis_cpg']),'n_trans_meQTL_CpGs':len(s['trans_cpg']),
            'n_unique_meQTL_SNPs':len(s['all_snp']),'n_cis_meQTL_SNPs':len(s['cis_snp']),'n_trans_meQTL_SNPs':len(s['trans_snp']),
            'n_meQTL_associations':s['n'],'n_cis_associations':s['cis_n'],'n_trans_associations':s['trans_n'],'n_unique_LD_clumps':len(s['ld']),
            'minimum_p_value':'' if not math.isfinite(s['min_p']) else s['min_p'],
            'minimum_FDR':'' if not math.isfinite(s['min_fdr']) else s['min_fdr'],
            'top_meQTL_CpG':s['top_cpg'],'top_meQTL_SNP':s['top_snp'],'top_meQTL_type':s['top_type'],'top_meQTL_beta':s['top_beta'],
            'top_meQTL_p_value':'' if not math.isfinite(s['top_p']) else s['top_p']})

    summary=os.path.join(args.out,'EPIGEN_DMR_meQTL_overlap_summary.csv')
    with open(summary,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    supp_fields=['DMR_region_id','DMR_chr','DMR_region_start','DMR_region_end','EPIGEN_meQTL_overlap','EPIGEN_cis_meQTL_overlap','EPIGEN_trans_meQTL_overlap',
                 'n_unique_meQTL_CpGs','n_cis_meQTL_CpGs','n_trans_meQTL_CpGs','n_unique_meQTL_SNPs','n_unique_LD_clumps','minimum_p_value','minimum_FDR']
    supp=os.path.join(args.out,'EPIGEN_DMR_meQTL_overlap_supplement.csv')
    with open(supp,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=supp_fields); w.writeheader(); [w.writerow({k:r[k] for k in supp_fields}) for r in rows]

    n_any=sum(r['EPIGEN_meQTL_overlap']=='Yes' for r in rows); n_cis=sum(r['EPIGEN_cis_meQTL_overlap']=='Yes' for r in rows); n_trans=sum(r['EPIGEN_trans_meQTL_overlap']=='Yes' for r in rows)
    headline=os.path.join(args.out,'EPIGEN_DMR_overlap_overall_summary.txt')
    with open(headline,'w',encoding='utf-8') as f:
        f.write('EPIGEN meQTL overlap with GAM-DMRs\n=================================\n')
        f.write(f'EPIGEN association rows scanned: {total:,}\nAssociation-DMR overlap rows: {matched:,}\nTotal DMRs: {len(dmrs):,}\n')
        f.write(f'DMRs with >=1 EPIGEN meQTL target CpG: {n_any:,} ({100*n_any/len(dmrs):.2f}%)\n')
        f.write(f'DMRs with >=1 cis-meQTL target CpG: {n_cis:,} ({100*n_cis/len(dmrs):.2f}%)\n')
        f.write(f'DMRs with >=1 trans-meQTL target CpG: {n_trans:,} ({100*n_trans/len(dmrs):.2f}%)\n')
    print('\nCompleted')
    print(f'DMRs with >=1 meQTL target CpG: {n_any}/{len(dmrs)} ({100*n_any/len(dmrs):.2f}%)')
    print(f'DMRs with >=1 cis target CpG: {n_cis}/{len(dmrs)} ({100*n_cis/len(dmrs):.2f}%)')
    print(f'DMRs with >=1 trans target CpG: {n_trans}/{len(dmrs)} ({100*n_trans/len(dmrs):.2f}%)')
    print('Summary:',summary); print('Supplement:',supp); print('Matches:',assoc_path); print('Headline:',headline)

if __name__=='__main__': main()
