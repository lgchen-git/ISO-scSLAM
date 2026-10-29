#!/usr/bin/env python3
"""Rebuild the ISO-SLAM event-level AS tables, figures, and HTML report."""
from __future__ import annotations
import argparse, base64, csv, gzip, html, json, math, re, shutil
from collections import defaultdict, Counter
from pathlib import Path

SAMPLES=("A-1","A-2","A-3","AB-1","AB-2","AB-3","B-1","B-2","B-3","NC-1","NC-2","NC-3")
GROUP={s:s.split('-',1)[0] for s in SAMPLES}
CONTRASTS=(("A_vs_NC","A","NC"),("B_vs_NC","B","NC"),("AB_vs_NC","AB","NC"),("A_vs_B","A","B"),("AB_vs_A","AB","A"),("AB_vs_B","AB","B"))
TYPES=("SE","RI","A3","A5","MX","AF","AL")

def tsv(p):
    with p.open(encoding='utf-8',newline='') as h: return list(csv.DictReader(h,delimiter='\t'))
def write(p, rows, fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def config(p):
    out={}
    for line in p.read_text(encoding='utf-8').splitlines():
        if ':' in line and line.strip() and not line.lstrip().startswith('#') and not line[:1].isspace():
            k,v=line.split(':',1); out[k.strip()]=v.split('#',1)[0].strip().strip("'\"")
    return out

def sample_sheet(p):
    with p.open(encoding='utf-8', newline='') as h:
        rows=list(csv.DictReader(h, delimiter='\t'))
    if not rows: raise RuntimeError(f'empty sample sheet: {p}')
    ids=tuple((r.get('sample_id') or '').strip() for r in rows)
    if any(not x for x in ids) or len(set(ids)) != len(ids): raise RuntimeError('sample sheet has invalid or duplicate sample IDs')
    groups={x:(r.get('group') or '').strip() for x,r in zip(ids,rows)}
    if any(not groups[x] for x in ids): raise RuntimeError('sample sheet has an empty group')
    return ids, groups
def num(x):
    try: return float(x)
    except (TypeError,ValueError): return None
def uri(p): return 'data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode('ascii')
def mean(xs): return sum(xs)/len(xs) if xs else None

def restyle_slam_fraction(source_delivery, out):
    """Regenerate the SLAM fraction chart with the legend above the bars."""
    import matplotlib.pyplot as plt
    import numpy as np
    rows=tsv(source_delivery/'tables/read_classification_summary.tsv')
    names=[r['sample_id'] for r in rows]
    new=np.array([float(r['new_rna_reads']) for r in rows]); old=np.array([float(r['old_rna_reads']) for r in rows]); amb=np.array([float(r['ambiguous_reads']) for r in rows]); total=new+old+amb
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(10,4.8)); x=np.arange(len(names)); oldf=old/total; newf=new/total; ambf=amb/total
    ax.bar(x,oldf,label='Old',color='#3C5488'); ax.bar(x,newf,bottom=oldf,label='New',color='#E64B35'); ax.bar(x,ambf,bottom=oldf+newf,label='Ambiguous',color='#BDBDBD')
    ax.set_ylabel('Fraction of classified reads'); ax.set_title('ISO-SLAM read classification by sample',pad=6); ax.set_ylim(0,1); ax.set_xticks(x,names,rotation=35,ha='right')
    ax.legend(frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(0.5,1.14),borderaxespad=0)
    fig.subplots_adjust(top=0.72,bottom=0.25,left=0.09,right=0.98); out.mkdir(parents=True,exist_ok=True); fig.savefig(out/'new_rna_fraction.png',dpi=300); plt.close(fig)

def copy_source(src, root):
    for rel in ('tables/total_rna','tables/new_rna','total_rna_provenance','qc'):
        p=src/rel
        if p.is_dir(): shutil.copytree(p,root/rel,dirs_exist_ok=True)
    if (src/'figures/total_rna').is_dir(): shutil.copytree(src/'figures/total_rna',root/'figures/total_rna',dirs_exist_ok=True)

def event_definitions(ioe_dir, event_ids, root):
    found={}
    for p in sorted(ioe_dir.glob('*.ioe')):
        with p.open(encoding='utf-8',errors='replace') as h:
            for r in csv.DictReader(h,delimiter='\t'):
                eid=r.get('event_id','')
                if eid in event_ids and eid not in found:
                    found[eid]={'event_id':eid,'event_type':event_ids[eid]['event_type'],'gene_id':event_ids[eid]['gene_id'],
                                'inclusion_transcripts':r.get('alternative_transcripts',''),'exclusion_transcripts':r.get('total_transcripts','')}
    if len(found)!=len(event_ids): raise RuntimeError(f'event definition mismatch: {len(event_ids)-len(found)} missing')
    write(root/'provenance/event_definition.tsv',list(found.values()),['event_id','event_type','gene_id','inclusion_transcripts','exclusion_transcripts'])
    return found

def make_metadata(root):
    rows=[]
    for s in SAMPLES:
        rows.append({'sample_id':s,'group':GROUP[s],'biological_replicate':s.rsplit('-',1)[1],'library_id':'NA','notes':'Sample identifiers are pseudonymous analysis labels.'})
    write(root/'provenance/sample_metadata.tsv',rows,['sample_id','group','biological_replicate','library_id','notes'])

def new_tables(root, src):
    support=tsv(src/'tables/new_rna/event_support.tsv'); psi=tsv(src/'tables/new_rna/event_psi.tsv')
    event_meta={r['event_id']:{'event_id':r['event_id'],'event_type':r['event_type'],'gene_id':r['gene_id']} for r in support}
    by_event=defaultdict(dict); support_by=defaultdict(dict)
    for r in support:
        by_event[r['event_id']][r['sample_id']]=num(r['psi'])
        support_by[r['event_id']][r['sample_id']]=r
    summary=[]
    for typ in TYPES:
        eids=[e for e in event_meta if event_meta[e]['event_type']==typ]
        valid=[sum(by_event[e].get(s) is not None for s in SAMPLES) for e in eids]
        cells=sum(valid); summary.append({'event_type':typ,'n_defined_events':len(eids),'n_events_any_valid':sum(x>=1 for x in valid),'n_events_ge2_samples':sum(x>=2 for x in valid),'n_events_all12_valid':sum(x==12 for x in valid),'n_event_sample_cells':len(eids)*12,'n_valid_event_sample_cells':cells,'valid_event_sample_percent':f'{100*cells/(len(eids)*12):.6f}' if eids else '0'})
    write(root/'tables/new_rna/as_event_summary.tsv',summary,['event_type','n_defined_events','n_events_any_valid','n_events_ge2_samples','n_events_all12_valid','n_event_sample_cells','n_valid_event_sample_cells','valid_event_sample_percent'])
    candidates={}
    fields=['event_id','event_type','gene_id','group_1','group_2','mean_psi_group_1','mean_psi_group_2','delta_psi','n_group_1','n_group_2','new_expected_support_group1_mean','new_expected_support_group2_mean','new_expected_support_group1_min','new_expected_support_group2_min']
    for name,g1,g2 in CONTRASTS:
        out=[]
        for eid in sorted(event_meta):
            a=[by_event[eid].get(s) for s in SAMPLES if GROUP[s]==g1 and by_event[eid].get(s) is not None]; b=[by_event[eid].get(s) for s in SAMPLES if GROUP[s]==g2 and by_event[eid].get(s) is not None]
            if len(a)<2 or len(b)<2: continue
            sr1=[num(support_by[eid][s]['expected_new_support']) for s in SAMPLES if GROUP[s]==g1]; sr2=[num(support_by[eid][s]['expected_new_support']) for s in SAMPLES if GROUP[s]==g2]
            d=mean(a)-mean(b)
            if abs(d)<.1: continue
            out.append({'event_id':eid,'event_type':event_meta[eid]['event_type'],'gene_id':event_meta[eid]['gene_id'],'group_1':g1,'group_2':g2,'mean_psi_group_1':f'{mean(a):.6f}','mean_psi_group_2':f'{mean(b):.6f}','delta_psi':f'{d:.6f}','n_group_1':len(a),'n_group_2':len(b),'new_expected_support_group1_mean':f'{mean(sr1):.6f}','new_expected_support_group2_mean':f'{mean(sr2):.6f}','new_expected_support_group1_min':f'{min(sr1):.6f}','new_expected_support_group2_min':f'{min(sr2):.6f}'})
        write(root/'tables/new_rna/contrasts'/f'{name}.event_candidates.tsv',out,fields); candidates[name]=out
    return event_meta,by_event,support_by,candidates,psi

def total_and_integration(root,src,event_meta,new_by,support_by):
    total_rows=tsv(src/'tables/total_rna/event_psi.tsv'); total={r['event_id']:r for r in total_rows}
    for name,g1,g2 in CONTRASTS:
        out=[]
        for eid in sorted(event_meta):
            tv=total.get(eid,{}); ta=[num(tv.get(s)) for s in SAMPLES if GROUP[s]==g1 and num(tv.get(s)) is not None]; tb=[num(tv.get(s)) for s in SAMPLES if GROUP[s]==g2 and num(tv.get(s)) is not None]
            na=[new_by[eid].get(s) for s in SAMPLES if GROUP[s]==g1 and new_by[eid].get(s) is not None]; nb=[new_by[eid].get(s) for s in SAMPLES if GROUP[s]==g2 and new_by[eid].get(s) is not None]
            if min(len(ta),len(tb),len(na),len(nb))<2: continue
            td=mean(ta)-mean(tb); nd=mean(na)-mean(nb); tc=abs(td)>=.1; nc=abs(nd)>=.1
            cls='Shared' if tc and nc and td*nd>0 else 'Opposite' if tc and nc else 'New_enriched' if nc else 'Total_enriched' if tc else 'Stable'
            sr1=[num(support_by[eid][s]['expected_new_support']) for s in SAMPLES if GROUP[s]==g1]; sr2=[num(support_by[eid][s]['expected_new_support']) for s in SAMPLES if GROUP[s]==g2]
            out.append({'event_id':eid,'event_type':event_meta[eid]['event_type'],'gene_id':event_meta[eid]['gene_id'],'total_delta_psi':f'{td:.6f}','new_delta_psi':f'{nd:.6f}','delta_delta_psi':f'{nd-td:.6f}','classification':cls,'total_n_group1':len(ta),'total_n_group2':len(tb),'new_n_group1':len(na),'new_n_group2':len(nb),'new_expected_support_group1_mean':f'{mean(sr1):.6f}','new_expected_support_group2_mean':f'{mean(sr2):.6f}','new_expected_support_group1_min':f'{min(sr1):.6f}','new_expected_support_group2_min':f'{min(sr2):.6f}'})
        write(root/'tables/integration'/f'{name}.total_vs_new.tsv',out,['event_id','event_type','gene_id','total_delta_psi','new_delta_psi','delta_delta_psi','classification','total_n_group1','total_n_group2','new_n_group1','new_n_group2','new_expected_support_group1_mean','new_expected_support_group2_mean','new_expected_support_group1_min','new_expected_support_group2_min'])
    # Rewrite Total SUPPA delta direction while preserving p values/padj.
    for name,g1,g2 in CONTRASTS:
        p=src/'tables/total_rna/contrasts'/f'{name}.suppa_diffsplice.tsv'
        if not p.is_file(): continue
        out=[]; tv=total
        for r in tsv(p):
            eid=r['event_id']; a=[num(tv.get(eid,{}).get(s)) for s in SAMPLES if GROUP[s]==g1 and num(tv.get(eid,{}).get(s)) is not None]; b=[num(tv.get(eid,{}).get(s)) for s in SAMPLES if GROUP[s]==g2 and num(tv.get(eid,{}).get(s)) is not None]
            if len(a)>=2 and len(b)>=2: r['delta_psi']=f'{mean(a)-mean(b):.9f}'; r['group_1']=g1; r['group_2']=g2
            out.append(r)
        write(root/'tables/total_rna/contrasts'/p.name,out,list(out[0]) if out else ['event_id','event_type','gene_id','group_1','group_2','delta_psi','pvalue','padj'])
    return total

def plots(root,src,event_meta,new_by,support_by,candidates,total):
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'axes.facecolor':'white'})
    colors={'A':'#E64B35','AB':'#00A087','B':'#4DBBD5','NC':'#3C5488'}
    fnew=root/'figures/new_rna'; fint=root/'figures/integration'; fcase=root/'figures/case_studies'; fnew.mkdir(parents=True,exist_ok=True); fint.mkdir(parents=True,exist_ok=True); fcase.mkdir(parents=True,exist_ok=True)
    # coverage
    cov={(t,s):[0,0] for t in TYPES for s in SAMPLES}
    for eid,m in event_meta.items():
        for s in SAMPLES: cov[(m['event_type'],s)][1]+=1; cov[(m['event_type'],s)][0]+=new_by[eid].get(s) is not None
    mat=np.array([[100*cov[(t,s)][0]/cov[(t,s)][1] for s in SAMPLES] for t in TYPES]); fig,ax=plt.subplots(figsize=(11,4.8)); im=ax.imshow(mat,vmin=0,vmax=100,aspect='auto',cmap='Blues'); ax.set_xticks(range(12),SAMPLES,rotation=45,ha='right'); ax.set_yticks(range(7),TYPES); ax.set_xlabel('Sample'); ax.set_ylabel('Event type'); ax.set_title('New RNA event PSI coverage'); fig.colorbar(im,ax=ax,label='Valid PSI (%)'); fig.savefig(fnew/'new_rna_event_coverage.png',dpi=300,bbox_inches='tight'); plt.close(fig)
    # candidate counts in one canonical ordered dataframe
    counts=[len(candidates[n]) for n,_,_ in CONTRASTS]; fig,ax=plt.subplots(figsize=(9,4.5)); x=np.arange(6); ax.bar(x,counts,color='#00A087'); ax.set_xticks(x,[n for n,_,_ in CONTRASTS],rotation=35,ha='right'); ax.set_ylabel('Effect-size candidates'); ax.set_title('New RNA AS effect-size candidates'); [ax.text(i,v,str(v),ha='center',va='bottom') for i,v in enumerate(counts)]; fig.savefig(fnew/'new_rna_as_candidate_counts.png',dpi=300,bbox_inches='tight'); plt.close(fig)
    # candidate count by event type
    type_mat=np.array([[sum(r['event_type']==t for r in candidates[n]) for t in TYPES] for n,_,_ in CONTRASTS]); fig,ax=plt.subplots(figsize=(10,5)); im=ax.imshow(type_mat,aspect='auto',cmap='YlGnBu'); ax.set_xticks(range(7),TYPES); ax.set_yticks(range(6),[n for n,_,_ in CONTRASTS]); ax.set_xlabel('Event type'); ax.set_ylabel('Comparison'); ax.set_title('New RNA candidates by event type'); fig.colorbar(im,ax=ax,label='Candidate count'); fig.savefig(fnew/'new_rna_event_type_candidate_counts.png',dpi=300,bbox_inches='tight'); plt.close(fig)
    # new delta heatmap
    allc=[r for n,_,_ in CONTRASTS for r in candidates[n]]; allc.sort(key=lambda r:abs(float(r['delta_psi'])),reverse=True); top=allc[:30]; fig,ax=plt.subplots(figsize=(10,8)); arr=np.array([[float(next((x['delta_psi'] for x in candidates[n] if x['event_id']==r['event_id']),'0')) for n,_,_ in CONTRASTS] for r in top]); im=ax.imshow(arr,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto'); ax.set_xticks(range(6),[n for n,_,_ in CONTRASTS],rotation=35,ha='right'); ax.set_yticks(range(len(top)),[f"{r['gene_id']} | {r['event_type']}" for r in top]); ax.set_xlabel('Comparison'); ax.set_ylabel('Gene | event type'); ax.set_title('Top New RNA event-level ΔPSI effects'); fig.colorbar(im,ax=ax,label='ΔPSI (group 1 − group 2)'); fig.savefig(fnew/'new_rna_delta_psi_heatmap.png',dpi=300,bbox_inches='tight'); plt.close(fig)
    # six-panel integration scatter
    fig,axs=plt.subplots(2,3,figsize=(14,8),constrained_layout=True); clscol={'Shared':'#00A087','New_enriched':'#E64B35','Total_enriched':'#3C5488','Opposite':'#7B2CBF','Stable':'#BDBDBD'}
    integration={n:tsv(root/'tables/integration'/f'{n}.total_vs_new.tsv') for n,_,_ in CONTRASTS}
    for ax,(n,_,_) in zip(axs.flat,CONTRASTS):
        rs=integration[n]
        for cl,col in clscol.items():
            z=[r for r in rs if r['classification']==cl]; ax.scatter([float(r['total_delta_psi']) for r in z],[float(r['new_delta_psi']) for r in z],s=5,alpha=.45,c=col,label=cl)
        ax.axhline(0,color='#CBD5E1',lw=.6); ax.axvline(0,color='#CBD5E1',lw=.6); ax.axline((0,0),slope=1,color='#64748B',ls='--',lw=.6); ax.axvline(.1,color='#CBD5E1',ls=':',lw=.6); ax.axvline(-.1,color='#CBD5E1',ls=':',lw=.6); ax.axhline(.1,color='#CBD5E1',ls=':',lw=.6); ax.axhline(-.1,color='#CBD5E1',ls=':',lw=.6); ax.set_title(f'{n}\nn={len(rs):,}  New={sum(r["classification"]=="New_enriched" for r in rs):,}  Shared={sum(r["classification"]=="Shared" for r in rs):,}  Opposite={sum(r["classification"]=="Opposite" for r in rs):,}',fontsize=9); ax.set_xlabel('Total RNA ΔPSI'); ax.set_ylabel('New RNA ΔPSI')
    axs.flat[0].legend(frameon=False,fontsize=7,loc='best'); fig.savefig(fint/'total_vs_new_delta_psi_6panel.png',dpi=300,bbox_inches='tight'); plt.close(fig)
    # paired heatmap and event-type comparison
    allrs=[{**r,'contrast':n} for n,_,_ in CONTRASTS for r in integration[n]]; allrs.sort(key=lambda r:abs(float(r['delta_delta_psi'])),reverse=True); top=allrs[:30]; arr=np.array([[float(r['total_delta_psi']),float(r['new_delta_psi'])] for r in top for _ in [0]])
    by_contrast={n:{x['event_id']:x for x in integration[n]} for n,_,_ in CONTRASTS}; fig,ax=plt.subplots(figsize=(12,8)); wide=np.array([[float(by_contrast[n].get(r['event_id'],{}).get(k,'nan')) for n,_,_ in CONTRASTS for k in ('total_delta_psi','new_delta_psi')] for r in top]); im=ax.imshow(wide,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto'); ax.set_xticks(range(12),[x for n,_,_ in CONTRASTS for x in (n+' Total',n+' New')],rotation=45,ha='right'); ax.set_yticks(range(len(top)),[f"{r['gene_id']} | {r['event_type']}" for r in top]); ax.set_title('Paired Total and New RNA ΔPSI'); fig.colorbar(im,ax=ax,label='ΔPSI (group 1 − group 2)'); fig.savefig(fint/'total_new_delta_psi_heatmap.png',dpi=300,bbox_inches='tight'); plt.close(fig)
    et=np.zeros((7,2),int)
    for t_i,t in enumerate(TYPES):
        et[t_i,0]=sum(abs(float(r['total_delta_psi']))>=.1 and r['event_type']==t for r in allrs); et[t_i,1]=sum(abs(float(r['new_delta_psi']))>=.1 and r['event_type']==t for r in allrs)
    fig,ax=plt.subplots(figsize=(9,4.8)); x=np.arange(7); ax.bar(x-.18,et[:,0],.36,label='Total RNA',color='#3C5488'); ax.bar(x+.18,et[:,1],.36,label='New RNA',color='#00A087'); ax.set_xticks(x,TYPES); ax.set_xlabel('Event type'); ax.set_ylabel('Effect-size event count'); ax.set_title('Total versus New RNA event-type comparison'); ax.legend(frameon=False); fig.savefig(fint/'event_type_comparison.png',dpi=300,bbox_inches='tight'); plt.close(fig)
    # case selection: top six distinct events, strong delta-delta, adequate support
    cases=[]; seen=set()
    for r in allrs:
        if r['event_id'] in seen or min(int(r['new_n_group1']),int(r['new_n_group2']))<2: continue
        seen.add(r['event_id']); cases.append(r)
        if len(cases)>=6: break
    case_rows=[]
    for i,r in enumerate(cases,1):
        eid=r['event_id']; tv=total[eid]; ev=event_meta[eid]; inc=ev.get('inclusion_transcripts','')
        case_rows.append({**r,'case_id':f'case_{i:02d}','inclusion_transcripts':inc,'exclusion_transcripts':ev.get('exclusion_transcripts','')})
        fig,axs=plt.subplots(2,1,figsize=(10,6.2),gridspec_kw={'height_ratios':[2,1]},constrained_layout=True)
        # Each sample occupies one x-axis position; samples are not a time series.
        xs=np.arange(12); offsets={'Total RNA':-0.12,'New RNA':0.12}
        for label,values,col in (('Total RNA',[num(tv.get(s)) for s in SAMPLES],'#3C5488'),('New RNA',[new_by[eid].get(s) for s in SAMPLES],'#00A087')):
            valid=[j for j,v in enumerate(values) if v is not None]
            axs[0].scatter(xs[valid]+offsets[label],[values[j] for j in valid],s=38,color=col,alpha=.85,label=label,zorder=3)
        axs[0].set_xticks(xs,SAMPLES,rotation=45,ha='right'); axs[0].set_ylabel('PSI'); axs[0].set_ylim(0,1); axs[0].set_title(f"case_{i:02d}  {r['gene_id']} | {r['event_type']}\n{eid}"); axs[0].legend(frameon=False); axs[0].text(.99,.02,'Each point is an independent sample; no temporal ordering',transform=axs[0].transAxes,ha='right',va='bottom',fontsize=8,color='#64748B')
        sp=[num(support_by[eid][s]['expected_new_support']) for s in SAMPLES]; axs[1].bar(np.arange(12),sp,color=[colors[GROUP[s]] for s in SAMPLES]); axs[1].set_xticks(np.arange(12),SAMPLES,rotation=45,ha='right'); axs[1].set_ylabel('Expected New support'); axs[1].set_title('Read-assignment support'); fig.savefig(fcase/f"case_{i:02d}.png",dpi=300,bbox_inches='tight'); plt.close(fig)
    write(root/'tables/case_studies.tsv',case_rows,['case_id','event_id','event_type','gene_id','contrast','total_delta_psi','new_delta_psi','delta_delta_psi','classification','total_n_group1','total_n_group2','new_n_group1','new_n_group2','new_expected_support_group1_mean','new_expected_support_group2_mean','inclusion_transcripts','exclusion_transcripts'])
    return case_rows

def report(root, event_meta, candidates, case_rows):
    img_specs=[('figures/slam_qc/new_rna_fraction.png','按样本划分的 New RNA fraction'),('figures/slam_qc/tc_conversion_qc.png','T-to-C 转换质量控制'),('figures/total_rna/highres/event_candidate_counts.png','Total RNA AS 候选事件'),('figures/total_rna/highres/delta_psi_heatmap.png','Total RNA ΔPSI 热图'),('figures/new_rna/new_rna_event_coverage.png','New RNA 事件 PSI 覆盖度'),('figures/new_rna/new_rna_as_candidate_counts.png','New RNA AS 效应候选事件'),('figures/new_rna/new_rna_event_type_candidate_counts.png','按事件类型统计的 New RNA 候选'),('figures/new_rna/new_rna_delta_psi_heatmap.png','New RNA 高效应 ΔPSI 热图'),('figures/integration/total_vs_new_delta_psi_6panel.png','Total 与 New RNA ΔPSI 对照'),('figures/integration/total_new_delta_psi_heatmap.png','Total 与 New RNA 配对 ΔPSI 热图'),('figures/integration/event_type_comparison.png','Total 与 New RNA 事件类型对照')]+[(f'figures/case_studies/{r["case_id"]}.png',f'{r["case_id"]}：{r["gene_id"]} {r["event_type"]}') for r in case_rows]
    def cards_for(pred):
        return ''.join(f'<article class="figure-card"><button class="figure-expand figure-zoom" type="button" data-hires="{uri(root/p)}" aria-label="Zoom {html.escape(t)}"><img src="{uri(root/p)}" alt="{html.escape(t)}"></button><figcaption>{html.escape(t)}</figcaption></article>' for p,t in img_specs if pred(p) and (root/p).is_file())
    slam_cards=cards_for(lambda p:p.startswith('figures/slam_qc/'))
    total_cards=cards_for(lambda p:p.startswith('figures/total_rna/'))
    new_cards=cards_for(lambda p:p.startswith('figures/new_rna/'))
    integration_cards=cards_for(lambda p:p.startswith('figures/integration/'))
    header_zh={'sample_id':'样本','group':'组别','biological_replicate':'生物学重复','event_type':'事件类型','n_defined_events':'定义事件数','n_events_any_valid':'至少1个有效样本的事件数','n_events_ge2_samples':'至少2个有效样本的事件数','n_events_all12_valid':'12个样本均有效的事件数','n_valid_event_sample_cells':'有效事件-样本单元数','valid_event_sample_percent':'有效比例','contrast':'比较','n_compared':'可比较事件数','n_shared':'Shared','n_new_enriched':'New 富集','n_total_enriched':'Total 富集','n_opposite':'方向相反','n_stable':'稳定','n_effect_candidates':'效应候选数','case_id':'案例','gene_id':'基因','classification':'分类','total_delta_psi':'Total ΔPSI','new_delta_psi':'New ΔPSI','delta_delta_psi':'ΔΔPSI','new_expected_support_group1_mean':'组1预期 New 支持均值','new_expected_support_group2_mean':'组2预期 New 支持均值'}
    def tab(rows,fields,limit=20):
        rows=rows[:limit]; return '<div class="table-wrap"><table><thead><tr>'+''.join(f'<th>{html.escape(header_zh.get(x,x))}</th>' for x in fields)+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{html.escape(str(r.get(x,"NA")))}</td>' for x in fields)+'</tr>' for r in rows)+'</tbody></table></div>'
    sample_rows=[{'sample_id':s,'group':GROUP[s],'biological_replicate':s.rsplit('-',1)[1]} for s in SAMPLES]
    summary=tsv(root/'tables/new_rna/as_event_summary.tsv'); comp=[]
    for n,_,_ in CONTRASTS:
        rs=tsv(root/'tables/integration'/f'{n}.total_vs_new.tsv'); comp.append({'contrast':n,'n_compared':len(rs),'n_shared':sum(x['classification']=='Shared' for x in rs),'n_new_enriched':sum(x['classification']=='New_enriched' for x in rs),'n_total_enriched':sum(x['classification']=='Total_enriched' for x in rs),'n_opposite':sum(x['classification']=='Opposite' for x in rs),'n_stable':sum(x['classification']=='Stable' for x in rs)})
    nav=''.join(f'<a href="#{i}">{t}</a>' for i,t in [('overview','项目概览'),('slam','3 SLAM QC'),('total','Total RNA AS'),('new','New RNA AS'),('integration','Total 与 New 对照'),('cases','其他案例分析'),('methods','方法与说明'),('references','参考文献')])
    js='''<script>document.querySelectorAll(".figure-zoom").forEach(function(b){b.addEventListener("click",function(){var m=document.getElementById("zoom-modal");m.querySelector("img").src=b.dataset.hires;m.classList.add("open");m.setAttribute("aria-hidden","false");m.querySelector("button").focus()})});var modal=document.getElementById("zoom-modal");modal.querySelector("button").addEventListener("click",function(){modal.classList.remove("open");modal.setAttribute("aria-hidden","true")});modal.addEventListener("click",function(e){if(e.target===modal){modal.classList.remove("open");modal.setAttribute("aria-hidden","true")}});document.addEventListener("keydown",function(e){if(e.key==="Escape"){modal.classList.remove("open");modal.setAttribute("aria-hidden","true")}});if(location.search.includes("browser_audit")){setTimeout(function(){var a=[].slice.call(document.querySelectorAll(".figure-card img")).map(function(i){return {complete:i.complete,naturalWidth:i.naturalWidth,renderedWidth:i.getBoundingClientRect().width}});var m=document.querySelector("main"),p=document.createElement("pre");p.id="browser-audit";p.textContent=JSON.stringify({imageWidths:a,contentWidth:m.clientWidth,clientWidth:document.documentElement.clientWidth,scrollWidth:document.documentElement.scrollWidth});document.body.appendChild(p)},900)}</script>'''
    css='''*{box-sizing:border-box}body{margin:0;overflow-x:hidden;background:#f6f8fb;color:#1f2937;font:14px/1.7 Arial,"Microsoft YaHei",sans-serif}nav{position:fixed;left:0;top:0;bottom:0;width:215px;background:#fff;border-right:1px solid #d9e2ec;padding:18px 12px;display:flex;flex-direction:column;gap:4px}nav a{color:#52606d;text-decoration:none;padding:7px 9px}main{margin-left:235px;max-width:1480px;padding:20px}header,section{background:#fff;border:1px solid #d9e2ec;border-radius:10px;padding:22px;margin-bottom:18px}h1,h2{color:#1d3557}.metric{display:inline-block;background:#f0fdfa;border-left:4px solid #00a087;padding:10px;margin:5px;min-width:145px}.metric b{display:block;font-size:22px}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}.figure-card{border:1px solid #d9e2ec;border-radius:8px;padding:12px;min-width:0}.figure-card img{display:block;width:100%;max-width:100%;height:auto}.figure-zoom{display:block;border:0;background:#fff;padding:0;width:100%;cursor:zoom-in}.figure-zoom:focus{outline:3px solid #4DBBD5}.table-wrap{overflow:auto}table{width:100%;border-collapse:collapse;font-size:12px}th,td{padding:7px;border-bottom:1px solid #d9e2ec;text-align:left;vertical-align:top}th{background:#eef3f8}.zoom-modal{display:none;position:fixed;inset:0;background:#000c;align-items:center;justify-content:center;padding:4vh 4vw;z-index:10}.zoom-modal.open{display:flex}.zoom-panel{width:92vw;height:86vh;max-width:92vw;max-height:86vh;background:#fff;padding:14px;overflow:hidden;position:relative;display:flex;align-items:center;justify-content:center}.zoom-panel img{display:block;max-width:calc(92vw - 28px);max-height:calc(86vh - 28px);width:auto;height:auto;object-fit:contain}.zoom-close{position:fixed;right:20px;top:15px;width:42px;height:42px;border:0;border-radius:50%;font-size:25px;cursor:pointer}@media(max-width:760px){nav{position:sticky;width:auto;height:auto;flex-direction:row;overflow:auto}main{margin-left:0}.grid{grid-template-columns:1fr}.zoom-panel{width:94vw;height:84vh;max-width:94vw;max-height:84vh}.zoom-panel img{max-width:calc(94vw - 28px);max-height:calc(84vh - 28px)}}@media print{nav,.zoom-modal,.figure-zoom{display:none!important}.figure-card{break-inside:avoid}.figure-card img{max-height:none}}'''
    def download(rel):
        p=root/rel; return 'data:text/tab-separated-values;base64,'+base64.b64encode(p.read_bytes()).decode('ascii') if p.is_file() else '#downloads'
    logo_path=Path(__file__).resolve().parents[1]/'assets'/'epi-logo.png'
    html_text=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>长读长可变剪切分析报告</title><style>{css}</style></head><body><nav class="sidebar">{nav}</nav><main><header><img src="{uri(logo_path)}" alt="Epi logo" style="width:180px;max-height:58px;object-fit:contain"><h1>长读长可变剪切分析报告</h1><p>ISO-SLAM-seq Total RNA、New RNA 与 Total-vs-New event-level alternative splicing。</p><div class="metric">Samples<b>12</b></div><div class="metric">Groups<b>4</b></div><div class="metric">Shared AS events<b>{len(event_meta):,}</b></div><div class="metric">New valid PSI cells<b>{sum(x['valid_psi']=='TRUE' for x in tsv(root/'tables/new_rna/event_support.tsv')):,}</b></div></header><section class="report-section" id="overview"><h2>1 项目概览</h2><p>本项目包含 A、AB、B、NC 四组，每组 3 个独立生物学重复。Total RNA、New RNA 和 integration 使用同一 21,505-event universe。</p><h3>Sample design</h3>{tab(sample_rows,['sample_id','group','biological_replicate'],20)}<h3>New RNA event summary</h3>{tab(summary,['event_type','n_defined_events','n_events_any_valid','n_events_ge2_samples','n_events_all12_valid','n_valid_event_sample_cells','valid_event_sample_percent'],20)}</section><section class="report-section" id="design"><h2>2 ISO-SLAM experimental design</h2><p>New RNA evidence is derived from read-level ISO-SLAM posterior assignments. Each group has three independent biological replicates.</p></section><section class="report-section" id="slam"><h2>3 SLAM QC</h2><p>New/Old/Ambiguous classification and T&gt;C posterior are inherited from upstream ISO-SLAM outputs. Ambiguous reads are not forced into New RNA.</p><div class="grid">{slam_cards}</div></section><section class="report-section" id="total"><h2>4 Total RNA AS</h2><p>Total RNA event PSI and SUPPA2 evidence are retained. All displayed Total RNA ΔPSI values use group 1 − group 2; formal P values are preserved when available.</p><div class="grid">{total_cards}</div>{tab(comp,['contrast','n_compared','n_shared','n_new_enriched','n_total_enriched','n_opposite','n_stable'],20)}</section><section class="report-section" id="new"><h2>5 New RNA AS</h2><p>New RNA PSI is aggregated from posterior probability-weighted transcript assignments. Valid PSI requires expected New event support ≥ 10; missing PSI is not imputed to zero.</p><div class="grid">{new_cards}</div><h3>New RNA comparison summary</h3>{tab([{'contrast':n,'n_effect_candidates':len(candidates[n])} for n,_,_ in CONTRASTS],['contrast','n_effect_candidates'],20)}</section><section class="report-section" id="integration"><h2>6 Total vs New integration</h2><p>New-enriched is a descriptive effect category: New RNA reaches the ΔPSI threshold while Total RNA shows a smaller change. It is not New-specific statistical significance. ΔΔPSI = New ΔPSI − Total ΔPSI and is not an interaction P value.</p><div class="grid">{integration_cards}</div>{tab(comp,['contrast','n_compared','n_shared','n_new_enriched','n_total_enriched','n_opposite','n_stable'],20)}</section><section class="report-section" id="cases"><h2>7 Long-read case studies</h2><p>Cases were selected using effect magnitude, ΔΔPSI, support and replicate consistency rather than a single low-support maximum.</p>{tab(case_rows,['case_id','gene_id','event_type','contrast','classification','total_delta_psi','new_delta_psi','delta_delta_psi','new_expected_support_group1_mean','new_expected_support_group2_mean'],20)}<div class="grid">{''.join(f'<article class="figure-card"><button class="figure-expand figure-zoom" type="button" data-hires="{uri(root/"figures/case_studies"/(r["case_id"]+".png"))}" aria-label="Zoom {r["case_id"]}"><img src="{uri(root/"figures/case_studies"/(r["case_id"]+".png"))}" alt="{r["case_id"]}"></button></article>' for r in case_rows)}</div></section><section class="report-section" id="usage"><h2>8 Total / New transcript usage</h2><p>Transcript usage is retained as auxiliary evidence and is interpreted separately from event-level AS.</p></section><section class="report-section" id="methods"><h2>9 Methods and statistical interpretation</h2><ul><li>Event universe: 21,505 events shared by Total and New RNA.</li><li>Canonical direction: group 1 − group 2.</li><li>Valid PSI: expected New event support ≥ 10.</li><li>Effect candidate: absolute ΔPSI ≥ 0.10 with at least 2 valid biological replicates per group.</li><li>Formal New RNA inference is <b>NOT_VALIDATED_FOR_FORMAL_INFERENCE</b> because the posterior-weighted expected-count framework has not completed method-specific distribution validation.</li></ul></section><section class="report-section" id="downloads"><h2>10 Downloads</h2><p><a download="new_event_psi.tsv" href="{download('tables/new_rna/event_psi.tsv')}">New RNA event PSI</a> · <a download="new_event_support.tsv" href="{download('tables/new_rna/event_support.tsv')}">New RNA event support</a> · <a download="A_vs_NC.event_candidates.tsv" href="{download('tables/new_rna/contrasts/A_vs_NC.event_candidates.tsv')}">New candidate table</a> · <a download="A_vs_NC.total_vs_new.tsv" href="{download('tables/integration/A_vs_NC.total_vs_new.tsv')}">Integration table</a> · <a download="event_definition.tsv" href="{download('provenance/event_definition.tsv')}">Event definitions</a></p></section><section class="report-section" id="references"><h2>References</h2><ol><li>Alamancos GP et al. SUPPA2: fast, accurate and uncertainty-aware differential splicing analysis.</li><li>Minoche AE et al. Modelling T-to-C conversions in SLAM-seq.</li><li>Voldman A et al. IsoQuant: transcript discovery and quantification from long reads.</li><li>Van den Berge K et al. DRIMSeq and stageR for differential transcript usage analysis.</li></ol></section><footer class="footer">Offline self-contained ISO-SLAM AS delivery · Figure text is English</footer></main><div id="zoom-modal" class="zoom-modal" aria-hidden="true"><div class="zoom-panel"><button class="zoom-close" type="button" aria-label="Close enlarged image">×</button><img alt="Enlarged figure"></div></div>{js}</body></html>'''
    translations={
        '长读长可变剪切分析报告':'ISO-SLAM-seq可变剪切分析报告',
        'ISO-SLAM-seq Total RNA、New RNA 与 Total-vs-New event-level alternative splicing。':'本报告分析 ISO-SLAM-seq 的 Total RNA、New RNA 及二者的事件级可变剪切对照。',
        'Samples':'样本数','Groups':'组数','Shared AS events':'共享 AS 事件数','New valid PSI cells':'New 有效 PSI 单元数',
        'Sample design':'样本设计','New RNA event summary':'New RNA 事件汇总','ISO-SLAM experimental design':'ISO-SLAM 实验设计',
        'New RNA evidence is derived from read-level ISO-SLAM posterior assignments. Each group has three independent biological replicates.':'New RNA 证据来自 ISO-SLAM read-level posterior 分配；每组包含 3 个独立生物学重复。',
        'New/Old/Ambiguous classification and T&gt;C posterior are inherited from upstream ISO-SLAM outputs. Ambiguous reads are not forced into New RNA.':'New、Old、Ambiguous 分类及 T&gt;C posterior 继承自上游 ISO-SLAM 输出；Ambiguous reads 不被强行归入 New RNA。',
        'Total RNA event PSI and SUPPA2 evidence are retained. All displayed Total RNA ΔPSI values use group 1 − group 2; formal P values are preserved when available.':'保留 Total RNA 事件 PSI 和 SUPPA2 证据；报告中的 Total RNA ΔPSI 统一为组1 − 组2，正式 P 值在可用时保留。',
        'New RNA PSI is aggregated from posterior-weighted transcript assignments. Valid PSI requires expected New event support ≥ 10; missing PSI is not imputed to zero.':'New RNA PSI 由 posterior 加权的转录本分配汇总得到；有效 PSI 要求预期 New event support ≥ 10，缺失 PSI 不填充为 0。',
        'New-enriched is a descriptive effect category: New RNA reaches the ΔPSI threshold while Total RNA shows a smaller change. It is not New-specific statistical significance. ΔΔPSI = New ΔPSI − Total ΔPSI and is not an interaction P value.':'New 富集是描述性效应类别：New RNA 达到 ΔPSI 阈值而 Total RNA 变化较小；这不是 New 特异性统计显著性。ΔΔPSI = New ΔPSI − Total ΔPSI，也不是 interaction P 值。',
        'Cases were selected using effect magnitude, ΔΔPSI, support and replicate consistency rather than a single low-support maximum.':'案例依据效应大小、ΔΔPSI、support 和重复一致性选择，而不是依据单个低支持度最大值选择。',
        'Transcript usage is retained as auxiliary evidence and is interpreted separately from event-level AS.':'转录本使用结果作为辅助证据保留，并与事件级 AS 分开解释。',
        'Event universe: 21,505 events shared by Total and New RNA.':'事件集合：Total RNA 与 New RNA 共享 21,505 个事件。','Canonical direction: group 1 − group 2.':'方向约定：组1 − 组2。','Valid PSI: expected New event support ≥ 10.':'有效 PSI：预期 New event support ≥ 10。','Effect candidate: absolute ΔPSI ≥ 0.10 with at least 2 valid biological replicates per group.':'效应候选：绝对 ΔPSI ≥ 0.10，且每组至少 2 个有效生物学重复。',
        'Formal New RNA inference is':'正式 New RNA 推断状态为','because the posterior-weighted expected-count framework has not completed method-specific distribution validation.':'，因为 posterior 加权预期计数框架尚未完成方法特异性的分布验证。','NOT_VALIDATED_FOR_FORMAL_INFERENCE':'未验证正式统计推断',
        'New RNA event PSI':'New RNA 事件 PSI','New RNA event support':'New RNA 事件 support','New candidate table':'New RNA 候选表','Integration table':'整合表','Event definitions':'事件定义',
        'References':'参考文献','Downloads':'下载','Long-read case studies':'长读长案例分析','Total / New transcript usage':'Total / New 转录本使用','Methods and statistical interpretation':'方法与统计说明',
        'Click to enlarge; press Escape to close.':'点击放大；按 Escape 关闭。','Zoom ':'放大：','Offline self-contained ISO-SLAM AS delivery · Figure text is English':'离线自包含 ISO-SLAM-seq 报告 · 图中文字为英文',
        'Total RNA AS':'Total RNA 可变剪切','New RNA AS':'New RNA 可变剪切','Total vs New integration':'Total 与 New RNA 整合','New RNA evidence is derived from read-level ISO-SLAM posterior assignments.':'New RNA 证据来自 ISO-SLAM read-level posterior 分配。',
        'Total RNA、New RNA 和 integration 使用同一 21,505-event universe。':'Total RNA、New RNA 和整合分析使用同一 21,505 个事件集合。','New RNA PSI is aggregated from posterior-weighted transcript assignments.':'New RNA PSI 由 posterior 加权的转录本分配汇总得到。','Total RNA event PSI and SUPPA2 evidence are retained.':'保留 Total RNA 事件 PSI 和 SUPPA2 证据。','New-enriched is a descriptive effect category':'New 富集是描述性效应类别',
        '7 Long-read case studies':'7 其他案例分析','8 Total / New transcript usage':'8 Total / New 转录本使用','9 Methods and statistical interpretation':'9 方法与统计说明','10 Downloads':'10 下载'
    }
    for old,new in translations.items(): html_text=html_text.replace(old,new)
    html_text=re.sub(r'<h3>样本设计</h3><div class="table-wrap">.*?</div>','',html_text,flags=re.S)
    html_text=html_text.replace('点击放大；按 Escape 关闭。','')
    for section_id in ('design','usage','downloads'):
        html_text=re.sub(r'<section class="report-section" id="'+section_id+r'">.*?</section>','',html_text,flags=re.S)
    html_text=html_text.replace('离线自包含 ISO-SLAM-seq 报告 · 图中文字为英文','ISO-SLAM-seq可变剪切分析报告')
    (root/'report').mkdir(parents=True,exist_ok=True); (root/'report/final_ISO_SLAM_AS_report.html').write_text(html_text,encoding='utf-8')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',type=Path,required=True); ap.add_argument('--samples',type=Path,required=True); ap.add_argument('--results-dir',type=Path,required=True); ap.add_argument('--work-dir',type=Path,required=True); a=ap.parse_args(); c=config(a.config); root=a.results_dir; src=Path(c['source_rebuild_root']); root.mkdir(parents=True,exist_ok=True); a.work_dir.mkdir(parents=True,exist_ok=True)
    global SAMPLES, GROUP
    SAMPLES, GROUP = sample_sheet(a.samples)
    if any(p for p in root.iterdir() if p.name not in {'tables','figures','reports','deliverables','editable','qc'}): raise RuntimeError(f'non-empty output: {root}')
    copy_source(src,root); make_metadata(root)
    support=tsv(src/'tables/new_rna/event_support.tsv'); event_meta,new_by,support_by,candidates,newpsi=new_tables(root,src)
    defs=event_definitions(Path(c['event_ioe_dir']),event_meta,root); event_meta.update(defs); total=total_and_integration(root,src,event_meta,new_by,support_by)
    slam_src=Path(c['slam_figures_root']); slam_dst=root/'figures/slam_qc'; slam_dst.mkdir(parents=True,exist_ok=True)
    for name in ('new_rna_fraction.png','tc_conversion_qc.png'):
        if (slam_src/name).is_file(): shutil.copy2(slam_src/name,slam_dst/name)
    source_delivery=c.get('source_delivery_root','').strip()
    if source_delivery and Path(source_delivery).is_dir(): restyle_slam_fraction(Path(source_delivery),slam_dst)
    cases=plots(root,src,event_meta,new_by,support_by,candidates,total)
    report(root,event_meta,candidates,cases)
    (root/'provenance/new_rna_psi_filtering.md').write_text('# New RNA PSI filtering\n\nValid PSI requires expected New event support >= 10; values below threshold remain NA.\n',encoding='utf-8')
    (root/'provenance/contrast_direction.md').write_text('# Contrast direction\n\nAll contrasts use group 1 - group 2: A-NC, B-NC, AB-NC, A-B, AB-A and AB-B.\n',encoding='utf-8')
    print(root)
if __name__=='__main__': main()
