"""Portable counts -> normalized bundle -> balanced evidence report.

No project paths, known accession IDs, targets or model names are hardcoded here.
"""
from pathlib import Path
import json,hashlib,datetime,gzip
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from .core import seeded_rng,assert_disjoint,response_weights,cosine_from_grams,summarize_scores
from .bundle import validate_bundle
from .metadata import validate_cell_metadata


def _save(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def _sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def _utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def _make(path):
 p=Path(path)
 if p.exists():raise FileExistsError('refuse to overwrite '+str(p))
 p.mkdir(parents=True);return p


def prepare_counts(config_path,output):
    """Create a source bundle from explicit raw-count and metadata files.

    If only a coordinate subset is supplied, normalization must use the supplied
    original full-feature totals; the subset's row sum is never silently substituted.
    """
    cp=Path(config_path);cfg=json.loads(cp.read_text());base=cp.parent
    cells=validate_cell_metadata(pd.read_csv(base/cfg['cells']));features=pd.read_csv(base/cfg['features'])
    x=sparse.load_npz(base/cfg['counts']).tocsr().astype(np.float64)
    if x.shape!=(len(cells),len(features)):raise ValueError('count/axis shape mismatch')
    if not cells.cell_id.is_unique or not features.feature_id.is_unique:raise ValueError('duplicate axes')
    if not np.isfinite(x.data).all()or(x.data<0).any()or not np.allclose(x.data,np.rint(x.data)):raise ValueError('input must be finite nonnegative integer counts')
    total=cells[cfg['full_count_total_column']].to_numpy(dtype=float)
    if not np.isfinite(total).all()or(total<=0).any()or(total+1e-8<np.asarray(x.sum(1)).ravel()).any():raise ValueError('invalid original full-count totals')
    mask=~features.gene_symbol.astype(str).isin(cfg['source_panel']).to_numpy()
    x=x.multiply((10000/total)[:,None]).tocsr();x.data=np.log1p(x.data);x=x[:,mask].tocsr();features=features.loc[mask].reset_index(drop=True)
    guide_min=int(cfg.get('min_source_cells_per_guide',20));sample_min=int(cfg.get('min_sampled_cells_per_guide',32))
    tg={};ledger=[]
    for target in cfg['target_ids']:
        groups={}
        sub=cells[(cells.target_id==target)&(~cells.is_control)]
        for guide,g in sub.groupby('guide_id',sort=True):
            ok=len(g)>=sample_min and g.source_guide_n.min()>=guide_min and not g.source_excluded_guide.any()
            if ok:groups[str(guide)]=g.index.astype(int).tolist()
        eligible=len(groups)>=2;ledger.append({'target':target,'qualified_guides':len(groups),'eligible':eligible})
        if eligible:tg[target]=groups
    if len(tg)<2:raise ValueError('fewer than two qualified target identities')
    reserve=int(cfg.get('pca_fit_cells_per_control_guide',20));ref_n=int(cfg.get('reference_cells_per_control_guide',16));ctrl={};fit=[]
    for q,g in cells[cells.is_control].groupby('guide_id',sort=True):
        if g.source_excluded_guide.any() or g.source_guide_n.min()<guide_min or len(g)<reserve+2*ref_n:continue
        ix=seeded_rng(cfg['dataset_id'],'pca_reserve',q).permutation(g.index).astype(int).tolist();fit+=ix[:reserve];ctrl[str(q)]=ix[reserve:]
    if len(ctrl)<4:raise ValueError('four qualified control guides required after source_excluded_guide, min_source_cells_per_guide and cell-budget checks; resolve source metadata or supply sufficient eligible controls')
    assert_disjoint(fit,*ctrl.values(),*[ids for gg in tg.values()for ids in gg.values()])
    out=_make(output);cells.to_csv(out/'cells.csv',index=False);features.to_csv(out/'features.csv',index=False);sparse.save_npz(out/'rna.npz',x)
    reps={'RNA':{'file':'rna.npz','type':'sparse_rna','dimensions':x.shape[1]}}
    if cfg.get('pca_dimensions',64)>0:
        dim=min(int(cfg['pca_dimensions']),len(fit)-1,x.shape[1]);pca=PCA(n_components=dim,svd_solver='full',whiten=False).fit(x[fit].toarray())
        (out/'representations').mkdir();proj=np.asarray(x@pca.components_.T)-pca.mean_@pca.components_.T;np.save(out/'representations/PCA.npy',proj)
        np.savez_compressed(out/'pca_fit.npz',components=pca.components_,mean=pca.mean_,fit_indices=fit)
        reps['PCA']={'file':'representations/PCA.npy','type':'array','dimensions':dim,'fit_indices':fit,'training':'independent control-only'}
    _save(out/'design.json',{'targets':list(tg),'target_groups':tg,'control_groups':ctrl,'pca_fit_indices':fit,'panel_excluded':cfg['source_panel'],'min_source_cells_per_guide':guide_min,'reference_guides_per_arm':2,'reference_cells_per_guide':ref_n,'seed':2026091601,'iterations':100,'treated_budgets':[16,8,12],'dataset_id':cfg['dataset_id']})
    pd.DataFrame(ledger).to_csv(out/'eligibility.csv',index=False)
    _save(out/'manifest.json',{'schema_version':'1.0','experiment':cfg['dataset_id'],'cell_count':len(cells),'rna_features':x.shape[1],'representations':reps,'source_provenance':cfg['provenance'],'coordinate_scope':cfg['coordinate_scope'],'normalization_total_column':cfg['full_count_total_column'],'input_hashes':{k:_sha(base/cfg[k])for k in ['counts','cells','features']},'config_sha256':_sha(cp),'prepared_utc':_utc()})
    validation=validate_bundle(out);_save(out/'validation.json',validation);return validation


def evaluate_measured(bundle_path,output,iterations=100,n=16):
    """Balanced measurement/representation evaluation; no native prediction implied."""
    if not 1<=n<=16 or iterations<1:raise ValueError('n must be 1..16 and iterations positive')
    b=Path(bundle_path);validation=validate_bundle(b);m=json.loads((b/'manifest.json').read_text());d=json.loads((b/'design.json').read_text());out=_make(output)
    ref_n=int(d['reference_cells_per_guide']);targets=d['targets'];draws=[]
    for it in range(iterations):
        rng=seeded_rng(m['experiment'],'matched',it);qs=rng.permutation(sorted(d['control_groups'])).tolist();ctrl={q:rng.permutation(d['control_groups'][q]).tolist()for q in qs[:4]}
        ca={q:ctrl[q][:ref_n]for q in qs[:2]};cb={q:ctrl[q][ref_n:2*ref_n]for q in qs[:2]};cd={q:ctrl[q][:ref_n]for q in qs[2:4]}
        assert_disjoint(d['pca_fit_indices'],sum(ca.values(),[]),sum(cb.values(),[]),sum(cd.values(),[]))
        ga={};tech={};cross={}
        for t in targets:
            qa,qb=rng.choice(sorted(d['target_groups'][t]),2,replace=False).tolist();ia=rng.permutation(d['target_groups'][t][qa]).tolist();ib=rng.permutation(d['target_groups'][t][qb]).tolist()
            if len(ia)<32 or len(ib)<32:raise ValueError('fixed n16-compatible guide universe required')
            ga[t]={qa:ia[:n]};tech[t]={qa:ia[16:16+n]};cross[t]={qb:ib[:n]};assert_disjoint(ia,ib)
        draws.append({'iteration':it,'query':ga,'technical':tech,'cross_guide':cross,'control_a':ca,'references':{'same_identity_disjoint_cells':cb,'shared_cells':ca,'different_guide_identities':cd}})
    with gzip.open(out/'draws.json.gz','wt')as f:json.dump(draws,f)
    rows=[]
    for name,spec in m['representations'].items():
        x=sparse.load_npz(b/spec['file'])if spec['type']=='sparse_rna'else np.asarray(np.load(b/spec['file']),dtype=float);g=x@x.T;g=g.toarray()if sparse.issparse(g)else g
        matrices={}
        for task in ['technical','cross_guide']:
            for arm in draws[0]['references']:matrices[(task,arm)]=[]
        for draw in draws:
            wa=response_weights(m['cell_count'],draw['query'],draw['control_a'])
            for task in ['technical','cross_guide']:
                for arm,ctrl in draw['references'].items():
                    wb=response_weights(m['cell_count'],draw[task],ctrl);s=cosine_from_grams(wa,wb,g);matrices[(task,arm)].append(s)
                    for direction,ss in [('A_to_B',s),('B_to_A',s.T)]:rows.append(dict(representation=name,task=task,reference=arm,n=n,iteration=draw['iteration'],direction=direction,**summarize_scores(ss)))
        for (task,arm),ss in matrices.items():np.savez_compressed(out/f'{name}_{task}_{arm}.npz',targets=np.array(targets),scores=np.array(ss))
    frame=pd.DataFrame(rows);frame.to_csv(out/'iteration_metrics.csv',index=False)
    metrics=['mean_correct_target_rank','top1_fraction_tie_adjusted','median_correct_cosine','median_correct_minus_other_cosine'];average=frame.groupby(['representation','task','reference','iteration'],as_index=False)[metrics].mean();summary=average.groupby(['representation','task','reference'],as_index=False)[metrics].mean();summary.to_csv(out/'summary.csv',index=False)
    report={'status':'completed','dataset':m['experiment'],'targets':targets,'n_per_treated_arm':n,'controls_per_arm':2*ref_n,'iterations':iterations,'validation':validation,'all_query_rows_finite':bool(frame.complete.all()),'scope':'measured response and supplied representation retrieval only; no native predictor executed','uncertainty':'repeated cell draws are technical conditional samples, not biological replication','bundle_manifest_sha256':_sha(b/'manifest.json'),'finished_utc':_utc()};_save(out/'report.json',report)
    lines=['# Layered measured-response report','',report['scope'],'',f"Dataset: {m['experiment']}; K={len(targets)}; n={n} per treated arm; {2*ref_n} control cells per arm; {iterations} iterations.",'','| Representation | Task | Reference | Mean normalized rank | Mean correct cosine |','|---|---|---|---:|---:|']
    for r in summary.itertuples():lines.append(f'| {r.representation} | {r.task} | {r.reference} | {r.mean_correct_target_rank:.4f} | {r.median_correct_cosine:.4f} |')
    lines+=['','Random-label rank reference is 0.5; ranks are not accuracy percentages. Missing rows remain missing. This report makes no biological root-cause diagnosis or claim of clinical relevance. Technical splits are not donor replication.']
    (out/'report.md').write_text('\n'.join(lines)+'\n');return report
