"""Model-agnostic evaluation of frozen predictions against observed responses.

This module does not train, run, or infer the semantics of any prediction model.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .core import cosine_matrix, score_rows, summarize_scores
from .pipeline import _make, _save, _sha, _utc

def _ids(value, base, name, allow_empty=False):
    ids=json.loads((base/value).read_text()) if isinstance(value,str) else value
    if not isinstance(ids,list) or (not ids and not allow_empty) or any(not isinstance(i,str) or not i for i in ids):
        raise ValueError(name+' requires a nonempty list of nonempty string IDs')
    if len(ids)!=len(set(ids)):raise ValueError(name+' has duplicate IDs')
    return ids

def _load(spec,base,name):
    for key in ['values','targets','features','output_space','readout','provenance']:
        if key not in spec or (not spec[key] and key!='targets'):raise ValueError(name+' requires '+key)
    for key in ['output_space','readout']:
        if not isinstance(spec[key],str):raise ValueError(name+' '+key+' must be a nonempty string')
    targets=_ids(spec['targets'],base,name+' targets',allow_empty=name=='predictions');features=_ids(spec['features'],base,name+' features')
    x=np.asarray(np.load(base/spec['values'],allow_pickle=False),dtype=np.float64)
    if x.shape!=(len(targets),len(features)):raise ValueError(name+' array/axis mismatch')
    return x,targets,features

def score_predictions(config_path,output):
    """Validate exact coordinate identity, retain missing predictions, save scores.

    Observed rows must cover the exact fixed gallery. Prediction rows can be a
    subset, with missing target rows explicitly unavailable rather than dropped.
    No target or coordinate intersection is inferred from data values.
    """
    cp=Path(config_path);base=cp.parent;cfg=json.loads(cp.read_text())
    for key in ['dataset_id','gallery_targets','predictions','observed','evaluation_provenance']:
        if key not in cfg or not cfg[key]:raise ValueError('config requires '+key)
    gallery=_ids(cfg['gallery_targets'],base,'fixed gallery')
    if len(gallery)<2:raise ValueError('fixed gallery needs at least two targets')
    p,pt,pf=_load(cfg['predictions'],base,'predictions');o,ot,of=_load(cfg['observed'],base,'observed')
    for key in ['output_space','readout']:
        if cfg['predictions'][key]!=cfg['observed'][key]:raise ValueError('incompatible '+key)
    if set(pf)!=set(of):raise ValueError('feature axes differ; intersection is not allowed')
    if set(ot)!=set(gallery):raise ValueError('observed targets must exactly cover fixed gallery')
    if not set(pt).issubset(gallery):raise ValueError('prediction targets outside fixed gallery')
    om={v:i for i,v in enumerate(ot)};fm={v:i for i,v in enumerate(of)};pm={v:i for i,v in enumerate(pt)}
    observed=o[[om[t] for t in gallery]][:,[fm[v] for v in pf]]
    prediction=np.full_like(observed,np.nan)
    for i,t in enumerate(gallery):
        if t in pm:prediction[i]=p[pm[t]]
    s=cosine_matrix(prediction,observed);r=score_rows(s);summary=summarize_scores(s)
    coverage=[]
    for i,t in enumerate(gallery):
        if t not in pm:status='missing_prediction'
        elif not np.isfinite(prediction[i]).all():status='nonfinite_prediction'
        elif np.linalg.norm(prediction[i])==0:status='zero_prediction_direction'
        elif not np.isfinite(s[i]).all():status='incomplete_observed_gallery'
        else:status='evaluable'
        observed_status='evaluable' if np.isfinite(observed[i]).all() and np.linalg.norm(observed[i])>0 else ('zero_observed_direction' if np.isfinite(observed[i]).all() else 'nonfinite_observed')
        coverage.append(dict(target_id=t,prediction_present=t in pm,prediction_status=status,observed_status=observed_status,**{key:values[i] for key,values in r.items()}))
    out=_make(output)
    np.savez_compressed(out/'score_matrix.npz',scores=s,query_targets=np.asarray(gallery),gallery_targets=np.asarray(gallery),feature_ids=np.asarray(pf))
    pd.DataFrame(coverage).to_csv(out/'per_target.csv',index=False)
    hashes={'config':_sha(cp)}
    for kind in ['predictions','observed']:
        for key in ['values','targets','features']:
            if isinstance(cfg[kind][key],str):hashes[kind+'.'+key]=_sha(base/cfg[kind][key])
    if isinstance(cfg['gallery_targets'],str):hashes['gallery_targets']=_sha(base/cfg['gallery_targets'])
    _save(out/'input_provenance.json',dict(configuration=cfg,input_sha256=hashes,observed_row_reordered=ot!=gallery,observed_features_reordered=of!=pf,prediction_row_reordered=pt!=gallery,coordinate_intersection_used=False))
    clean={key:(None if isinstance(v,float) and not np.isfinite(v) else v) for key,v in summary.items()}
    report=dict(status='completed',dataset=cfg['dataset_id'],task='frozen_prediction_to_observed_response_retrieval',output_space=cfg['predictions']['output_space'],readout=cfg['predictions']['readout'],summary=clean,predictions_supplied=len(pt),fixed_gallery_size=len(gallery),missing_prediction_targets=[t for t in gallery if t not in pm],scope='Evaluation of supplied frozen arrays only; no model inference, fitting, causal diagnosis, or biological replication implied.',uncertainty='One supplied frozen comparison; no inferential P value or resampling confidence interval computed.',finished_utc=_utc())
    _save(out/'report.json',report)
    rank=clean['mean_correct_target_rank']
    (out/'report.md').write_text('\n'.join(['# Frozen prediction evaluation','',report['scope'],'',f"Dataset: {cfg['dataset_id']}",f"Output space/readout: {report['output_space']} / {report['readout']}",f"Fixed candidates: {len(gallery)}; supplied predictions: {len(pt)}; fully evaluable query rows: {summary['n_finite_query_rows']}.",f"Mean normalized correct-target rank: {rank if rank is not None else 'unavailable'}.",'','Rank reference under random labels is 0.5; rank is not percent accuracy. Zero directions and incomplete candidate-score rows remain unavailable. Coverage must accompany any available-case aggregate. Inputs are explicitly aligned by unique target and feature IDs; no coordinate intersection is taken.','',report['uncertainty']])+'\n')
    return report
