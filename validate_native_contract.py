"""Portable metadata validator for an explicitly budgeted native-response task.

Standard-library only. It checks identities and declared computational inputs,
not biological validity of a token operation or hidden model training exposure.
"""
from pathlib import Path
import argparse,csv,json,hashlib,datetime

def ids(values,label):
    if not isinstance(values,list) or any(not isinstance(x,str) or not x for x in values):
        raise ValueError(label+': expected string ID list')
    if len(set(values))!=len(values): raise ValueError(label+': duplicate IDs')
    return set(values)

def validate(cfg,cells,sequences=None):
    for k in ['source_id','source_id_column','special_token_ids','gallery','native_groups','reference_groups','pca_fit_cells','native_reserve_cells','budget']:
        if k not in cfg: raise ValueError('missing '+k)
    if not isinstance(cfg['source_id'],str) or not cfg['source_id'] or not isinstance(cfg['source_id_column'],str) or not cfg['source_id_column']:raise ValueError('source identity and metadata column must be nonempty strings')
    if not cells or any(r.get(cfg['source_id_column'])!=cfg['source_id'] for r in cells):raise ValueError('source ID does not match metadata source identity')
    if not isinstance(cfg['special_token_ids'],list) or any(isinstance(t,bool) or not isinstance(t,int) for t in cfg['special_token_ids']):raise ValueError('special token IDs must be an explicit integer list')
    index={r['cell_id']:r for r in cells}
    if len(index)!=len(cells):raise ValueError('duplicate metadata cell IDs')
    if any(r['is_control'] not in ('True','False','true','false','1','0') for r in cells):
        raise ValueError('is_control must be an explicit boolean')
    control={r['cell_id'] for r in cells if r['is_control'] in ('True','true','1')}
    gallery=ids(cfg['gallery'],'gallery')
    if len(gallery)<2:raise ValueError('need at least two gallery targets')
    if set(cfg['native_groups'])!=gallery:raise ValueError('native target axis must exactly equal declared gallery; no silent intersection')
    for k in ['native_guides','native_cells_per_guide','reference_min_cells_per_guide','reference_min_guides']:
        if k not in cfg['budget'] or isinstance(cfg['budget'][k],bool) or not isinstance(cfg['budget'][k],int) or cfg['budget'][k]<1:raise ValueError('invalid budget: '+k)
    reserve=ids(cfg['native_reserve_cells'],'native reserve');fit=ids(cfg['pca_fit_cells'],'PCA fit');reference=set()
    for guide,group in cfg['reference_groups'].items():
        g=ids(group,'reference '+guide)
        if len(g)<cfg['budget']['reference_min_cells_per_guide']:raise ValueError('reference guide capacity below contract')
        if g&reference:raise ValueError('reference cells duplicated across guides')
        for cell in g:
            if cell not in index or index[cell]['guide_id']!=guide:raise ValueError('reference ID/guide mismatch')
        reference|=g
    if len(cfg['reference_groups'])<cfg['budget']['reference_min_guides']:raise ValueError('too few reference guides')
    if not (reserve|fit|reference)<=control:raise ValueError('role pool contains unknown or treated cell')
    if reserve&fit or reserve&reference or fit&reference:raise ValueError('native/reference/PCA role pools overlap')
    selected=set();coverage=[]
    for target in cfg['gallery']:
        groups=cfg['native_groups'][target]
        if len(groups)!=cfg['budget']['native_guides']:raise ValueError('native guide budget mismatch for '+target)
        allcells=set()
        for guide,group in groups.items():
            g=ids(group,'native '+target+'/'+guide)
            if len(g)!=cfg['budget']['native_cells_per_guide']:raise ValueError('native cell budget mismatch for '+target)
            if allcells&g:raise ValueError('native cell counted in multiple guides')
            if not g<=reserve:raise ValueError('native start outside reserved role pool')
            for cell in g:
                if index[cell]['guide_id']!=guide:raise ValueError('native ID/guide mismatch')
            allcells|=g
        if sequences is not None:
            if 'target_tokens' not in cfg or target not in cfg['target_tokens']:raise ValueError('missing target token identity')
            token=cfg['target_tokens'][target]
            if isinstance(token,bool) or not isinstance(token,int) or token<0 or token in cfg['special_token_ids']:raise ValueError('target token must be a nonnegative gene token, not a special token or boolean')
            if any(cell not in sequences or sequences[cell].count(token)!=1 for cell in allcells):raise ValueError('target token not uniquely visible in native start')
        selected|=allcells
        coverage.append({'target':target,'native_guides':len(groups),'native_cells':len(allcells)})
    return {'status':'passed','source':cfg['source_id'],'source_metadata_identity_checked':True,'fixed_gallery_size':len(gallery),'unique_native_starts':len(selected),'native_start_uses':sum(x['native_cells'] for x in coverage),'native_reserve_cells':len(reserve),'reference_cells':len(reference),'pca_fit_cells':len(fit),'role_pools_pairwise_disjoint':True,'declared_token_visibility_checked':sequences is not None,'target_to_biological_gene_mapping_validated':False,'coverage':coverage,'scope':'Declared cell/guide/budget/token-visibility contract only; target-to-biological-gene mapping needs its separate fixed dictionary audit. No inference or biological validity claim.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('config');p.add_argument('--output',required=True);a=p.parse_args();cp=Path(a.config);cfg=json.loads(cp.read_text());base=cp.parent
    with (base/cfg['cells_csv']).open() as f:cells=list(csv.DictReader(f))
    seq=json.loads((base/cfg['cell_sequences_json']).read_text()) if cfg.get('cell_sequences_json') else None
    report=validate(cfg,cells,seq);report['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();report['input_sha256']={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [cp,base/cfg['cells_csv']]+([base/cfg['cell_sequences_json']] if seq is not None else [])}
    out=Path(a.output)
    if out.exists():raise FileExistsError('refuse to overwrite existing audit')
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
