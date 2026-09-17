"""Portable source-bundle validation, independent of dataset names and local project paths."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from .core import assert_disjoint
from .metadata import validate_cell_metadata


def validate_bundle(path, check_array_values=True):
    """Validate explicit axes, source labels, control exclusion and frozen guide groups.

    Bundles contain manifest.json, design.json, cells.csv, features.csv, and relative
    representation paths. RNA may be sparse; other representations are dense arrays.
    This function validates data identity and design, not biological correctness.
    """
    path=Path(path);m=json.loads((path/'manifest.json').read_text());d=json.loads((path/'design.json').read_text())
    cells=validate_cell_metadata(pd.read_csv(path/'cells.csv'));features=pd.read_csv(path/'features.csv');n=len(cells)
    if not features.feature_id.is_unique or features.feature_id.isna().any():raise ValueError('feature IDs must be unique')
    if m['cell_count']!=n:raise ValueError('manifest cell count mismatch')
    if list(d['target_groups'])!=d['targets']:raise ValueError('target axis mismatch')
    fit=list(map(int,d['pca_fit_indices']));reference=[];treated=[]
    def check_indices(ix):
        if any(i<0 or i>=n for i in ix):raise ValueError('out-of-bounds cell index')
        if len(ix)!=len(set(ix)):raise ValueError('duplicate group cell')
    for q,ix in d['control_groups'].items():
        check_indices(ix);reference.extend(ix)
        if not (cells.guide_id.iloc[ix]==q).all()or not cells.is_control.iloc[ix].all():raise ValueError('control identity mismatch')
    for g,groups in d['target_groups'].items():
        for q,ix in groups.items():
            check_indices(ix);treated.extend(ix)
            if not (cells.guide_id.iloc[ix]==q).all()or not(cells.target_id.iloc[ix]==g).all()or cells.is_control.iloc[ix].any():raise ValueError('treated identity mismatch')
    check_indices(fit)
    if not cells.is_control.iloc[fit].all():raise ValueError('PCA fit includes treated cells')
    assert_disjoint(fit,reference,treated)
    # Legacy bundles did not save the setting and used the documented default 20.
    guide_min=int(d.get('min_source_cells_per_guide',20))
    for role,ix in [('PCA fit',fit),('reference',reference),('treated',treated)]:
        assigned=cells.iloc[ix]
        excluded=assigned.loc[assigned.source_excluded_guide,'guide_id'].unique().tolist()
        if excluded:raise ValueError(f'{role} includes source_excluded_guide=True guide(s) {excluded!r}; rebuild the design without these guides')
        low=assigned.loc[assigned.source_guide_n<guide_min,'guide_id'].unique().tolist()
        if low:raise ValueError(f'{role} includes guide(s) {low!r} below min_source_cells_per_guide={guide_min}; rebuild the design with eligible guides')
    reps={}
    for name,spec in m['representations'].items():
        rel=Path(spec['file'])
        if rel.is_absolute()or '..'in rel.parts:raise ValueError('representation path must be contained relative path')
        x=sparse.load_npz(path/rel)if spec['type']=='sparse_rna'else np.load(path/rel,mmap_mode='r')
        if x.ndim!=2 or x.shape[0]!=n or x.shape[1]!=spec['dimensions']:raise ValueError('representation shape mismatch')
        if spec['type']=='sparse_rna'and x.shape[1]!=len(features):raise ValueError('RNA feature axis mismatch')
        values=x.data if sparse.issparse(x)else x
        if check_array_values and not np.isfinite(values).all():raise ValueError('nonfinite representation')
        reps[name]={'shape':list(x.shape),'finite_checked':check_array_values}
    return {'status':'passed','cell_count':n,'target_count':len(d['targets']),'fit_controls':len(fit),'reference_pool_cells':len(reference),'treated_pool_cells':len(treated),'source_exclusion_checked':True,'source_min_cells_checked':guide_min,'source_min_cells_origin':'design' if 'min_source_cells_per_guide' in d else 'legacy_default_20','representations':reps,'scope':'identity, source eligibility and split validation; not biological correctness'}
