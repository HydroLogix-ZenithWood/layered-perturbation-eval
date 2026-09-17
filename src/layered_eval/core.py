"""Identity-safe, missing-aware kernels for layered perturbation evaluation."""
from __future__ import annotations
import hashlib
import numpy as np
from scipy import sparse


def seeded_rng(*keys, seed=2026091601):
    digest=hashlib.sha256(('/'.join(map(str,(seed,*keys)))).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], 'little'))


def align_rows(values, input_ids, required_ids):
    """Align explicitly by unique IDs; never assume matching array row order."""
    a=list(map(str,input_ids));b=list(map(str,required_ids))
    if len(a)!=len(set(a)) or len(b)!=len(set(b)): raise ValueError('duplicate IDs')
    if len(values)!=len(a): raise ValueError('array/ID length mismatch')
    look={v:i for i,v in enumerate(a)}
    missing=[x for x in b if x not in look]
    if missing: raise ValueError(f'missing IDs: {missing[:5]}')
    return values[[look[x] for x in b]]


def assert_disjoint(*groups):
    seen=set()
    for group in groups:
        g=list(group)
        if len(g)!=len(set(g)): raise ValueError('duplicate within split')
        overlap=seen.intersection(g)
        if overlap: raise ValueError(f'split/reference leakage: {sorted(overlap)[:5]}')
        seen.update(g)


def cosine_matrix(x,y):
    x=np.asarray(x,dtype=np.float64);y=np.asarray(y,dtype=np.float64)
    if x.ndim!=2 or y.ndim!=2 or x.shape[1]!=y.shape[1]:raise ValueError('incompatible coordinate axes')
    xn=np.linalg.norm(x,axis=1);yn=np.linalg.norm(y,axis=1)
    with np.errstate(invalid='ignore',divide='ignore'):
        s=(x@y.T)/xn[:,None]/yn[None,:]
    s[(~np.isfinite(x).all(1))|(xn==0),:]=np.nan
    s[:,(~np.isfinite(y).all(1))|(yn==0)]=np.nan
    return np.clip(s,-1,1)


def score_rows(scores, correct_columns=None):
    s=np.asarray(scores,dtype=np.float64)
    if s.ndim!=2 or s.shape[1]<2: raise ValueError('need >=2 gallery candidates')
    if correct_columns is None:
        if s.shape[0]!=s.shape[1]:raise ValueError('supply correct columns for non-square scores')
        correct_columns=np.arange(len(s))
    corr=np.asarray(correct_columns,dtype=int)
    if len(corr)!=len(s) or (corr<0).any() or (corr>=s.shape[1]).any():raise ValueError('invalid target mapping')
    k=s.shape[1];rank=np.full(len(s),np.nan);top=rank.copy();cos=rank.copy();margin=rank.copy()
    for i in np.flatnonzero(np.isfinite(s).all(1)):
        v=s[i,corr[i]];others=np.delete(s[i],corr[i]);rank[i]=(np.sum(others<v)+.5*np.sum(others==v))/(k-1)
        top[i]=1/np.sum(s[i]==v) if v==np.max(s[i]) else 0
        cos[i]=v;margin[i]=v-others.mean()
    return dict(rank=rank,top1=top,correct_cosine=cos,correct_minus_other=margin)


def summarize_scores(scores, correct_columns=None):
    r=score_rows(scores,correct_columns);finite=np.isfinite(r['rank']);n=int(finite.sum())
    def avg(key,median=False):
        return float((np.median if median else np.mean)(r[key][finite])) if n else float('nan')
    return dict(n_targets=len(r['rank']),n_candidates=np.asarray(scores).shape[1],n_finite_query_rows=n,
                complete=bool(finite.all()),mean_correct_target_rank=avg('rank'),
                top1_fraction_tie_adjusted=avg('top1'),median_correct_cosine=avg('correct_cosine',True),
                median_correct_minus_other_cosine=avg('correct_minus_other',True),
                chance_mean_rank=.5,chance_top1=1/np.asarray(scores).shape[1])


def aggregate_responses(values, target_groups, control_groups):
    """Guide-equal target minus guide-equal reference, preserving target insertion order."""
    def means(groups):
        if not groups or any(not len(v) for v in groups.values()):raise ValueError('empty guide')
        return np.mean([np.asarray(values[v].mean(axis=0)).ravel() for v in groups.values()],axis=0)
    control=means(control_groups)
    return np.asarray([means(groups)-control for groups in target_groups.values()])


def response_weights(n_cells, groups, controls):
    """Sparse guide-equal aggregation design; each row sums to zero."""
    rr=[];cc=[];vv=[]
    if not controls or any(not len(ix) for ix in controls.values()):raise ValueError('empty controls')
    for row, tg in enumerate(groups.values()):
        if not tg or any(not len(ix) for ix in tg.values()):raise ValueError('empty target guide')
        for q,ix in tg.items():
            rr.extend([row]*len(ix));cc.extend(ix);vv.extend([1/(len(tg)*len(ix))]*len(ix))
        for q,ix in controls.items():
            rr.extend([row]*len(ix));cc.extend(ix);vv.extend([-1/(len(controls)*len(ix))]*len(ix))
    w=sparse.csr_matrix((vv,(rr,cc)),shape=(len(groups),n_cells),dtype=np.float64)
    if not np.allclose(np.asarray(w.sum(1)).ravel(),0,atol=1e-12):raise ValueError('unbalanced response')
    return w


def cosine_from_grams(wa,wb,gaa,gbb=None,gab=None):
    """Equivalent to cosine(wa @ X, wb @ Y), avoiding repeated high-dimensional work."""
    if gbb is None:gbb=gaa
    if gab is None:gab=gaa
    ag=wa@gaa;bg=wb@gbb
    an=np.asarray(wa.multiply(ag).sum(1)).ravel();bn=np.asarray(wb.multiply(bg).sum(1)).ravel()
    numerator=np.asarray((wb@(wa@gab).T).T)
    with np.errstate(invalid='ignore',divide='ignore'):
        s=numerator/np.sqrt(np.maximum(an,0))[:,None]/np.sqrt(np.maximum(bn,0))[None,:]
    s[an<=1e-12,:]=np.nan;s[:,bn<=1e-12]=np.nan
    return np.clip(s,-1,1)


def conditional_label_null(scores,n=1000,seed=2026091601):
    """Finite, fixed-gallery label calibration; invalid rows may not yield a P value."""
    s=np.asarray(scores,dtype=np.float64)
    if s.shape[0]!=s.shape[1] or len(s)<2 or not np.isfinite(s).all():raise ValueError('label null needs fully finite square matrix')
    rng=np.random.default_rng(seed);obs=summarize_scores(s)['mean_correct_target_rank']
    values=np.asarray([summarize_scores(s,rng.permutation(len(s)))['mean_correct_target_rank']for _ in range(n)])
    return dict(observed=obs,p=float((1+np.sum(values>=obs-1e-12))/(n+1)),null=values)
