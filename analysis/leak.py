"""Does propagation still help once the perturbed TF's own mRNA is masked? (rules out own-mRNA leaking via upstream regulons)"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/leak.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd
from propagate import A, propagate, tfs
from controls import shuffle_preserving_degree, Ashuf

adata = ad.read_h5ad("../data/knocktf.h5ad"); net = pd.read_parquet("../data/collectri.parquet")
X = adata.to_df()
for exp, tf in adata.obs.source.items():
    if tf in X.columns: X.loc[exp, tf] = 0                      # erase the knockdown's own mRNA drop
masked = ad.AnnData(X, obs=adata.obs)
dc.mt.ulm(masked, net, tmin=5); S = masked.obsm["score_ulm"].reindex(columns=tfs)

def evaluate(S):
    fr = []
    for exp, tf in adata.obs.source.items():
        if tf not in S.columns: continue
        s = adata.obs.loc[exp, "type_p"] * S.loc[exp]; fr.append(s.rank(ascending=False)[tf] / S.shape[1])
    return np.array(fr)
fr0 = evaluate(S)
print(f"\n{'variant (own mRNA masked)':32s} top10  top5  better worse")
_rows = []                                   # F-029: collected for ../results/leak.csv, prints unchanged
def show(name, S2):
    fr = evaluate(S2); _rows.append(dict(variant=name, seeds=1, top10=np.mean(fr<=.1), top5=np.mean(fr<=.05),
        top10_min=np.mean(fr<=.1), top10_max=np.mean(fr<=.1), better=np.mean(fr<fr0), worse=np.mean(fr>fr0)))
    print(f"{name:32s} {np.mean(fr<=.1):.2f}   {np.mean(fr<=.05):.2f}  {np.mean(fr<fr0):.2f}   {np.mean(fr>fr0):.2f}")
show("ulm", S)
def show_multi(name, Slist):
    r = np.array([[np.mean(evaluate(S2) <= .1), np.mean(evaluate(S2) <= .05)] for S2 in Slist])
    _rows.append(dict(variant=name, seeds=len(Slist), top10=r[:,0].mean(), top5=r[:,1].mean(),
                      top10_min=r[:,0].min(), top10_max=r[:,0].max(), better=np.nan, worse=np.nan))
    print(f"{name:32s} {r[:,0].mean():.2f}   {r[:,1].mean():.2f}  "
          f"[{r[:,0].min():.2f}-{r[:,0].max():.2f}] [{r[:,1].min():.2f}-{r[:,1].max():.2f}]")
for norm, alpha in [("sum", .5), ("mean", 1)]:
    show(f"up {norm} a={alpha}", pd.DataFrame(propagate(S, A, alpha, "up", norm), index=S.index, columns=tfs))
    show_multi(f"up {norm} a={alpha} SHUF x10 mean", [pd.DataFrame(propagate(S, B, alpha, "up", norm), index=S.index, columns=tfs) for B in Ashuf])

pd.DataFrame(_rows).to_csv("../results/leak.csv", index=False)   # F-029: variants + per-experiment ranks
pd.DataFrame(dict(exp=[e for e, t in adata.obs.source.items() if t in S.columns],
                  tf=[t for t in adata.obs.source if t in S.columns],
                  masked_frac=fr0)).to_csv("../results/leak_rows.csv", index=False)

r, c = np.array([(i, X.columns.get_loc(t)) for i, t in enumerate(adata.obs.source) if t in X.columns]).T
assert (X.values[r, c] == 0).all() and (adata.to_df().values[r, c] != 0).all(), "perturbed TF's own gene not masked"
