"""ULM + one-step propagation along CollecTRI TF->TF edges. alpha=0 is plain ULM."""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/propagate.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd

adata = ad.read_h5ad("../data/knocktf.h5ad")
net = pd.read_parquet("../data/collectri.parquet")
dc.mt.ulm(adata, net, tmin=5)
S = adata.obsm["score_ulm"]; tfs = S.columns
tt = net[net.source.isin(tfs) & net.target.isin(tfs)]
A = tt.pivot_table(index="source", columns="target", values="weight", aggfunc="first").reindex(index=tfs, columns=tfs).fillna(0)
A = A.mask(np.eye(len(tfs), dtype=bool), 0)        # drop autoregulation
deg = ((A != 0) | (A.T != 0)).sum(1)

def propagate(S, A, alpha, direction, norm):
    """S' = S + alpha * S @ P.  P[b, a] = signed edge from a's neighbour b -> contribution to a."""
    P = {"down": A.T, "up": A, "both": A.T + A}[direction]   # column a receives from rows b
    if norm == "mean": P = P / np.maximum((P != 0).sum(0), 1)
    return S + alpha * (S.values @ P.values)

def evaluate(S):
    fr = []
    for exp, tf in adata.obs.source.items():
        if tf not in S.columns: continue
        s = adata.obs.loc[exp, "type_p"] * S.loc[exp]
        fr.append(s.rank(ascending=False)[tf] / S.shape[1])
    fr = np.array(fr); return np.mean(fr <= .1), np.mean(fr <= .05), fr

base10, base5, fr0 = evaluate(S)
print(f"ulm            top10%={base10:.2f} top5%={base5:.2f}")
rows = []
for direction in ["down", "up", "both"]:
    for norm in ["sum", "mean"]:
        for alpha in [0.25, 0.5, 1, 2]:
            S2 = pd.DataFrame(propagate(S, A, alpha, direction, norm), index=S.index, columns=tfs)
            t10, t5, fr = evaluate(S2)
            hub = np.corrcoef(np.abs(S2).mean(0), deg)[0, 1]
            rows.append(dict(direction=direction, norm=norm, alpha=alpha, top10=t10, top5=t5,
                             better=np.mean(fr < fr0), worse=np.mean(fr > fr0), hub_corr=hub))
res = pd.DataFrame(rows).sort_values("top10", ascending=False)
print(f"ulm hub_corr (mean |score| vs TF-TF degree) = {np.corrcoef(np.abs(S).mean(0), deg)[0,1]:.2f}")
print(res.round(2).to_string(index=False))
res.to_csv("../results/propagate.csv", index=False)          # F-029: the printed grid, persisted
assert (res.alpha > 0).all() and res.top10.between(0, 1).all()
