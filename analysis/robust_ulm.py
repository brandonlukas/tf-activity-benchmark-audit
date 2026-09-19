"""ULM vs Huber-IRLS ULM on knockTF. Same regression, robust loss; residual weights = 'edge usage'."""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/robust_ulm.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

from concurrent.futures import ThreadPoolExecutor
import anndata as ad, decoupler as dc, numpy as np, pandas as pd

adata = ad.read_h5ad("../data/knocktf.h5ad")
net = pd.read_parquet("../data/collectri.parquet")
dc.mt.ulm(adata, net, tmin=5)
S = adata.obsm["score_ulm"]; tfs = S.columns
genes = adata.var_names
W = net[net.source.isin(tfs)].pivot_table(index="target", columns="source", values="weight", aggfunc="first").reindex(index=genes, columns=tfs).fillna(0).values
Y = np.asarray(adata.X, dtype=float)
Y = np.nan_to_num(Y)

def huber_ulm(y, W, k=1.345, iters=10):
    """Per-TF slope t-stat, Huber IRLS. k in units of MAD-scaled residuals. Starts at OLS (=ULM)."""
    yc = y - y.mean(); n = len(y)
    Wc = W - W.mean(0)
    beta = (Wc * yc[:, None]).sum(0) / (Wc**2).sum(0)               # OLS = ULM slope
    w = np.ones_like(W)
    for _ in range(iters):
        R = yc[:, None] - Wc * beta                                 # genes x TFs residuals
        s = 1.4826 * np.median(np.abs(R), 0) + 1e-12
        w = np.minimum(1, k * s / (np.abs(R) + 1e-12))              # Huber weights
        beta = (w * Wc * yc[:, None]).sum(0) / (w * Wc**2).sum(0)
    R = yc[:, None] - Wc * beta
    sigma2 = (w * R**2).sum(0) / (w.sum(0) - 2)
    se = np.sqrt(sigma2 / (w * Wc**2).sum(0))
    return beta / se, w

def huber_all(Y, workers=24, cells=None):
    """huber_ulm over every experiment (row of Y). Rows are independent and numpy releases the GIL, so threads do.
    cells[i] = (gene, TF) index pair or None; when given, also return that row's FINAL Huber weight at that
    cell (T-025: the weight on the perturbed TF's own gene in its own regression, V-008's 0.094)."""
    # ponytail: memory-bandwidth bound, ~7x on 24 cores, same as processes. Threads because a process pool started
    # while this module is being imported deadlocks on the import lock when it pickles the worker function.
    def one(i):
        t, w = huber_ulm(Y[i], W)
        return t, np.nan if cells is None or cells[i] is None else w[cells[i]]
    with ThreadPoolExecutor(workers) as ex:
        T, wc = zip(*ex.map(one, range(len(Y))))
    return np.vstack(T) if cells is None else (np.vstack(T), np.array(wc))

gi = {g: i for i, g in enumerate(genes)}; ti = {t: i for i, t in enumerate(tfs)}
own_cell = [(gi[t], ti[t]) if t in gi and t in ti else None for t in adata.obs.source]
T, w_own = huber_all(Y, cells=own_cell)        # w_own: Huber weight on the knockdown's own gene, unmasked
S2 = pd.DataFrame(T, index=S.index, columns=tfs)
assert np.corrcoef(S.values.ravel(), S2.values.ravel())[0, 1] > 0.5   # same thing in the bulk case

def evaluate(S, name):
    fr, wrong = [], []
    for exp, tf in adata.obs.source.items():
        if tf not in S.columns: continue
        sign = adata.obs.loc[exp, "type_p"]; s = sign * S.loc[exp]
        fr.append(s.rank(ascending=False)[tf] / S.shape[1]); wrong.append(s[tf] < 0)
    fr = np.array(fr)
    print(f"{name:12s} n={len(fr)} top10%={np.mean(fr<=.1):.2f} top5%={np.mean(fr<=.05):.2f} "
          f"median_frac={np.median(fr):.3f} wrong_sign={np.mean(wrong):.2f}")
    return fr
f1 = evaluate(S, "ulm"); f2 = evaluate(S2, "huber_ulm")
print("paired: huber better in", np.mean(f2 < f1).round(2), "worse in", np.mean(f2 > f1).round(2))
# F-029: persist. One row per scoreable experiment, both arms' rank fractions and the Huber
# weight this script's whole argument turns on. Written after the last print, changes none.
_sc = [(e, t) for e, t in adata.obs.source.items() if t in S.columns]
pd.DataFrame(dict(exp=[e for e, _ in _sc], tf=[t for _, t in _sc], ulm_frac=f1, huber_frac=f2,
                  huber_weight_own_gene=w_own[[i for i, t in enumerate(adata.obs.source)
                                               if t in S.columns]])
             ).to_csv("../results/robust_ulm.csv", index=False)
