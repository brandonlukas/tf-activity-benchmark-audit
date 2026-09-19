"""Controls for propagate.py: shuffled-graph null, and own-mRNA augmentation as the trivial competitor."""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/controls.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd
from propagate import adata, S, tfs, A, propagate, evaluate, fr0   # reuses loaded data, reruns the grid print

def shuffle_preserving_degree(A, n_swaps=20000, seed=0):
    """Edge swaps keep every TF's in/out degree; destroys who-regulates-whom.
    Duplicate check runs against the LIVE edge set: checking the original A let swaps
    create duplicates that collapsed on fill (7296 edges -> 6825, degrees off by up to 39)."""
    rng = np.random.default_rng(seed)
    e = np.argwhere(A.values != 0); w = A.values[e[:, 0], e[:, 1]]
    live = {(int(a), int(b)) for a, b in e}
    for _ in range(n_swaps):
        i, j = rng.integers(len(e), size=2)
        a, b = map(int, e[i]); c, d = map(int, e[j])
        if a == d or c == b or (a, d) in live or (c, b) in live: continue
        live.difference_update({(a, b), (c, d)}); live.update({(a, d), (c, b)})
        e[i], e[j] = (a, d), (c, b)
    B = np.zeros_like(A.values); B[e[:, 0], e[:, 1]] = w
    return pd.DataFrame(B, index=A.index, columns=A.columns)

own = adata.to_df().reindex(columns=tfs).fillna(0)    # TF's own logFC, 0 if not measured
own = own / own.values.std()
print(f"\n{'variant':32s} top10  top5  better worse")
_rows = []                                   # F-029: collected for ../results/controls.csv, prints unchanged
def show(name, S2):
    t10, t5, fr = evaluate(S2)
    _rows.append(dict(variant=name, seeds=1, top10=t10, top5=t5, top10_min=t10, top10_max=t10,
                      better=np.mean(fr < fr0), worse=np.mean(fr > fr0)))
    print(f"{name:32s} {t10:.2f}   {t5:.2f}  {np.mean(fr<fr0):.2f}   {np.mean(fr>fr0):.2f}")
show("ulm", S)
# V-006: signs travel with the source's edge slot, so signed OUT-degree is preserved, signed IN-degree is not
Ashuf = [shuffle_preserving_degree(A, seed=s) for s in range(10)]   # T-001: 10 draws, not one seed-0 draw
def show_multi(name, Slist):
    """mean [min-max] over seeds of the shuffled null."""
    r = np.array([evaluate(S2)[:2] for S2 in Slist])
    _rows.append(dict(variant=name, seeds=len(Slist), top10=r[:,0].mean(), top5=r[:,1].mean(),
                      top10_min=r[:,0].min(), top10_max=r[:,0].max(), better=np.nan, worse=np.nan))
    print(f"{name:32s} {r[:,0].mean():.2f}   {r[:,1].mean():.2f}  "
          f"[{r[:,0].min():.2f}-{r[:,0].max():.2f}] [{r[:,1].min():.2f}-{r[:,1].max():.2f}]")
for norm, alpha in [("sum", .5), ("mean", 1)]:
    show(f"up {norm} a={alpha}", pd.DataFrame(propagate(S, A, alpha, "up", norm), index=S.index, columns=tfs))
    show_multi(f"up {norm} a={alpha} SHUF x10 mean", [pd.DataFrame(propagate(S, B, alpha, "up", norm), index=S.index, columns=tfs) for B in Ashuf])
for a in [0.5, 1, 2]:
    show(f"ulm + {a}*own mRNA", S + a * own)
show("ulm + 1*own + up mean a=1", pd.DataFrame(propagate(S + own, A, 1, "up", "mean"), index=S.index, columns=tfs))

pd.DataFrame(_rows).to_csv("../results/controls.csv", index=False)   # F-029: the printed table, persisted

E = A.values != 0
assert all(((B.values != 0).sum() == E.sum() and ((B.values != 0).sum(1) == E.sum(1)).all()
            and ((B.values != 0).sum(0) == E.sum(0)).all()) for B in Ashuf), "shuffle is not degree-preserving"
