"""T-003: is the 57->28 collapse specific to masking the perturbed TF, or does masking ANY gene do it?

Zero one random gene per experiment instead of the perturbed TF's own gene, then redo ULM and
up/sum/alpha=0.5 propagation. Two variants: (a) any measured gene, (b) a measured gene matched on
how many CollecTRI regulons it belongs to (the leak path in F-004 is regulon membership, so an
unmatched random gene could be too weak a control on its own).

PREDICTION (stated before running, T-003): the collapse does NOT reproduce. Random-masked
top-10% stays at about 0.45 for ULM and about 0.57 for up/sum/a=0.5, against 0.34/0.28 when the
perturbed TF itself is masked. Matched masking may cost a point or two more than unmatched but
nothing like 29 points. Shuffled-graph null stays near 0.31.
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/randmask.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd
from propagate import A, propagate, tfs          # reruns the propagation grid
from controls import shuffle_preserving_degree, Ashuf   # reruns unmasked table (incl. 10-seed shuffle)
import leak                                      # reruns the perturbed-TF-masked table

adata = ad.read_h5ad("../data/knocktf.h5ad"); net = pd.read_parquet("../data/collectri.parquet")
X0 = adata.to_df()
memb = net.groupby("target").size().reindex(X0.columns).fillna(0)   # regulons each measured gene sits in

def evaluate(S):                                  # same statistic as propagate.py / leak.py
    fr = []
    for exp, tf in adata.obs.source.items():
        if tf not in S.columns: continue
        s = adata.obs.loc[exp, "type_p"] * S.loc[exp]; fr.append(s.rank(ascending=False)[tf] / S.shape[1])
    return np.array(fr)

def pick(tf, rng, matched):
    """A measured gene that is not the perturbed TF; matched => similar regulon membership."""
    if not matched: return rng.choice(memb.index[memb.index != tf])
    c = memb.drop(tf, errors="ignore").sub(memb.get(tf, 0)).abs().nsmallest(50).index
    return rng.choice(c)

rows, picks = [], []
for matched in [False, True]:
    for seed in range(5):
        rng = np.random.default_rng(seed)
        X = X0.copy()
        for exp, tf in adata.obs.source.items():
            g = pick(tf, rng, matched); picks.append((tf, g)); X.loc[exp, g] = 0
        a = ad.AnnData(X, obs=adata.obs); dc.mt.ulm(a, net, tmin=5)
        S = a.obsm["score_ulm"].reindex(columns=tfs)
        fr0 = evaluate(S)
        P = pd.DataFrame(propagate(S, A, .5, "up", "sum"), index=S.index, columns=tfs)
        frp = evaluate(P)
        shuf = [evaluate(pd.DataFrame(propagate(S, B, .5, "up", "sum"), index=S.index, columns=tfs)) for B in Ashuf]
        tag = "matched" if matched else "any"
        rows.append(dict(variant=f"{tag} ulm", seed=seed, top10=np.mean(fr0 <= .1), top5=np.mean(fr0 <= .05),
                         better=0.0, worse=0.0))
        rows.append(dict(variant=f"{tag} up sum a=0.5", seed=seed, top10=np.mean(frp <= .1), top5=np.mean(frp <= .05),
                         better=np.mean(frp < fr0), worse=np.mean(frp > fr0)))
        rows.append(dict(variant=f"{tag} up sum a=0.5 SHUF", seed=seed,
                         top10=np.mean([np.mean(f <= .1) for f in shuf]), top5=np.mean([np.mean(f <= .05) for f in shuf]),
                         better=np.mean([np.mean(f < fr0) for f in shuf]), worse=np.mean([np.mean(f > fr0) for f in shuf])))

df = pd.DataFrame(rows)
g = df.groupby("variant", sort=False).agg(top10=("top10", "mean"), t10_lo=("top10", "min"), t10_hi=("top10", "max"),
                                          top5=("top5", "mean"), better=("better", "mean"), worse=("worse", "mean"))
print("\nrandom-gene mask, 5 seeds (better/worse are paired vs the SAME-seed random-masked ulm)")
print(g.round(3).to_string())
print("\nreference rows above: unmasked ulm 0.45 / up sum a=0.5 0.57 ; perturbed-TF-masked ulm 0.34 / 0.28")
mm = pd.DataFrame(picks[len(picks)//2:], columns=["tf", "gene"])   # matched half
print(f"matched picks: mean |regulons(pick) - regulons(TF gene)| = "
      f"{(memb.reindex(mm.gene).values - memb.reindex(mm.tf).fillna(0).values).__abs__().mean():.2f}, "
      f"median regulons(TF gene) = {memb.reindex(mm.tf).fillna(0).median():.0f}")
g.to_csv("../results/randmask.csv")                 # F-029: the printed 5-seed table
df.to_csv("../results/randmask_seeds.csv", index=False)   # per-seed rows behind it

assert (pd.DataFrame(picks, columns=["tf", "gene"]).tf != pd.DataFrame(picks, columns=["tf", "gene"]).gene).all() \
    and (X.values[np.arange(len(X)), [X.columns.get_loc(g) for _, g in picks[-len(X):]]] == 0).all(), \
    "masked gene must be a non-perturbed measured gene, actually zeroed"
