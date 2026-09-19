"""T-004 / F-010: by which path does the perturbed TF's own mRNA reach the scores?

ULM is univariate, so gene X's logFC can enter TF X's OWN score only through a CollecTRI
self-edge X->X. It enters every OTHER TF's score whenever X sits in that TF's regulon (the
path F-004 says propagation exploits). Zeroing the gene (leak.py) closes both at once.
Here the two are separated on the network side, which leaves the gene in the matrix:

  ulm              full network, unmasked matrix                       (F-001: 0.45)
  masked           full network, own gene zeroed in its own row        (F-005: 0.34)
  no self-edge     network minus the edge X->X only; gene still visible to every other regulon
  no other regulon network minus every edge Y->X (Y != X); X->X kept
  gene off network both of the above at once (network-side analogue of masked)

Each condition is a per-experiment network: ULM is univariate, so dropping X->X changes only
column X, and dropping gene X from other regulons changes only the other columns. Both are done
by zeroing entries of decoupler's own adjacency matrix and calling decoupler's own ULM kernel
(dc.mt.ulm.func), so no ULM is reimplemented; the final assert checks the kernel reproduces
propagate.py's score matrix exactly on the untouched network.

PREDICTIONS (stated before running, from the task):
  (a) in the no-self-edge group, masking moves ULM top-10% by under 2 points;
  (b) removing the self-edge alone reproduces most of the 0.45->0.34 drop: predict 0.35;
  (c) for up/sum/alpha=0.5 propagation (0.57 unmasked -> 0.28 masked), removing the self-edge
      alone does NOT reproduce the collapse (predict it stays near 0.50), while removing the
      gene from the other TFs' regulons does.
If (a) fails, the competitor-score path matters for plain ULM too, not only for propagation.
No gain is claimed here, so no shuffled network: the masked condition is the control.
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/selfedge.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd, scipy.sparse as sps
from propagate import adata, S as S_ulm, A, propagate, tfs   # reruns the propagation grid
import leak                                                  # reruns controls.py + the masked table

net = pd.read_parquet("../data/collectri.parquet")

# decoupler's own preprocessing, so the adjacency we edit is the one it would have scored
mat, obs, var = dc.pp.extract(adata)
if sps.issparse(mat): mat = mat.toarray()
pnet = dc.pp.prune(features=var, net=net, tmin=5)
sources, feats, adjm = dc.pp.adjmat(features=var, net=pnet)
adjm_orig = adjm.copy()                       # T-019: the final assert checks the whole adjacency is restored
src_i = {s: i for i, s in enumerate(sources)}
var_i = {v: i for i, v in enumerate(var)}
obs_i = {o: i for i, o in enumerate(obs)}

def score(m, adj, rows=None):
    return dc.mt.ulm.func(m if rows is None else m[rows], adj)[0]

es_full = pd.DataFrame(score(mat, adjm), index=obs, columns=sources).reindex(columns=tfs)

# --- condition: every self-edge removed. ULM is univariate, so this changes column Y only via Y->Y
self_rc = [(var_i[s], src_i[s]) for s in sources if s in var_i and adjm[var_i[s], src_i[s]] != 0]
adjm_ns = adjm.copy()
for r, c in self_rc: adjm_ns[r, c] = 0
es_ns = pd.DataFrame(score(mat, adjm_ns), index=obs, columns=sources).reindex(columns=tfs)

exps = [(e, tf) for e, tf in adata.obs.source.items() if tf in tfs]
has_self = np.array([tf in var_i and adjm[var_i[tf], src_i[tf]] != 0 for _, tf in exps])
in_mat = np.array([tf in var_i for _, tf in exps])

# --- condition: gene X dropped from every OTHER TF's regulon (one network per perturbed TF)
S_noother, S_nogene = es_full.copy(), es_full.copy()
for tf in sorted({t for _, t in exps}):
    rows = [obs_i[e] for e, t in exps if t == tf]
    if tf not in var_i: continue
    g, keep = var_i[tf], adjm[var_i[tf]].copy()
    adjm[g] = 0; adjm[g, src_i[tf]] = keep[src_i[tf]]          # keep the self-edge
    r = pd.DataFrame(score(mat, adjm, rows), index=obs[rows], columns=sources).reindex(columns=tfs)
    adjm[g] = keep                                             # restore
    S_noother.loc[r.index] = r
    S_nogene.loc[r.index] = r
S_noself = es_full.copy()
for e, tf in exps:                          # only the perturbed TF's own column loses its self-edge
    S_noself.loc[e, tf] = es_ns.loc[e, tf]
    S_nogene.loc[e, tf] = es_ns.loc[e, tf]  # no other regulon AND no self-edge

def evaluate(S):                            # copied from propagate.py / leak.py
    fr = []
    for exp, tf in adata.obs.source.items():
        if tf not in S.columns: continue
        s = adata.obs.loc[exp, "type_p"] * S.loc[exp]
        fr.append(s.rank(ascending=False)[tf] / S.shape[1])
    return np.array(fr)

cond = {"ulm (full net)": es_full, "no self-edge": S_noself, "no other regulon": S_noother,
        "gene off network": S_nogene, "masked (leak.py)": leak.S.reindex(columns=tfs)}
fr = {k: evaluate(v) for k, v in cond.items()}
fr0 = fr["ulm (full net)"]
prop = {k: evaluate(pd.DataFrame(propagate(v, A, .5, "up", "sum"), index=v.index, columns=tfs))
        for k, v in cond.items()}

def line(name, f, ref, n=None):
    n = len(f) if n is None else n
    print(f"{name:34s} {np.mean(f<=.1):.3f} {np.mean(f<=.1)*len(f):5.0f}/{len(f):<4d} "
          f"{np.mean(f<=.05):.3f}  {np.mean(f<ref):.2f}   {np.mean(f>ref):.2f}")

print(f"\n{len(exps)} experiments, {len({t for _, t in exps})} TFs; "
      f"{has_self.sum()} experiments ({np.mean(has_self):.0%}) whose perturbed TF has a CollecTRI "
      f"self-edge surviving tmin=5, {(~has_self).sum()} without; "
      f"{(~in_mat).sum()} whose gene is not in the scored matrix at all")
print(f"self-edges in the pruned network: {len(self_rc)} of {len(sources)} TFs")

print("\n(1) plain ULM unmasked vs own-gene-masked, split by self-edge "
      "(better/worse = masked vs unmasked ULM, within group)")
print(f"{'group':34s} top10   hits/n    top5  better worse")
for lab, m in [("self-edge: unmasked", has_self), ("self-edge: MASKED", has_self),
               ("no self-edge: unmasked", ~has_self), ("no self-edge: MASKED", ~has_self)]:
    f = fr["masked (leak.py)"][m] if "MASKED" in lab else fr0[m]
    line(lab, f, fr0[m])

print("\n(2) decomposition of the 0.45 -> 0.34 drop, all experiments "
      "(better/worse vs unmasked full-network ULM)")
print(f"{'variant':34s} top10   hits/n    top5  better worse")
for k, f in fr.items(): line(k, f, fr0)
print("    same, restricted to the self-edge group")
for k, f in fr.items(): line("  " + k, f[has_self], fr0[has_self])
print("    same, restricted to the NO-self-edge group")
for k, f in fr.items(): line("  " + k, f[~has_self], fr0[~has_self])

print("\n(3) same decomposition under up/sum/alpha=0.5 propagation "
      "(0.57 unmasked -> 0.28 masked in F-004; better/worse vs unmasked plain ULM)")
print(f"{'variant + propagation':34s} top10   hits/n    top5  better worse")
for k, f in prop.items(): line(k, f, fr0)
print("    same, restricted to the self-edge group")
for k, f in prop.items(): line("  " + k, f[has_self], fr0[has_self])
print("    same, restricted to the NO-self-edge group")
for k, f in prop.items(): line("  " + k, f[~has_self], fr0[~has_self])

# F-029: persist. Per-experiment rank fractions for all five conditions plain and propagated,
# plus the group table the three blocks above print. No printed line changes.
_rows = pd.DataFrame(dict(exp=[e for e, _ in exps], tf=[t for _, t in exps], has_self_edge=has_self,
                          in_matrix=in_mat))
for k in cond:
    _rows[k.split(" (")[0].replace(" ", "_")] = fr[k]
    _rows["prop_" + k.split(" (")[0].replace(" ", "_")] = prop[k]
_rows.to_csv("../results/selfedge_rows.csv", index=False)
pd.DataFrame([dict(block=blk, group=g, variant=k, n=int(msk.sum()),
                   top10=float(np.mean(d[k][msk] <= .1)), hits=int(np.sum(d[k][msk] <= .1)),
                   top5=float(np.mean(d[k][msk] <= .05)),
                   better=float(np.mean(d[k][msk] < fr0[msk])), worse=float(np.mean(d[k][msk] > fr0[msk])))
              for blk, d in [("ulm", fr), ("propagated up/sum/a=0.5", prop)]
              for g, msk in [("all", np.ones(len(fr0), bool)), ("self-edge", has_self),
                             ("no self-edge", ~has_self)]
              for k in d]).to_csv("../results/selfedge.csv", index=False)

assert np.abs(es_full.values - S_ulm.reindex(columns=tfs).values).max() < 1e-8 \
    and np.array_equal(adjm, adjm_orig), \
    "ULM kernel does not reproduce propagate.py's scores, or the adjacency was left edited"
