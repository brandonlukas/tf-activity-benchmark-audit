"""Huber-IRLS ULM vs plain ULM with the perturbed TF's own gene masked (T-002, from V-003).

V-003 left F-003 ("Huber is 13 points worse than ULM") ambiguous: both arms were unmasked, so part of
Huber's deficit could be the robust loss clipping the own-mRNA abundance leak rather than losing activity
signal. This runs the same comparison with the leak.py mask (X.loc[exp, tf] = 0) applied before scoring.

Prediction, stated before the run:
  masked ULM   top-10% = 0.34  (known from leak.py / F-005)
  masked Huber top-10% = 0.26  (i.e. the deficit shrinks but does not close; Huber still clearly worse)
If masked Huber >= masked ULM, the README line "signal lives in the strongly responding targets" needs
revising, robust losses reopen, and a degree-preserving shuffled-network control is then owed (not run
here: no gain is being claimed, and controls.py's shuffle is under repair in T-001).

Imports robust_ulm, which reruns the unmasked comparison (~4 min threaded, was ~25 serial) -- repo convention, and it gives the
unmasked side of the table for free. The known defects in robust_ulm's Huber (unweighted intercept,
uncentered MAD) are deliberately left alone: V-003 showed they do not drive the result, and the point is
to compare against the same implementation F-003 used.

T-025 adds the mechanism V-008 measured from scratch: the 197/82 self-edge split of all four arms, and
the final Huber weight on the perturbed TF's own gene in its own regression (V-008: median 0.094,
IQR 0.051-0.153 over the unmasked self-edge experiments). robust_ulm.huber_all now takes `cells` and
returns that one weight per experiment instead of discarding the whole weight matrix.
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/masked_huber.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd
from robust_ulm import (adata, net, huber_all, W, Y, tfs, S as S_ulm, S2 as S_hub, evaluate,
                        f1 as fr_ulm, f2 as fr_hub, w_own, genes)

# --- mask the perturbed TF's own gene, exactly as leak.py does -------------------------------------
X = adata.to_df()
hits = [(X.index.get_loc(exp), X.columns.get_loc(tf))
        for exp, tf in adata.obs.source.items() if tf in X.columns]
for exp, tf in adata.obs.source.items():
    if tf in X.columns: X.loc[exp, tf] = 0
rows = np.array([r for r, _ in hits]); cols = np.array([c for _, c in hits])

masked = ad.AnnData(X, obs=adata.obs)
dc.mt.ulm(masked, net, tmin=5)
S_ulm_m = masked.obsm["score_ulm"].reindex(columns=tfs)

Ym = np.nan_to_num(np.asarray(X.values, dtype=float))          # same nan handling as robust_ulm
T = huber_all(Ym)
S_hub_m = pd.DataFrame(T, index=S_ulm.index, columns=tfs)

fr_ulm_m = evaluate(S_ulm_m, "ulm masked")
fr_hub_m = evaluate(S_hub_m, "huber masked")

def row(name, fr, base):
    print(f"{name:22s} {np.mean(fr<=.1):.2f}   {np.mean(fr<=.05):.2f}  "
          f"{'-' if base is None else format(np.mean(fr<base), '.2f'):>6s} "
          f"{'-' if base is None else format(np.mean(fr>base), '.2f'):>6s}")
print(f"\n{'variant':22s} top10  top5  better  worse   (better/worse = vs ULM in the same mask condition)")
row("ulm unmasked", fr_ulm, None)
row("huber unmasked", fr_hub, fr_ulm)
row("ulm MASKED", fr_ulm_m, None)
row("huber MASKED", fr_hub_m, fr_ulm_m)
print(f"\nn={len(fr_ulm)} experiments scored; {len(hits)} of {adata.n_obs} experiments had their own gene zeroed")
print(f"unmasked deficit (ulm - huber, top10) = {np.mean(fr_ulm<=.1) - np.mean(fr_hub<=.1):+.2f}; "
      f"masked deficit = {np.mean(fr_ulm_m<=.1) - np.mean(fr_hub_m<=.1):+.2f}")
if np.mean(fr_hub_m <= .1) >= np.mean(fr_ulm_m <= .1):
    print("NOTE: masked Huber >= masked ULM. A degree-preserving shuffled-network control is now owed.")

# --- T-025: where the deficit lives, and the weight that puts it there ------------------------------
gi = {g: i for i, g in enumerate(genes)}
exps = [(e, t) for e, t in adata.obs.source.items() if t in tfs]        # the 279 evaluate() scores, same order
has_self = np.array([t in gi and W[gi[t], tfs.get_loc(t)] != 0 for _, t in exps])
scored = np.array([t in tfs for t in adata.obs.source])

print(f"\nself-edge split of the four arms (V-010's condition: W[gene X, column X] != 0), top-10% hits")
print(f"{'group':22s} {'n':>4s}  {'ulm':>5s} {'huber':>6s} {'ulm_m':>6s} {'hub_m':>6s}")
for lab, m in [("self-edge", has_self), ("no self-edge", ~has_self)]:
    print(f"{lab:22s} {m.sum():4d}  " + " ".join(f"{np.sum(f[m] <= .1):5d} " for f in
          (fr_ulm, fr_hub, fr_ulm_m, fr_hub_m)))
ws = w_own[scored][has_self]
print(f"final Huber weight on the perturbed TF's own gene in its own regression, unmasked, self-edge "
      f"experiments (n={len(ws)}): median {np.median(ws):.3f}, IQR {np.percentile(ws, 25):.3f}-"
      f"{np.percentile(ws, 75):.3f}; below 1 in {np.mean(ws < 1):.0%}, below 0.5 in {np.mean(ws < 0.5):.0%}")
print(f"  same weight in the {(~has_self).sum()} no-self-edge experiments (own gene is a non-target there): "
      f"median {np.median(w_own[scored][~has_self]):.3f}")

# the mask is the whole experiment: it must hit exactly one cell per scoreable row, the TF's own column,
# the Huber arm must see the same masked matrix the ULM arm saw, and (T-025) the masked Huber scores must
# actually differ from robust_ulm's unmasked ones wherever the masked cell can reach the TF's own score.
# F-029: persist. One row per scoreable experiment, all four arms' rank fractions and the
# self-edge split they are read on. Written before the assert; no printed line changes.
# NOT YET GENERATED: this script takes about 8 minutes and was not rerun when F-029 was filed.
pd.DataFrame(dict(exp=[e for e, _ in exps], tf=[t for _, t in exps], has_self_edge=has_self,
                  ulm=fr_ulm, huber=fr_hub, ulm_masked=fr_ulm_m, huber_masked=fr_hub_m,
                  huber_weight_own_gene=w_own[scored])
             ).to_csv("../results/masked_huber.csv", index=False)

own_m = np.array([S_hub_m.loc[e, t] for e, t in exps]); own_u = np.array([S_hub.loc[e, t] for e, t in exps])
assert (Ym[rows, cols] == 0).all() and len(hits) == len(set(rows)) and \
       np.array_equal(np.argwhere(Ym != np.nan_to_num(Y)), np.c_[rows, cols][np.nan_to_num(Y)[rows, cols] != 0]) \
       and (own_m[has_self] != own_u[has_self]).all()
