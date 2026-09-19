"""T-017 / F-011: does the SIGN of the perturbed TF's CollecTRI self-edge decide who the
unmasked knockTF benchmark rewards?

V-010 left this open. Every knockTF knockdown has own logFC < 0 (V-005) and the benchmark ranks
type_p * score with type_p = -1, so through the self-edge X->X:
  weight +1  ->  the own-mRNA drop lowers TF X's ULM score  ->  -1 * score goes UP  -> rewarded
  weight -1  ->  the own-mRNA drop raises TF X's ULM score  ->  -1 * score goes DOWN -> penalised
i.e. positive autoregulators should be inflated by the leak and negative autoregulators deflated
by it. Splitting F-010's 197 self-edge experiments by sign tests that directly: the "no self-edge"
condition (network minus X->X only, gene left in the matrix) is exactly the leak's removal.

Conditions and evaluation are selfedge.py's, imported, not reimplemented; this script only regroups
its per-experiment rank fractions and adds regulon size. No gain is claimed, so no shuffled network:
the masked condition is the control, alongside "no self-edge".

PREDICTIONS (stated before running):
  (a) positive-self-edge group: removing the self-edge LOSES hits, roughly 0.56 -> 0.38 top-10%;
  (b) negative-self-edge group: removing it GAINS or is flat, predict about +3 points, i.e. the
      unmasked benchmark penalises negative autoregulators;
  (c) masked tracks "no self-edge" in both groups (F-010: for plain ULM they are the same edge);
  (d) V-010 says the self-edge effect shrinks with regulon size, so both effects should be larger
      in the regulon <= 60 targets subset.
If the negative group is under ~20 experiments, counts are reported and percentages are not to be
read as rates.

Runtime ~30 s (V-011 measured it; the earlier "about 2.5 min" note in F-011 is wrong).
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/selfsign.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import numpy as np, pandas as pd
from selfedge import (adata, adjm_orig, es_full, es_ns, exps, fr, has_self, src_i, tfs,
                      var_i)   # reruns F-010's tables

w_self = np.array([adjm_orig[var_i[tf], src_i[tf]] if tf in var_i else 0.0 for _, tf in exps])
pos, neg = w_self > 0, w_self < 0
reg_size = np.array([(adjm_orig[:, src_i[tf]] != 0).sum() for _, tf in exps])
type_p = np.array([adata.obs.loc[e, "type_p"] for e, _ in exps])

VAR = ["ulm (full net)", "no self-edge", "masked (leak.py)"]
ref = fr["ulm (full net)"]

def block(title, m):
    n_tf = len({t for (_, t), k in zip(exps, m) if k})
    print(f"\n{title}: n={m.sum()} experiments, {n_tf} TFs, median regulon "
          f"{np.median(reg_size[m]):.0f} targets, type_p values {sorted(set(type_p[m]))}")
    if m.sum() == 0: return
    print(f"  {'variant':22s} top10  hits/n     top5  hits/n    better worse")
    for k in VAR:
        f, r = fr[k][m], ref[m]
        print(f"  {k:22s} {np.mean(f<=.1):.3f} {np.sum(f<=.1):3d}/{m.sum():<4d} "
              f"{np.mean(f<=.05):.3f} {np.sum(f<=.05):3d}/{m.sum():<4d} "
              f"{np.mean(f<r):.2f}   {np.mean(f>r):.2f}")

print(f"\n{'='*78}\n(4) T-017: the {has_self.sum()} self-edge experiments split by self-edge SIGN "
      f"(better/worse = that variant vs unmasked full-network ULM, within the group)")
print(f"distinct self-edge weights: {sorted(set(w_self[has_self]))}; "
      f"{pos.sum()} positive, {neg.sum()} negative, "
      f"{(~has_self).sum()} experiments with no self-edge (unchanged by the conditions)")
block("POSITIVE self-edge", pos)
block("NEGATIVE self-edge", neg)
block("no self-edge (reference)", ~has_self)

print(f"\n    within-sign, regulon <= 60 targets (V-010: the self-edge effect shrinks with size)")
block("POSITIVE self-edge, regulon <= 60", pos & (reg_size <= 60))
block("NEGATIVE self-edge, regulon <= 60", neg & (reg_size <= 60))
if (neg & (reg_size <= 60)).sum() < 20:
    print("    NOTE: negative-sign cells below ~20 experiments; read the counts, not the rates")

# T-022: the sign split must partition F-010's self-edge group AND the self-edge must move the
# perturbed TF's own score the way its sign predicts -- down for +1, up for -1, in every experiment.
# F-029: persist. One row per experiment with the sign split and the three variants' ranks.
pd.DataFrame(dict(exp=[e for e, _ in exps], tf=[t for _, t in exps], self_edge_weight=w_self,
                  has_self_edge=has_self, sign=np.where(pos, "+1", np.where(neg, "-1", "none")),
                  regulon_size=reg_size, type_p=type_p,
                  **{k.split(" (")[0].replace(" ", "_"): fr[k] for k in VAR})
             ).to_csv("../results/selfsign.csv", index=False)

own_full = np.array([es_full.loc[e, tf] for e, tf in exps])
own_ns = np.array([es_ns.loc[e, tf] for e, tf in exps])
assert ((pos ^ neg) == has_self).all() and (own_full[pos] < own_ns[pos]).all() \
    and (own_full[neg] > own_ns[neg]).all(), \
    "sign split does not partition F-010's self-edge group, or a self-edge moves the TF's own score the wrong way"
