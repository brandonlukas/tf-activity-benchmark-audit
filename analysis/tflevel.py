"""T-018 / T-009 / T-020 / F-014: everything so far at ONE ROW PER TF.

279 knockTF experiments cover 155 TFs, so every rate in F-001..F-011 treats repeated TFs as
independent (V-002, V-010, V-011 all flag it). Here each TF gets the MEDIAN rank fraction over
its experiments and hit = median fraction <= 0.10 (top-5%: <= 0.05). Four blocks:

  (1) TF-level baselines: ULM full network / self-edge removed / own gene masked.
  (2) F-002's stratification at TF level: miss rate by regulon-size tercile, by max-Jaccard
      overlap tercile, and size x overlap (V-002's "would flip if" was exactly this).
  (3) size-matched self-edge comparison AFTER self-edge removal: each no-self-edge TF matched
      to its nearest self-edge TF on log10 target count (with replacement), difference in hit
      rate with a paired bootstrap CI over the matched pairs (numpy, seed 0).
  (4) the 12 TFs whose self-edge is -1 (V-011): target counts, and whether a size-matched
      comparison against the +1 TFs is answerable at all.

Fractions, conditions and the ULM kernel are selfedge.py's, unchanged (imported, which reruns
propagate.py/controls.py/leak.py); regulon size and max Jaccard are gapcheck.py's, unchanged.

PREDICTIONS (stated before running):
  (1) TF-level ULM 0.40 -> 0.28 masked (experiment level is 0.45 -> 0.34; V-010's scratch
      TF-level run gave 0.398 -> 0.277, self-edge TFs 0.524 -> 0.328, others 0.194).
  (2) overlap terciles flat or FALLING at TF level (no rise: V-002's flip condition is not met);
      the size effect persists (small regulons miss more).
  (3) the size-matched self-edge-removed difference is within +/-0.07 of zero.
  (4) 3 or fewer -1 TFs per size bin, i.e. unanswerable in knockTF.
No gain over ULM is claimed, so no shuffled network; masked is in the table (F-010's "no
self-edge" is the same thing for plain ULM, and both are shown).

--- T-029 / F-019 additions (V-014 asked for all four; blocks (5)-(8), appended so every block
above prints byte-identically to the run V-014 checked) -------------------------------------
MEDIAN CONVENTION: pandas' median averages the two middle values for even n. 26 TFs have
exactly 2 experiments, and for those "median <= 0.10" is in effect "both experiments hit" when
the two straddle the line (V-014 (a): 11 straddling TFs under full-network ULM, all counted as
misses). Block (8) prints the lower-median variant so the convention's cost is visible.

  (5) TIE-AVERAGED matcher. Target counts are integers and 25 of the 59 no-self-edge TFs have
      two or three equally near self-edge TFs; np.argmin in block (3) breaks those ties in
      alphabetical TF order, and that choice decides the SIGN of the -0.017 (V-014 (c): drop
      BCLAF1 and it becomes +0.136). Here each no-self-edge TF is matched to ALL self-edge TFs
      at the minimum |log10 size| distance and their hit indicators are averaged.
  (6) TWO-POOL BOOTSTRAP. Block (3) resamples the 59 matched pairs only, so it carries neither
      the self-edge pool's sampling nor the match. Here the 59 no-self-edge and 96 self-edge TFs
      are resampled independently and RE-MATCHED inside each draw (numpy, seed 0, 5000 draws),
      for the matched difference, the unmatched difference, and their difference.
  (7) F-002 at the EXPERIMENT-level cut points. TF-level terciles cut at 27 / 90.7 targets where
      experiments cut at 48.3 / 110, so F-014's "gets stronger for size" compared different
      strata (V-014 (b)). Same cut points for both levels here.
  (8) lower-median collapse, the other side of the n=2 convention.

PREDICTIONS for (5)-(8), stated before running, all from V-014:
  (5) tie-averaged matched 0.203 vs 0.186, difference +0.017 (block (3)'s argmin gives -0.017).
  (6) tie-averaged matched difference CI about [-0.169, +0.233]; unmatched difference +0.125,
      CI [-0.011, +0.260]; unmatched minus matched +0.109, CI [-0.079, +0.271], so the CI
      contains 0 and the matched CI contains the unmatched gap.
  (7) TF level at the experiment cut points: miss 0.736 (n=72) / 0.565 (n=46) / 0.541 (n=37),
      low minus high 0.195 against 0.204 for experiments at the same cuts -- the SAME gradient,
      not a stronger one; overlap 0.694 / 0.700 / 0.488, still no rise.
  (8) lower median: ulm 70/155 (0.452), no self-edge 46/155, masked 46/155.
The claim this backs (F-019) is the WEAKENED one: the matched difference is about 0, its CI
contains both 0 and the unmatched 0.126, and 59 no-self-edge TFs with 11 hits cannot settle
whether regulon size explains the residual self-edge gap. NOT "matching removes the gap".
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/tflevel.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import numpy as np, pandas as pd
import selfedge as se                  # reruns propagate.py, controls.py, leak.py, selfedge.py
import gapcheck as gc                  # regulon size + max pairwise Jaccard per TF

# ---- per-experiment frame, then collapse to one row per TF -------------------------------
exps = se.exps                                             # [(experiment, perturbed TF)], 279
df = pd.DataFrame({"exp": [e for e, _ in exps], "tf": [t for _, t in exps],
                   "ulm": se.fr["ulm (full net)"], "noself": se.fr["no self-edge"],
                   "masked": se.fr["masked (leak.py)"], "has_self": se.has_self})
df["sign"] = [np.sign(se.adjm_orig[se.var_i[t], se.src_i[t]]) if t in se.var_i else 0.0
              for _, t in exps]
g = gc.df.set_index("exp")
df = df.join(g[["size", "max_jac", "frac"]].rename(columns={"frac": "gc_frac"}), on="exp")

VAR = ["ulm", "noself", "masked"]
tl = df.groupby("tf").agg(n=("exp", "size"), size=("size", "first"), max_jac=("max_jac", "first"),
                          has_self=("has_self", "first"), sign=("sign", "first"),
                          **{v: (v, "median") for v in VAR})
for v in VAR:
    tl["hit10_" + v] = tl[v] <= .10
    tl["hit5_" + v] = tl[v] <= .05

def row(tab, label, v, ref="ulm"):
    """top-10%, hits/n, top-5%, paired better/worse vs ULM on the same TFs."""
    hit10, hit5 = tab["hit10_" + v], tab["hit5_" + v]
    better, worse = np.mean(tab[v] < tab[ref]), np.mean(tab[v] > tab[ref])
    print(f"{label:34s} {hit10.mean():.3f} {hit10.sum():4d}/{len(tab):<4d} {hit5.mean():.3f}"
          f"  {better:.2f}   {worse:.2f}")

print(f"\n=== TF LEVEL: {len(tl)} TFs over {len(df)} experiments "
      f"(median experiments/TF {tl.n.median():.0f}, max {tl.n.max()}; "
      f"{(tl.n > 1).sum()} TFs appear more than once) ===")

print("\n(1) baselines, one row per TF (median rank fraction; better/worse vs ULM, per TF)")
print(f"{'variant':34s} top10   hits/n    top5  better worse")
for v in VAR: row(tl, {"ulm": "ulm (full net)", "noself": "no self-edge", "masked": "masked"}[v], v)
print("    self-edge TFs only")
for v in VAR: row(tl[tl.has_self], "  " + v, v)
print("    no-self-edge TFs only")
for v in VAR: row(tl[~tl.has_self], "  " + v, v)
print("  experiment level for reference (F-001/F-010): "
      f"ulm {np.mean(df.ulm <= .1):.3f}  no self-edge {np.mean(df.noself <= .1):.3f}  "
      f"masked {np.mean(df.masked <= .1):.3f}")
alt = {v: df.assign(h=df[v] <= .1).groupby("tf").h.mean().mean() for v in VAR}
print("  other collapse rule (mean over TFs of that TF's experiment-level hit rate, which is "
      f"what V-010's scratch run reported): ulm {alt['ulm']:.3f}  no self-edge {alt['noself']:.3f}"
      f"  masked {alt['masked']:.3f}")

# ---- (2) F-002 stratification at TF level ------------------------------------------------
print("\n(2) F-002 stratification at TF level (miss = median fraction > 0.10; terciles over TFs)")
for col in ["size", "max_jac"]:
    tl[col + "_bin"] = pd.qcut(tl[col], 3, labels=["low", "mid", "high"], duplicates="drop")
    t = tl.groupby(col + "_bin", observed=True).agg(
        n=("n", "size"), exps=("n", "sum"),
        miss_ulm=("hit10_ulm", lambda s: 1 - s.mean()),
        miss_noself=("hit10_noself", lambda s: 1 - s.mean()),
        miss_masked=("hit10_masked", lambda s: 1 - s.mean()),
        med=(col, "median"))
    print(f"\nby {col} (TF level):\n{t.round(3).to_string()}")
    e = gc.df.assign(b=pd.qcut(gc.df[col], 3, labels=["low", "mid", "high"], duplicates="drop"))
    print("  experiment level (F-002/V-002): ",
          e.groupby("b", observed=True).miss.mean().round(2).to_dict())
piv = tl.pivot_table(index="size_bin", columns="max_jac_bin", values="hit10_ulm",
                     aggfunc=lambda s: 1 - s.mean(), observed=True)
cnt = tl.pivot_table(index="size_bin", columns="max_jac_bin", values="n", aggfunc="size", observed=True)
print(f"\nsize x overlap, ULM miss rate (TF level):\n{piv.round(2).to_string()}"
      f"\nn TFs per cell:\n{cnt.to_string()}")
print("  Spearman(median fraction, max_jac) = "
      f"{pd.Series(tl.ulm).corr(tl.max_jac, method='spearman'):+.3f}; "
      f"vs size = {pd.Series(tl.ulm).corr(tl['size'], method='spearman'):+.3f}")

# ---- (3) size-matched self-edge vs no-self-edge, after self-edge removal ------------------
# Judgment call: nearest neighbour on log10(target count), WITH replacement, each of the 59
# no-self-edge TFs matched to one self-edge TF. Size bins would leave the top bin with 4 TFs
# (V-010); NN keeps every no-self-edge TF and makes the match gap reportable.
se_tf, ns_tf = tl[tl.has_self], tl[~tl.has_self]
ls, ln = np.log10(se_tf["size"].values), np.log10(ns_tf["size"].values)
idx = np.abs(ln[:, None] - ls[None, :]).argmin(1)
gap = np.abs(ln - ls[idx])
m_hit = se_tf["hit10_noself"].values[idx]              # matched self-edge TFs, self-edge REMOVED
n_hit = ns_tf["hit10_noself"].values                   # identical to their ULM row by construction
d = m_hit.mean() - n_hit.mean()
rng = np.random.default_rng(0)
b = rng.integers(len(n_hit), size=(10000, len(n_hit)))
boot = m_hit[b].mean(1) - n_hit[b].mean(1)
lo, hi = np.percentile(boot, [2.5, 97.5])
print(f"\n(3) size-matched, self-edge REMOVED (nearest neighbour on log10 target count, "
      f"with replacement)\n    no-self-edge TFs n={len(ns_tf)} median {ns_tf['size'].median():.0f} "
      f"targets, hit10 {n_hit.mean():.3f}"
      f"\n    matched self-edge TFs (distinct {len(set(idx))} of {len(se_tf)}) median "
      f"{se_tf['size'].values[idx].mean():.0f} mean / {np.median(se_tf['size'].values[idx]):.0f} "
      f"median targets, hit10 {m_hit.mean():.3f}"
      f"\n    difference {d:+.3f}  bootstrap 95% CI over the {len(n_hit)} matched pairs "
      f"[{lo:+.3f}, {hi:+.3f}] (seed 0, 10000 draws)"
      f"\n    match quality: median |log10 size gap| {np.median(gap):.3f} "
      f"(x{10**np.median(gap):.2f}), max {gap.max():.3f} (x{10**gap.max():.2f}); "
      f"{np.mean(gap > np.log10(2)):.0%} of pairs differ by more than 2x"
      f"\n    unmatched for contrast: self-edge TFs with self-edge removed "
      f"{se_tf['hit10_noself'].mean():.3f} (n={len(se_tf)}), unmasked {se_tf['hit10_ulm'].mean():.3f}")

# ---- (4) the -1 self-edge TFs ------------------------------------------------------------
neg = tl[tl.sign < 0].sort_values("size", ascending=False)
pos = tl[tl.sign > 0]
print(f"\n(4) the {len(neg)} TFs with a -1 CollecTRI self-edge (V-011), vs {len(pos)} with +1")
print(neg[["n", "size", "max_jac", "ulm", "noself", "masked"]].round(3).to_string())
print(f"    -1 median {neg['size'].median():.0f} targets, +1 median {pos['size'].median():.0f}, "
      f"no self-edge median {ns_tf['size'].median():.0f}")
bins = [0, 60, 150, 400, 10 ** 9]
lab = ["<=60", "61-150", "151-400", ">400"]
tl["sbin"] = pd.cut(tl["size"], bins, labels=lab)
tab = pd.DataFrame({"neg_n": tl[tl.sign < 0].groupby("sbin", observed=False).size(),
                    "pos_n": tl[tl.sign > 0].groupby("sbin", observed=False).size(),
                    "noself_n": tl[~tl.has_self].groupby("sbin", observed=False).size(),
                    "neg_hit10_noself": tl[tl.sign < 0].groupby("sbin", observed=False).hit10_noself.mean(),
                    "pos_hit10_noself": tl[tl.sign > 0].groupby("sbin", observed=False).hit10_noself.mean()})
print(f"\nsize bins (targets):\n{tab.round(2).to_string()}")
worst = tab.neg_n.max()
print(f"    largest -1 cell = {worst} TFs; V-011 expects 3 or fewer per bin -> "
      + ("ANSWERABLE" if worst >= 10 else "UNANSWERABLE IN KNOCKTF: no size bin holds enough -1 "
         "TFs to match against the +1 TFs; counts above, no number forced."))

# ---- (5) tie-averaged matcher (T-029; argmin's alphabetical tie-break decided (3)'s sign) ----
se_hit = se_tf["hit10_noself"].values.astype(float)


def tie_avg(ln_, ls_, hit_):
    """Match each no-self-edge TF to ALL nearest self-edge TFs on log10 size; average their hits.
    Returns (per-TF averaged hit indicator, tie weights, |log10 size gap| of the matches)."""
    D = np.abs(ln_[:, None] - ls_[None, :])
    tie = D == D.min(1, keepdims=True)
    W = tie / tie.sum(1, keepdims=True)
    return W @ hit_, W, D.min(1)


m_avg, W, gap_avg = tie_avg(ln, ls, se_hit)
d_avg = m_avg.mean() - n_hit.mean()
nties = W.astype(bool).sum(1)
print(f"\n(5) TIE-AVERAGED size match, self-edge REMOVED (each no-self-edge TF matched to ALL "
      f"self-edge TFs at the minimum |log10 target count| distance)"
      f"\n    ties per no-self-edge TF: median {np.median(nties):.0f}, max {nties.max()}; "
      f"{np.mean(nties > 1):.0%} of the {len(ln)} TFs have more than one nearest self-edge TF"
      f"\n    distinct self-edge TFs used {len(np.flatnonzero(W.sum(0)))} of {len(ls)} "
      f"(argmin in (3) used {len(set(idx))}); median |log10 gap| {np.median(gap_avg):.3f} "
      f"(x{10 ** np.median(gap_avg):.2f}), max {gap_avg.max():.3f} (x{10 ** gap_avg.max():.2f})"
      f"\n    matched self-edge {m_avg.mean():.3f} vs no-self-edge {n_hit.mean():.3f}, "
      f"difference {d_avg:+.3f}   ((3)'s argmin: {d:+.3f}; unmatched {se_hit.mean():.3f} vs "
      f"{n_hit.mean():.3f} = {se_hit.mean() - n_hit.mean():+.3f})")

# ---- (6) two-pool re-matching bootstrap ---------------------------------------------------
rng2 = np.random.default_rng(0)
NDRAW = 5000
bm, bu = np.empty(NDRAW), np.empty(NDRAW)
for k in range(NDRAW):
    i = rng2.integers(len(ls), size=len(ls))          # self-edge pool, resampled
    j = rng2.integers(len(ln), size=len(ln))          # no-self-edge pool, resampled
    mh, _, _ = tie_avg(ln[j], ls[i], se_hit[i])       # re-matched inside the draw
    bm[k] = mh.mean() - n_hit[j].mean()
    bu[k] = se_hit[i].mean() - n_hit[j].mean()
ci = lambda v: np.percentile(v, [2.5, 97.5])
lo_m, hi_m = ci(bm); lo_u, hi_u = ci(bu); lo_g, hi_g = ci(bu - bm)
print(f"\n(6) two-pool bootstrap (both TF pools resampled independently and RE-MATCHED inside "
      f"each draw; tie-averaged matcher, numpy seed 0, {NDRAW} draws)"
      f"\n    matched difference     {d_avg:+.3f}  95% CI [{lo_m:+.3f}, {hi_m:+.3f}]  "
      f"(block (3)'s pairs-only CI was [{lo:+.3f}, {hi:+.3f}], width {hi - lo:.3f} vs "
      f"{hi_m - lo_m:.3f})"
      f"\n    unmatched difference   {se_hit.mean() - n_hit.mean():+.3f}  95% CI "
      f"[{lo_u:+.3f}, {hi_u:+.3f}]"
      f"\n    unmatched minus matched {bu.mean() - bm.mean():+.3f}  95% CI "
      f"[{lo_g:+.3f}, {hi_g:+.3f}]   P(<= 0) = {np.mean(bu - bm <= 0):.2f}"
      f"\n    P(matched draw >= the unmatched {se_hit.mean() - n_hit.mean():.3f}) = "
      f"{np.mean(bm >= se_hit.mean() - n_hit.mean()):.2f}; matched CI contains 0: "
      f"{bool(lo_m <= 0 <= hi_m)}, contains the unmatched gap: "
      f"{bool(lo_m <= se_hit.mean() - n_hit.mean() <= hi_m)}"
      f"\n    -> the matched difference is about 0, but these data cannot tell 'size matching "
      f"removes the gap' from 'size matching does nothing to it' (V-014 (d))")

# ---- (7) F-002 at the EXPERIMENT-level cut points (V-014 (b): TF terciles cut elsewhere) ----
print("\n(7) same stratification at the EXPERIMENT-level tercile cut points, TF rows vs "
      "experiment rows (F-014's (2) used TF-level cut points, which are different strata)")
for col in ["size", "max_jac"]:
    q = np.quantile(gc.df[col].values, [1 / 3, 2 / 3])
    edges = [-np.inf, q[0], q[1], np.inf]
    tl[col + "_ebin"] = pd.cut(tl[col], edges, labels=["low", "mid", "high"])
    e = gc.df.assign(b=pd.cut(gc.df[col], edges, labels=["low", "mid", "high"]))
    t = tl.groupby(col + "_ebin", observed=True).agg(
        n_tf=("n", "size"), miss_ulm=("hit10_ulm", lambda s: 1 - s.mean()),
        miss_noself=("hit10_noself", lambda s: 1 - s.mean()), med=(col, "median"))
    t["n_exp"] = e.groupby("b", observed=True).size()
    t["miss_exp"] = e.groupby("b", observed=True).miss.mean()
    print(f"\nby {col}, cut at {q[0]:.3g} / {q[1]:.3g} (experiment-level terciles):\n"
          f"{t.round(3).to_string()}")
    lmh_tf = t.miss_ulm.iloc[0] - t.miss_ulm.iloc[-1]
    lmh_ex = t.miss_exp.iloc[0] - t.miss_exp.iloc[-1]
    own = 1 - tl.groupby(col + "_bin", observed=True).hit10_ulm.mean()   # block (2)'s TF terciles
    print(f"  low minus high: TF level {lmh_tf:+.3f}, experiment level {lmh_ex:+.3f} -> "
          + ("same gradient" if abs(lmh_tf - lmh_ex) < 0.05 else
             f"gradients differ by {lmh_tf - lmh_ex:+.3f}")
          + f"; at block (2)'s TF-level cut points it was {own.iloc[0] - own.iloc[-1]:+.3f}")

# ---- (8) lower-median collapse (the other side of the n=2 convention) ----------------------
lower = df.groupby("tf")[VAR].agg(lambda s: np.sort(s.values)[(len(s) - 1) // 2])
for v in VAR:
    lower["hit10_" + v] = lower[v] <= .10
    lower["hit5_" + v] = lower[v] <= .05
print("\n(8) lower-median collapse (ties among a TF's experiments go the OTHER way; pandas' "
      "median averages the two middle values, so 26 two-experiment TFs need both to hit)")
print(f"{'variant':34s} top10   hits/n    top5  better worse")
for v in VAR:
    row(lower, {"ulm": "lower-median ulm", "noself": "lower-median no self-edge",
                "masked": "lower-median masked"}[v], v)
print(f"    for comparison, the median rule above: ulm {tl.hit10_ulm.sum()}/155, no self-edge "
      f"{tl.hit10_noself.sum()}/155, masked {tl.hit10_masked.sum()}/155")

# F-029: persist. The per-TF collapse (155 rows) and the per-experiment table behind it.
tl.to_csv("../results/tflevel.csv")
df.to_csv("../results/tflevel_rows.csv", index=False)

assert len(tl) == 155 and np.abs(df.ulm.values - df.gc_frac.values).max() < 1e-12 \
    and (tl.loc[~tl.has_self, "ulm"] == tl.loc[~tl.has_self, "noself"]).all() \
    and len(idx) == len(ns_tf) \
    and se_tf["has_self"].values[idx].all() and gap.max() < np.log10(2) \
    and W.sum(1).round(12).tolist() == [1.0] * len(ln) \
    and (se_tf["has_self"].values[W.astype(bool).any(0)]).all() \
    and gap_avg.max() < np.log10(2), \
    "TF collapse lost rows, selfedge/gapcheck fractions disagree, self-edge removal touched a " \
    "no-self-edge TF, or the matching is broken (a matched TF without a self-edge, a size gap " \
    "beyond 2x, or tie weights that do not sum to 1)"
