"""T-033 / F-021: do CollecTRI's SIGNS carry information for masked ULM, and do the default-filled ones?

V-017's skeptic ran a scratch probe on MASKED ULM (the perturbed TF's own gene zeroed in its own
row; baseline 96/279 top-10%, 78/279 top-5%) and found that erasing every sign costs 15 hits while
deleting or "correcting" the default-activation fills costs nothing. That probe was not a committed
experiment, had no stated prediction, and -- the reason for this script -- did NOT re-prune the
network to tmin=5 after removing edges, so its "default edges removed" arm scored some TFs on
regulons of 1-4 measured targets that decoupler would never have scored.

Part A. Five network interventions, every one scored on leak.py's own-gene-masked matrix with
decoupler's own ULM kernel and RE-PRUNED to tmin=5 whenever edges are removed:
  1 sign-erased      every weight set to +1
  2 no-default       default-activation edges dropped
  3 flip-default-rep default edges set to -1 for TFs whose PMID-signed edges are majority
                     repressing (V-016's definition: PMID -1 count > PMID +1 count), others left
  4 drop-default-rep default edges dropped for those TFs only
  5a rand-flip       control for 1: a RANDOM set of edges of the size of the -1 set (5822) has its
                     sign flipped, 5 seeds -- does erasing the real signs cost more than scrambling
                     an equal number of arbitrary ones?
  5b rand-drop-perTF control for 2: a random set of edges of the size of the default set, matched
                     PER TF (the same number removed from each TF's regulon), 5 seeds -- so
                     "removing the defaults" is compared with "removing that many edges".
Re-pruning changes the scored source set, so each variant reports (i) hits of all 279 with
unscoreable experiments counted as misses and (ii) the paired common scoreable set, with
better/worse vs masked ULM, exact McNemar on hit transitions and a TF-cluster bootstrap CI
(numpy, seed 0, 4000 resamples) on the top-10% difference.

PREDICTIONS for Part A (stated before running; 1-4 are V-017's un-re-pruned probe numbers, so the
re-pruned values may differ and that difference is the point of the exercise):
  1  81/279, 31 hits lost / 16 gained
  2  89/279
  3  96/279, 7 lost / 7 gained
  4  99/279
  5a random sign flips of the same count: 85/279. (That random flips could cost MORE than erasing
     the real signs is not obvious, so this is a committed number, not a hedge: erasing signs sets
     every edge to the majority sign, which is right for 86% of edges, whereas flipping 5822 random
     edges makes ~4996 of them actively wrong. I predict random flips hurt LESS than sign erasure
     is FALSE -- 85 < 81 is the prediction, i.e. random flips hurt slightly LESS. )
  5b per-TF-matched random removal: 90/279, i.e. indistinguishable from removing the defaults.

Part B. The post hoc lead in V-017 (7): TFs whose PMID-signed edges are majority repressing AND
that have default edges score 5/43 masked top-10% vs 66/154 for PMID-majority-activating TFs with
default edges (25 vs 85 TFs; TF-label permutation p = 0.0013, 0.0054 within size bin; but rank-based
TF-level Mann-Whitney p = 0.14). It was found by looking, so it gets a pre-stated test here.
  PRIMARY  : TF-level comparison of the MEDIAN masked rank fraction of the two classes, permutation
             of the class label among TFs WITHIN regulon-size bin (<=60, 61-150, >150), two-sided.
  SECONDARY: the experiment-level top-10% rate, same permutation scheme.
  PREDICTION (primary): repressor-majority TFs rank worse, TF-level median rank fraction 0.39 vs
  0.32, p about 0.14, i.e. NOT significant on the primary test. The secondary is expected to stay
  significant, which would mean the lead lives in the top-10% cut and not in the ranking.
  Confounds checked: knockdown efficiency (the perturbed TF's own logFC in the two classes) and
  regulon size. Mechanistic check: what fraction of a repressor-majority TF's scored regulon is
  default +1 vs PMID -1, and does scoring those TFs on PMID-signed edges only (re-pruned, common
  scoreable set) change their masked hit rate?

Controls: the own gene is masked throughout -- that IS the baseline here, no unmasked number is
quoted as a result. The degree-preserving shuffled network from controls.py is added for
interventions 1 and 2 only if one of them comes out ABOVE masked ULM's 96/279 by more than its
bootstrap CI; otherwise no gain is claimed and the row would be meaningless.

Loading, the ULM kernel and evaluate() are selfedge.py's / leak.py's, imported or copied verbatim,
not reimplemented. Runtime ~4 min (the import chain reruns propagate/controls/leak/selfedge).
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/signprobe.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import numpy as np, pandas as pd, decoupler as dc
from scipy import stats
from selfedge import (adata, mat, obs, var, obs_i, var_i, src_i, adjm_orig, exps, tfs, pnet)
import leak                                    # already imported by selfedge; masked-ULM baseline

net = pd.read_parquet("../data/collectri.parquet")
N_NEG = int((net.weight < 0).sum())
DEF = net.sign_decision == "default activation"

# ---------------------------------------------------------------- masked matrix (leak.py's)
mat_m = mat.copy()
for e, tf in adata.obs.source.items():
    if tf in var_i: mat_m[obs_i[e], var_i[tf]] = 0        # erase the knockdown's own mRNA drop

def build_score(net2):
    """prune to tmin=5 on the INTERVENED network, then decoupler's own adjacency + ULM kernel."""
    pn = dc.pp.prune(features=var, net=net2, tmin=5)
    s, _f, a = dc.pp.adjmat(features=var, net=pn)
    return pd.DataFrame(dc.mt.ulm.func(mat_m, a)[0], index=obs, columns=s)

def evaluate(S):
    """copied from selfedge.py / leak.py; returns exp -> rank fraction of the true TF."""
    out = {}
    for exp, tf in adata.obs.source.items():
        if tf not in S.columns: continue
        s = adata.obs.loc[exp, "type_p"] * S.loc[exp]
        out[exp] = s.rank(ascending=False)[tf] / S.shape[1]
    return out

base = evaluate(build_score(net))
EXPS = [e for e, _ in exps]                                # the 279 scoreable under the real network
TF_OF = {e: t for e, t in exps}
fr_b = np.array([base[e] for e in EXPS])
tf_arr = np.array([TF_OF[e] for e in EXPS])
hit_b10, hit_b5 = fr_b <= .1, fr_b <= .05

# ---------------------------------------------------------------- statistics
def mcnemar(h1, h0):
    b, c = int(np.sum(h1 & ~h0)), int(np.sum(~h1 & h0))
    p = 1.0 if b + c == 0 else stats.binomtest(b, b + c, 0.5).pvalue
    return b, c, p

def boot_ci(h1, h0, clusters, n=4000, seed=0):
    """TF-cluster bootstrap 95% CI for the top-10% difference (variant minus baseline)."""
    rng = np.random.default_rng(seed)
    u = np.unique(clusters); g = [np.flatnonzero(clusters == t) for t in u]
    d = np.empty(n)
    for i in range(n):
        idx = np.concatenate([g[k] for k in rng.integers(len(u), size=len(u))])
        d[i] = h1[idx].mean() - h0[idx].mean()
    return np.percentile(d, [2.5, 97.5])

def report(name, fr_v_dict):
    """one table row; paired stats on the common scoreable set, plus the all-279 count."""
    common = [e for e in EXPS if e in fr_v_dict]
    lost_exp = [e for e in EXPS if e not in fr_v_dict]
    f1 = np.array([fr_v_dict[e] for e in common])
    m = np.array([e in fr_v_dict for e in EXPS])
    f0 = fr_b[m]; cl = tf_arr[m]
    h1_10, h0_10, h1_5, h0_5 = f1 <= .1, f0 <= .1, f1 <= .05, f0 <= .05
    b, c, p = mcnemar(h1_10, h0_10)
    lo, hi = boot_ci(h1_10.astype(float), h0_10.astype(float), cl)
    print(f"  {name:20s} {len(common):4d} {len(lost_exp):4d}  "
          f"{h1_10.sum():3d}/279  {h1_10.sum():3d}/{len(common):<4d} {h1_10.mean():.3f}  "
          f"{h1_5.sum():3d}/{len(common):<4d} {h1_5.mean():.3f}  "
          f"{np.mean(f1 < f0):.2f} {np.mean(f1 > f0):.2f}  {b:3d} {c:3d} {p:6.3f}  "
          f"{h1_10.mean()-h0_10.mean():+.3f} [{lo:+.3f},{hi:+.3f}]")
    return dict(name=name, n=len(common), unscoreable=lost_exp, top10=int(h1_10.sum()),
                top5=int(h1_5.sum()), gained=b, lost=c, p=p, ci=(lo, hi),
                diff=h1_10.mean() - h0_10.mean(), fr=fr_v_dict)

HDR = (f"  {'variant':20s} {'n':>4s} {'uns':>4s}  {'top10/279':>9s}  {'top10 common':>15s}  "
       f"{'top5 common':>14s}  {'bet':>4s} {'wor':>4s}  {'gn':>3s} {'ls':>3s} {'McNem':>6s}  "
       f"{'d(top10) [TF-cluster 95% CI]':>28s}")

print(f"\n{'='*140}\nPART A -- sign interventions on MASKED ULM (own gene zeroed throughout), "
      f"re-pruned to tmin=5 where edges are removed")
print(f"  network: {len(net)} edges, {N_NEG} at -1, {int(DEF.sum())} default activation; "
      f"masked-ULM baseline {hit_b10.sum()}/279 top-10%, {hit_b5.sum()}/279 top-5%")
print(f"  'uns' = experiments that become unscoreable (TF drops below tmin=5); 'top10/279' counts "
      f"those as misses; everything right of it is the PAIRED common scoreable set.")
print(f"  gn/ls = hits gained/lost vs masked ULM, McNem = exact two-sided McNemar on those.\n")
print(HDR)
res = {}
res["masked ULM"] = report("0 masked ULM", base)

# --- 1 signs erased -------------------------------------------------------
res["sign-erased"] = report("1 sign-erased", evaluate(build_score(net.assign(weight=1.0))))
# --- 2 default edges removed ----------------------------------------------
res["no-default"] = report("2 no-default", evaluate(build_score(net[~DEF])))
# --- 3/4 the PMID-majority-repressing TFs ---------------------------------
pm = net[net.sign_decision == "PMID"]
cnt = pm.groupby("source").weight.agg(pos=lambda s: (s > 0).sum(), neg=lambda s: (s < 0).sum())
REP = set(cnt.index[cnt.neg > cnt.pos])
rep_def = DEF & net.source.isin(REP)
print(f"\n  [3/4] {len(REP)} TFs have PMID-majority-repressing signs; "
      f"{net.loc[rep_def, 'source'].nunique()} of them also carry {int(rep_def.sum())} "
      f"default-activation edges.\n")
flipped = net.copy(); flipped.loc[rep_def, "weight"] = -1.0
res["flip-default-rep"] = report("3 flip-default-rep", evaluate(build_score(flipped)))
res["drop-default-rep"] = report("4 drop-default-rep", evaluate(build_score(net[~rep_def])))

# --- 5a random sign flips of the same count as the -1 set ------------------
print()
def multi(label, nets, key):
    rows = [report(f"{label} s{s}", evaluate(build_score(n2))) for s, n2 in enumerate(nets)]
    t10 = np.array([r["top10"] for r in rows]); t5 = np.array([r["top5"] for r in rows])
    print(f"  {label+' MEAN':20s} {'':4s} {'':4s}  {t10.mean():5.1f}/279  "
          f"[{t10.min()}-{t10.max()}] top10, [{t5.min()}-{t5.max()}] top5 over {len(rows)} seeds")
    res[key] = rows
    return rows

nets_flip = []
for s in range(5):
    rng = np.random.default_rng(100 + s)
    w = net.weight.values.copy()
    w[rng.choice(len(net), size=N_NEG, replace=False)] *= -1
    nets_flip.append(net.assign(weight=w))
rows_flip = multi("5a rand-flip", nets_flip, "rand-flip")

print()
def_per_tf = net[DEF].groupby("source").size()
nets_drop, n_dropped = [], []
for s in range(5):
    rng = np.random.default_rng(200 + s)
    keep = np.ones(len(net), bool)
    for src, idx in net.groupby("source").indices.items():
        k = int(def_per_tf.get(src, 0))
        if k: keep[rng.choice(idx, size=k, replace=False)] = False
    n_dropped.append(int((~keep).sum()))
    nets_drop.append(net[keep])
rows_drop = multi("5b rand-drop-perTF", nets_drop, "rand-drop")

for r in [res["sign-erased"], res["no-default"], res["drop-default-rep"]]:
    if r["unscoreable"]:
        print(f"\n  unscoreable under {r['name']}: {len(r['unscoreable'])} experiments, TFs "
              f"{sorted({TF_OF[e] for e in r['unscoreable']})}; "
              f"{sum(hit_b10[EXPS.index(e)] for e in r['unscoreable'])} of them were masked-ULM hits")

# ---------------------------------------------------------------- shuffled-network gate
gain = [r for r in [res["sign-erased"], res["no-default"], res["flip-default-rep"],
                    res["drop-default-rep"]] if r["ci"][0] > 0]
print(f"\n  shuffled-network control: ", end="")
if gain:
    from controls import shuffle_preserving_degree
    print(f"{[r['name'] for r in gain]} beat masked ULM with a CI excluding 0; adding the "
          f"degree-preserving null for interventions 1 and 2")
    for lab, n2 in [("1 sign-erased SHUF", net.assign(weight=1.0)), ("2 no-default SHUF", net[~DEF])]:
        pn = dc.pp.prune(features=var, net=n2, tmin=5)
        s, _f, a = dc.pp.adjmat(features=var, net=pn)
        B = shuffle_preserving_degree(pd.DataFrame(a), n_swaps=200000, seed=0).values
        report(lab, evaluate(pd.DataFrame(dc.mt.ulm.func(mat_m, B)[0], index=obs, columns=s)))
else:
    print("NO intervention exceeds masked ULM's 96/279 by more than its bootstrap CI, so no gain "
          "is claimed and no shuffled-network row is run (it would be a null for a null result).")

# ================================================================= PART B
print(f"\n{'='*140}\nPART B -- pre-stated test of V-017's post hoc lead: PMID-majority-REPRESSING "
      f"TFs with default edges rank worse under masked ULM")
has_def = set(net.loc[DEF, "source"])
def cls(tf):
    if tf not in has_def: return "no default edges"
    if tf in REP: return "repressor-majority"
    p, n_ = (cnt.loc[tf, "pos"], cnt.loc[tf, "neg"]) if tf in cnt.index else (0, 0)
    return "activator-majority" if p > n_ else "tie / no PMID"
lab = np.array([cls(t) for t in tf_arr])
reg_size = np.array([(adjm_orig[:, src_i[t]] != 0).sum() for t in tf_arr])
own_lfc = np.array([mat[obs_i[e], var_i[t]] if t in var_i else np.nan for e, t in exps])

print(f"\n  {'class':22s} {'exps':>5s} {'TFs':>4s} {'med.size':>9s} {'top10':>11s} {'top5':>11s} "
      f"{'med.rank':>9s} {'med own logFC':>14s}")
for g in ["repressor-majority", "activator-majority", "no default edges", "tie / no PMID"]:
    m = lab == g
    print(f"  {g:22s} {m.sum():5d} {len(set(tf_arr[m])):4d} {np.median(reg_size[m]):9.0f} "
          f"{hit_b10[m].sum():3d}/{m.sum():<3d} {hit_b10[m].mean():.3f} "
          f"{hit_b5[m].sum():3d}/{m.sum():<3d} {hit_b5[m].mean():.3f} {np.median(fr_b[m]):9.3f} "
          f"{np.nanmedian(own_lfc[m]):14.3f}")

# TF-level table
tfl = pd.DataFrame({"tf": tf_arr, "fr": fr_b, "hit": hit_b10, "size": reg_size,
                    "lfc": own_lfc, "cls": lab}).groupby("tf").agg(
    cls=("cls", "first"), size=("size", "first"), n=("fr", "size"),
    med_fr=("fr", "median"), mean_fr=("fr", "mean"), hits=("hit", "sum"), med_lfc=("lfc", "median"))
R, A = tfl[tfl.cls == "repressor-majority"], tfl[tfl.cls == "activator-majority"]
print(f"\n  TF-level ({len(R)} repressor-majority vs {len(A)} activator-majority TFs):")
print(f"    median of per-TF median masked rank fraction : {R.med_fr.median():.3f} vs {A.med_fr.median():.3f}")
print(f"    mean   of per-TF mean   masked rank fraction : {R.mean_fr.mean():.3f} vs {A.mean_fr.mean():.3f}")
print(f"    median regulon size                          : {R['size'].median():.0f} vs {A['size'].median():.0f}")
print(f"    median own logFC (knockdown efficiency)      : {R.med_lfc.median():.3f} vs {A.med_lfc.median():.3f}"
      f"   Mann-Whitney p = {stats.mannwhitneyu(R.med_lfc.dropna(), A.med_lfc.dropna()).pvalue:.3f}")
print(f"    Mann-Whitney on per-TF median rank fraction  : p = "
      f"{stats.mannwhitneyu(R.med_fr, A.med_fr).pvalue:.3f}")

SB = [("<=60", None), ("61-150", None), (">150", None)]
sz = tfl["size"].values
bin_of = np.where(sz <= 60, "<=60", np.where(sz <= 150, "61-150", ">150"))
two = tfl.cls.isin(["repressor-majority", "activator-majority"]).values

def perm_within_bin(stat, n=20000, seed=0):
    """permute the repressor/activator label among TFs within regulon-size bin; two-sided p."""
    rng = np.random.default_rng(seed)
    isR = (tfl.cls.values == "repressor-majority")[two]
    b = bin_of[two]
    obs_stat = stat(isR)
    null = np.empty(n)
    for i in range(n):
        p = isR.copy()
        for bb in np.unique(b):
            k = b == bb
            p[k] = rng.permutation(isR[k])
        null[i] = stat(p)
    return obs_stat, float(np.mean(np.abs(null) >= abs(obs_stat) - 1e-12)), null

tf_medfr = tfl.med_fr.values[two]
primary = lambda isR: np.median(tf_medfr[isR]) - np.median(tf_medfr[~isR])
o1, p1, _ = perm_within_bin(primary)
tf_of_exp = {t: i for i, t in enumerate(tfl.index)}
exp_isR_src = np.array([tf_of_exp[t] for t in tf_arr])
keep_exp = two[exp_isR_src]
def secondary(isR):
    memb = np.zeros(len(tfl), bool); memb[np.flatnonzero(two)] = isR
    e = memb[exp_isR_src][keep_exp]; h = hit_b10[keep_exp]
    return h[e].mean() - h[~e].mean()
o2, p2, _ = perm_within_bin(secondary)
print(f"\n  PRIMARY  (pre-stated): TF-level median masked rank fraction, repressor minus activator")
print(f"    observed {o1:+.3f} (medians {R.med_fr.median():.3f} vs {A.med_fr.median():.3f}); "
      f"permutation within size bin, two-sided p = {p1:.4f}")
print(f"  SECONDARY: experiment-level top-10% rate, repressor minus activator")
print(f"    observed {o2:+.3f} ({hit_b10[lab=='repressor-majority'].sum()}/{(lab=='repressor-majority').sum()}"
      f" vs {hit_b10[lab=='activator-majority'].sum()}/{(lab=='activator-majority').sum()}); "
      f"permutation within size bin, two-sided p = {p2:.4f}")
print(f"\n  per size bin (repressor vs activator, top-10% hits/exps):")
for bl, _ in SB:
    cells = []
    for g in ["repressor-majority", "activator-majority"]:
        m = (lab == g) & np.isin(tf_arr, tfl.index[(bin_of == bl) & (tfl.cls.values == g)])
        cells.append(f"{hit_b10[m].sum():3d}/{m.sum():<4d} ({len(set(tf_arr[m])):2d} TFs, "
                     f"med rank {np.median(fr_b[m]) if m.sum() else np.nan:.3f})")
    print(f"    {bl:8s} {cells[0]:34s} {cells[1]:34s}")

# --- mechanistic check: weight mass and PMID-only scoring ------------------
scored = pnet[pnet.source.isin(src_i)].merge(net[["source", "target", "sign_decision"]],
                                             on=["source", "target"], how="left")
scored["bucket"] = scored.sign_decision + " " + np.where(scored.weight > 0, "+1", "-1")
print(f"\n  mechanistic: composition of the SCORED regulon (|weight| is 1 for every CollecTRI edge, "
      f"so weight mass == edge count)")
print(f"    {'class':22s} " + " ".join(f"{b:>22s}" for b in
      ["default activation +1", "PMID +1", "PMID -1", "regulon +1/-1"]))
for g in ["repressor-majority", "activator-majority"]:
    srcs = set(tfl.index[tfl.cls == g])
    sub = scored[scored.source.isin(srcs)]
    frac = sub.groupby("source").bucket.value_counts(normalize=True).unstack(fill_value=0)
    for b in ["default activation +1", "PMID +1", "PMID -1", "regulon +1", "regulon -1"]:
        if b not in frac: frac[b] = 0.0
    print(f"    {g:22s} " + " ".join(f"{v:>22.3f}" for v in
          [frac["default activation +1"].median(), frac["PMID +1"].median(),
           frac["PMID -1"].median(), (frac["regulon +1"] + frac["regulon -1"]).median()]))

fr_pmid = evaluate(build_score(net[net.sign_decision == "PMID"]))
print(f"\n  PMID-signed edges only (re-pruned to tmin=5), masked, on the COMMON scoreable set:")
print(f"    {'class':22s} {'n common':>9s} {'masked ULM':>12s} {'PMID-only':>12s} {'gained':>7s} {'lost':>6s}")
for g in ["repressor-majority", "activator-majority", "no default edges", "tie / no PMID", "ALL"]:
    m = np.ones(len(EXPS), bool) if g == "ALL" else (lab == g)
    m = m & np.array([e in fr_pmid for e in EXPS])
    h1 = np.array([fr_pmid[e] for e in np.array(EXPS)[m]]) <= .1
    h0 = hit_b10[m]
    print(f"    {g:22s} {m.sum():9d} {h0.sum():5d}/{m.sum():<6d} {h1.sum():5d}/{m.sum():<6d} "
          f"{int(np.sum(h1&~h0)):7d} {int(np.sum(~h1&h0)):6d}")

# F-029: persist. Per-experiment masked-ULM baseline with the PMID/default class labels, and the
# per-TF aggregate the repressor-vs-activator block prints. No printed line changes.
pd.DataFrame(dict(exp=EXPS, tf=tf_arr, masked_frac=fr_b, hit10=hit_b10, hit5=hit_b5,
                  regulon_size=reg_size, own_logfc=own_lfc, sign_class=lab)
             ).to_csv("../results/signprobe.csv", index=False)
tfl.to_csv("../results/signprobe_tf.csv")

assert (len(EXPS) == 279 and hit_b10.sum() == 96 and hit_b5.sum() == 78
        and np.abs(fr_b - leak.evaluate(leak.S)).max() < 1e-9
        and (mat_m[[obs_i[e] for e, t in exps], [var_i[t] for e, t in exps]] == 0).all()
        and not (net[~DEF].sign_decision == "default activation").any()
        and all(d == int(DEF.sum()) for d in n_dropped)), \
    "masked-ULM baseline does not reproduce leak.py's 96/78 of 279, the own gene is not masked, " \
    "or a random-removal control did not remove exactly as many edges as the default set"
