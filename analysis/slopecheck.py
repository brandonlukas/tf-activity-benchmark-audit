"""T-030 / F-018: is the own-logFC slope of F-015 identifiable in knockTF, or is it a constant offset?

V-015 ran this as a scratch PREVIEW and refused to verify F-015 on the strength of it. This is
the committed version. The quantity under test is the premise of README's "Where this points":
"a slope calibrated on knockdown data predicts the abundance-driven part of activity".

  a_te = self-edge-dropped ULM score of TF t in experiment e  (selfedge.py's es_ns, RAW, not
         type_p x a: type_p = -1 for every knockTF row, so multiplying flips the sign convention
         V-015 fixed -- b = a/m > 0 means "knockdown lowers activity")
  m_te = logFC of gene t in experiment e, read from the scored matrix

CONTROLS. The masked control is degenerate here: the input is already the self-edge-dropped
score, and F-010 says self-edge deletion IS masking for ULM (96/279 both), so zeroing gene t
moves a_te by a median 0.002 t (V-015 (B)). No rank gain is claimed anywhere in this script --
there is no top-10% number to inflate -- so the degree-preserving shuffle has nothing to
protect. The control that bites is the PLACEBO: own logFC permuted across experiments within
regulon-size tercile. Every knockTF m is negative, so a through-origin slope is a rescaled
mean of a and a placebo regressor with the right marginal distribution recovers most of it.

PREDICTIONS (stated before running, from V-015's preview):
  (1) corr(a, m) = 0.10; OLS with intercept: slope 0.41, intercept -0.90;
      within-TF slope (TF fixed intercepts, 176 experiments) 0.08 +/- 0.34.
  (2) through-origin pooled slope: real mu 0.82, placebo mu 0.75 (5-95%: 0.69-0.83).
      In-sample RMSE ratio vs a = 0: real 0.905, placebo 0.919, constant mean(a) 0.904.
  (3) LOO RMSE on the 176: null 4.50, global slope 3.98, global mean 3.96, per-TF slope 3.19,
      per-TF mean 2.74 (median |error| 1.62 / 1.66 / 1.73 / 1.65 / 1.39).
  (4) negative raw per-TF slopes: 50 of 155 (32.3%), against a model-implied expectation of
      about 0.318 -- i.e. F-015's kill condition K4 ("more than 30% negative") is miscalibrated,
      it fires at the value the model itself predicts.
  (5) m is negative in all 279 (median about -1.74) and few TFs have any within-TF dose range.

CLAIM UNDER TEST: "in knockTF the own-logFC slope is not distinguishable from a constant
knockdown offset". WHAT WOULD OVERTURN IT: the within-TF slope CI (TF-cluster bootstrap)
excluding 0, or the per-TF slope beating the per-TF mean out of sample with a CI excluding 0.
Either would mean the dose-response is identified and F-015's b_t x m term earns its place.

JUDGMENT CALLS: (a) LOO folds re-estimate the global quantities on the other 278 experiments
(singletons inform the global fit, as F-015 designs) and the per-TF quantities on that TF's
remaining experiments; (b) regulon size = nonzero column count of the self-edge-dropped pruned
adjacency, the network actually scored; (c) 2000 TF-cluster bootstrap draws, seed 0.
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/slopecheck.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import numpy as np, pandas as pd
from scipy.stats import norm
from selfedge import es_ns, adata, mat, var_i, obs_i, src_i, adjm_ns, exps  # reruns propagate/controls/leak/selfedge

rng = np.random.default_rng(0)
NB = 2000
obs_logfc = adata.obs["logFC"] if "logFC" in adata.obs else None

tf_e = np.array([t for _, t in exps])
a = np.array([es_ns.loc[e, t] for e, t in exps])                       # raw self-edge-dropped score
m = np.array([mat[obs_i[e], var_i[t]] if t in var_i else np.nan for e, t in exps])
n_missing = int(np.isnan(m).sum())
if n_missing:                                                          # fall back to knockTF's own field
    for i, (e, t) in enumerate(exps):
        if np.isnan(m[i]): m[i] = obs_logfc.loc[e]
size = np.array([(adjm_ns[:, src_i[t]] != 0).sum() for _, t in exps])  # regulon size, self-edge already gone

tfs_u = np.array(sorted(set(tf_e)))
idx = {t: np.flatnonzero(tf_e == t) for t in tfs_u}
multi = np.array([t for t in tfs_u if len(idx[t]) >= 2])
in176 = np.isin(tf_e, multi)

def thr(mm, aa):                      # through-origin slope
    return float((mm * aa).sum() / (mm * mm).sum())

print(f"\n{'='*78}\nT-030 / F-018  own-logFC slope vs constant offset, self-edge-dropped ULM")
print(f"{len(a)} experiments, {len(tfs_u)} TFs; {len(multi)} TFs with >= 2 experiments "
      f"= {in176.sum()} experiments; {len(tfs_u) - len(multi)} singletons; "
      f"{n_missing} own-logFC values taken from obs instead of the matrix")

# ---------------------------------------------------------------- (1) correlation and slopes
r = float(np.corrcoef(a, m)[0, 1])
sl, ic = np.polyfit(m, a, 1)
am = a[in176] - pd.Series(a[in176]).groupby(tf_e[in176]).transform("mean").values   # TF fixed intercepts
mm = m[in176] - pd.Series(m[in176]).groupby(tf_e[in176]).transform("mean").values
within = thr(mm, am)

def boot_tfs(pool, stat):
    """TF-cluster bootstrap: resample TFs with replacement, recompute stat on the pooled rows."""
    out = []
    for _ in range(NB):
        draw = rng.choice(pool, size=len(pool), replace=True)
        rows = np.concatenate([idx[t] for t in draw])
        out.append(stat(rows, draw))
    return np.array(out)

def within_stat(rows, draw):
    aa, mmm = a[rows], m[rows]
    lab = np.concatenate([[i] * len(idx[t]) for i, t in enumerate(draw)])           # per-DRAW group, not per-TF
    ac = aa - pd.Series(aa).groupby(lab).transform("mean").values
    mc = mmm - pd.Series(mmm).groupby(lab).transform("mean").values
    return thr(mc, ac) if (mc * mc).sum() > 0 else np.nan

wb = boot_tfs(multi, within_stat); wb = wb[~np.isnan(wb)]
lo, hi = np.percentile(wb, [5, 95])
print(f"\n(1) corr(a, m) = {r:.3f}   OLS with intercept: slope {sl:.3f}, intercept {ic:.3f}")
print(f"    within-TF slope (TF fixed intercepts, {in176.sum()} exps, {len(multi)} TFs) = {within:.3f}"
      f"  TF-cluster bootstrap 90% CI [{lo:.2f}, {hi:.2f}]  (sd {wb.std():.2f})")
print(f"    -> CI {'EXCLUDES' if lo > 0 or hi < 0 else 'INCLUDES'} 0")

# ---------------------------------------------------------------- (2) pooled slope vs permuted placebo
terc = pd.qcut(pd.Series(size), 3, labels=False).values
mu_real = thr(m, a)
rmse = lambda err: float(np.sqrt(np.mean(err ** 2)))
r_null = rmse(a)
ratio_real = rmse(a - mu_real * m) / r_null
ratio_const = rmse(a - a.mean()) / r_null
perm_mu, perm_ratio = [], []
for _ in range(500):
    mp = m.copy()
    for g in np.unique(terc):
        j = np.flatnonzero(terc == g); mp[j] = m[rng.permutation(j)]
    b = thr(mp, a); perm_mu.append(b); perm_ratio.append(rmse(a - b * mp) / r_null)
perm_mu, perm_ratio = np.array(perm_mu), np.array(perm_ratio)
print(f"\n(2) through-origin pooled slope, 500 permutations of own logFC within size tercile")
print(f"{'fit':34s} {'slope':>8s} {'in-sample RMSE ratio vs a=0':>30s}")
print(f"{'real m':34s} {mu_real:8.3f} {ratio_real:30.3f}")
print(f"{'placebo (m permuted in tercile)':34s} {perm_mu.mean():8.3f} {perm_ratio.mean():30.3f}")
print(f"{'  5-95% of placebo':34s} {np.percentile(perm_mu,5):.2f}-{np.percentile(perm_mu,95):.2f}  "
      f"{np.percentile(perm_ratio,5):>21.3f}-{np.percentile(perm_ratio,95):.3f}")
print(f"{'constant mean(a) = %.2f' % a.mean():34s} {'--':>8s} {ratio_const:30.3f}")
print(f"    real slope is above {np.mean(perm_mu < mu_real):.1%} of placebos; "
      f"real RMSE ratio below {np.mean(perm_ratio > ratio_real):.1%} of them")

# ---------------------------------------------------------------- (3) leave-one-experiment-out
loo = {k: np.full(len(a), np.nan) for k in
       ["null (0)", "global slope x m", "global mean", "per-TF slope x m", "per-TF mean"]}
for i in np.flatnonzero(in176):
    o = np.ones(len(a), bool); o[i] = False                 # all other 278, singletons included
    t = np.flatnonzero((tf_e == tf_e[i]) & o)               # this TF's remaining experiments
    loo["null (0)"][i] = 0.0
    loo["global slope x m"][i] = thr(m[o], a[o]) * m[i]
    loo["global mean"][i] = a[o].mean()
    loo["per-TF slope x m"][i] = thr(m[t], a[t]) * m[i]
    loo["per-TF mean"][i] = a[t].mean()
err = {k: a - v for k, v in loo.items()}
no_rest = in176 & (tf_e != "REST")
print(f"\n(3) leave-one-experiment-out on the {in176.sum()} multi-experiment experiments "
      f"(everything re-estimated inside the fold)")
print(f"{'predictor of held-out a':22s} {'RMSE':>7s} {'med|err|':>9s} | {'RMSE-noREST':>11s} {'med|err|':>9s}")
for k, v in err.items():
    print(f"{k:22s} {rmse(v[in176]):7.2f} {np.median(np.abs(v[in176])):9.2f} | "
          f"{rmse(v[no_rest]):11.2f} {np.median(np.abs(v[no_rest])):9.2f}")
rest_share = (err["null (0)"][in176 & (tf_e == 'REST')] ** 2).sum() / (err["null (0)"][in176] ** 2).sum()
print(f"    REST: {(in176 & (tf_e=='REST')).sum()} experiments, {rest_share:.0%} of the null's squared error")

def dstat(rows, draw, keep=None):
    rr = rows if keep is None else rows[keep[rows]]
    if len(rr) < 2: return np.nan
    return rmse(err["per-TF slope x m"][rr]) - rmse(err["per-TF mean"][rr])
d_all = rmse(err["per-TF slope x m"][in176]) - rmse(err["per-TF mean"][in176])
d_nr = rmse(err["per-TF slope x m"][no_rest]) - rmse(err["per-TF mean"][no_rest])
b_all = boot_tfs(multi, dstat)
b_nr = boot_tfs(multi, lambda rows, draw: dstat(rows, draw, no_rest)); b_nr = b_nr[~np.isnan(b_nr)]
for lab, d, b in [("all 176", d_all, b_all), ("without REST", d_nr, b_nr)]:
    l, h = np.percentile(b, [5, 95])
    print(f"    RMSE(per-TF slope) - RMSE(per-TF mean), {lab:12s} = {d:+.2f}  "
          f"90% TF-cluster CI [{l:+.2f}, {h:+.2f}]  -> {'EXCLUDES' if l > 0 or h < 0 else 'INCLUDES'} 0"
          f"{', slope WORSE' if d > 0 else ', slope better'}")

# ---------------------------------------------------------------- (4) negative slopes vs model expectation
b_t = np.array([thr(m[idx[t]], a[idx[t]]) for t in tfs_u])
res = np.concatenate([a[idx[t]] - b_t[i] * m[idx[t]] for i, t in enumerate(tfs_u)])
sigma2 = float((res ** 2).sum() / (len(a) - len(tfs_u)))                  # dof 279 - 155
s2 = np.array([sigma2 / (m[idx[t]] ** 2).sum() for t in tfs_u])
w = 1 / s2
mu_fe = float((w * b_t).sum() / w.sum())
Q = float((w * (b_t - mu_fe) ** 2).sum())
tau2 = max(0.0, (Q - (len(tfs_u) - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))   # DerSimonian-Laird
p_neg = norm.cdf(-mu_fe / np.sqrt(tau2 + s2))
print(f"\n(4) raw per-TF through-origin slopes, {len(tfs_u)} TFs")
print(f"    negative: {(b_t < 0).sum()} of {len(tfs_u)} = {np.mean(b_t < 0):.3f}   median b_t {np.median(b_t):.2f}")
print(f"    model-implied expectation at mu = {mu_fe:.3f}, tau = {np.sqrt(tau2):.2f}, sigma = {np.sqrt(sigma2):.2f}: "
      f"{p_neg.mean():.3f} +/- {np.sqrt((p_neg*(1-p_neg)).sum())/len(tfs_u):.3f}")
print(f"    F-015's K4 line (30%) sits {(0.30 - p_neg.mean())/ (np.sqrt((p_neg*(1-p_neg)).sum())/len(tfs_u)):+.1f} sd "
      f"from what the fitted model itself predicts -> K4 is a miscalibrated threshold, not a test")

# ---------------------------------------------------------------- (5) dose range
q = np.percentile(m, [0, 25, 50, 75, 100])
rng_tf = np.array([m[idx[t]].max() - m[idx[t]].min() for t in multi])
print(f"\n(5) dose. m over {len(m)} experiments: min {q[0]:.2f}  q25 {q[1]:.2f}  median {q[2]:.2f}  "
      f"q75 {q[3]:.2f}  max {q[4]:.2f}; negative in {(m < 0).sum()}/{len(m)}")
print(f"    within-TF range of m over the {len(multi)} multi-experiment TFs: "
      f"min {rng_tf.min():.2f}  median {np.median(rng_tf):.2f}  max {rng_tf.max():.2f}")
print(f"    TFs with within-TF range > 1 logFC: {(rng_tf > 1).sum()} of {len(multi)} "
      f"({(rng_tf > 1).sum()} of {len(tfs_u)} TFs overall); > 2 logFC: {(rng_tf > 2).sum()}")
print(f"    a slope needs within-TF variation in dose; {(rng_tf <= 1).sum()} of {len(multi)} TFs have under 1 logFC of it")
print(f"{'='*78}")

# F-029: persist. One row per experiment: the self-edge-dropped score a, the own logFC dose m,
# regulon size and the TF label -- the inputs to every slope in this script. Prints unchanged.
pd.DataFrame(dict(exp=[e for e, _ in exps], tf=tf_e, a=a, m=m, regulon_size=size,
                  multi_tf=in176)).to_csv("../results/slopecheck.csv", index=False)

# through-origin algebra: a singleton TF's residual is identically 0, and every knockdown is a knockdown
sing = [t for t in tfs_u if len(idx[t]) == 1]
assert (len(a) == 279 and len(tfs_u) == 155 and in176.sum() == 176 and (m < 0).all()
        and max(abs(a[idx[t]][0] - b_t[list(tfs_u).index(t)] * m[idx[t]][0]) for t in sing) < 1e-9
        and np.isfinite([r, sl, within, mu_real, tau2]).all()), \
    "design broke: wrong counts, a non-negative own logFC, or the through-origin fit is not exact on singletons"
