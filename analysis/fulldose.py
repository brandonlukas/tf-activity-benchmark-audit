"""T-034 / F-022: the slope check and the benchmark, on the UNFILTERED knockTF (dose-resolved).

V-018 said knockTF cannot calibrate an abundance->activity slope because decoupler's default
`dc.ds.knocktf(thr_fc=-1)` keeps only experiments whose TF's own logFC is below -1: every
experiment is a deep knockdown, no dose range, no data where the method would be applied.
`dc.ds.knocktf(thr_fc=None)` is the same source with the filter off: data/knocktf_full.h5ad,
907 experiments over 456 TFs, own logFC from -8.5 to +9.9, 379 in (-1, 0), 140 at or above 0,
184 TFs with >= 2 experiments, 76 of them with more than 1 logFC of within-TF range.

  a_te = self-edge-dropped ULM score of TF t in experiment e (RAW, not type_p x a; type_p = -1
         everywhere, so b = a/m > 0 reads "knockdown lowers activity", the V-015 convention)
  m_te = logFC of gene t in experiment e, read from the scored matrix (checked against obs.logFC)

Built the same way selfedge.py builds it (dc.pp.extract / prune tmin=5 / adjmat, every TF's
self-edge zeroed, dc.mt.ulm.func) but on the full file, which selfedge.py cannot load. The var
axis of the two files is identical (21985 genes), and ULM scores rows independently, so the
logFC < -1 subset of this pipeline must reproduce slopecheck.py exactly: 279 scored, 155 TFs,
corr 0.095, OLS slope 0.413 / intercept -0.901, within-TF slope 0.079. The final assert checks it.

PREDICTIONS (the main session's, stated blind before running):
  A1  within-TF slope 0.45 t-units per logFC, 95% TF-cluster CI excluding 0
  A1  OLS intercept -0.4 (a dose-response through the origin predicts ~0, a pure knockdown
      offset predicts ~-1.7; -0.4 sits between the two models)
  A2  mean a by own-logFC bin (< -2 / -2..-1 / -1..-0.5 / -0.5..0 / >= 0):
      -2.3 / -1.6 / -0.9 / -0.4 / -0.1, i.e. monotone and through ~0 at m >= 0
  A4  (per-TF mean + pooled within-TF slope) beats (per-TF mean) out of sample by 0.05 RMSE,
      CI excluding 0
  B   masked ULM top-10% by the same bins: 0.36 / 0.33 / 0.22 / 0.15 / 0.11, i.e. falling to
      the 10% chance rate as m -> 0
  control  shuffled network: within-TF slope 0.0 +/- 0.1 and a flat binned curve

KILL / SUPPORT (from the task): the direction dies if the within-TF slope CI sits inside
+/- 0.2 on this design AND the binned curve is flat. It survives if the binned curve is
monotone through about 0 at m >= 0 and the within-TF slope CI excludes 0 and stays out under
the Knock.Method adjustment and the leverage check.

CONTROLS. Masked-gene: built in twice over. Part A's input is the self-edge-dropped score, and
F-010 showed self-edge deletion IS own-gene masking for ULM; Part B scores the masked matrix
the way leak.py does. Degree-preserving shuffle: the TF-TF swapper in controls.py does not
apply to a TF-gene bipartite graph, so the null here redraws each TF's targets at random from
the measured genes (excluding its own) keeping that TF's count of +1 and of -1 targets, 10
seeds. That destroys which genes a TF regulates while preserving regulon size and sign mix, so
a regulon-specific dose-response must vanish under it while a global "deep knockdowns perturb
everything, and ULM t-values pick that up" artefact would survive. It does NOT preserve each
gene's in-degree (a gene's number of regulators), so a hub-gene artefact is not controlled.

JUDGMENT CALLS: (a) bins are the five given in the task, closed on the left; (b) TF-cluster
bootstraps resample TFs with replacement over ALL scored experiments and recompute the
statistic on the drawn rows, so a TF spanning two bins keeps its rows together; 2000 draws,
seed 0, 95% CIs (slopecheck used 90%); (c) LOO folds re-estimate every quantity inside the
fold, global ones on the other 906 experiments (singletons included), per-TF ones on that TF's
remaining experiments; (d) regulon size = nonzero column count of the self-edge-dropped pruned
adjacency; (e) the Knock.Method adjustment is TF dummies + method dummies + m by lstsq.
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/fulldose.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd, scipy.sparse as sps

rng = np.random.default_rng(0)
NB = 2000
BIN_EDGES = [-2.0, -1.0, -0.5, 0.0]
BIN_NAMES = ["m < -2", "-2 <= m < -1", "-1 <= m < -0.5", "-0.5 <= m < 0", "m >= 0"]

# ------------------------------------------------------------------ load, exactly as selfedge.py
adata = ad.read_h5ad("../data/knocktf_full.h5ad")
net = pd.read_parquet("../data/collectri.parquet")
mat, obs, var = dc.pp.extract(adata)
if sps.issparse(mat): mat = mat.toarray()
pnet = dc.pp.prune(features=var, net=net, tmin=5)
sources, feats, adjm = dc.pp.adjmat(features=var, net=pnet)
src_i = {s: i for i, s in enumerate(sources)}
var_i = {v: i for i, v in enumerate(var)}
obs_i = {o: i for i, o in enumerate(obs)}

def score(mm, adj):
    return dc.mt.ulm.func(mm, adj)[0]

adjm_ns = adjm.copy()                                   # every TF's self-edge deleted
self_rc = [(var_i[s], src_i[s]) for s in sources if s in var_i and adjm[var_i[s], src_i[s]] != 0]
for r, c in self_rc: adjm_ns[r, c] = 0
es_ns = pd.DataFrame(score(mat, adjm_ns), index=obs, columns=sources)

exps = [(e, t) for e, t in adata.obs.source.items() if t in src_i]
tf_e = np.array([t for _, t in exps])
a = np.array([es_ns.loc[e, t] for e, t in exps])
m = np.array([mat[obs_i[e], var_i[t]] if t in var_i else np.nan for e, t in exps])
m_obs = np.array([adata.obs.loc[e, "logFC"] for e, _ in exps], float)
n_nomat = int(np.isnan(m).sum())
m_gap = float(np.nanmax(np.abs(m - m_obs)))             # matrix logFC vs obs.logFC
m = np.where(np.isnan(m), m_obs, m)
size = np.array([(adjm_ns[:, src_i[t]] != 0).sum() for _, t in exps])
method = adata.obs.loc[[e for e, _ in exps], "Knock.Method"].values
platform = adata.obs.loc[[e for e, _ in exps], "Platform"].values

tfs_u = np.array(sorted(set(tf_e)))
code = pd.Categorical(tf_e, categories=tfs_u).codes
idx = {t: np.flatnonzero(tf_e == t) for t in tfs_u}
multi = np.array([t for t in tfs_u if len(idx[t]) >= 2])
inm = np.isin(tf_e, multi)
rng_tf = {t: m[idx[t]].max() - m[idx[t]].min() for t in multi}
bin_e = np.digitize(m, BIN_EDGES)

rmse = lambda err: float(np.sqrt(np.mean(err ** 2)))
def thr(mm, aa):                                        # through-origin slope
    return float((mm * aa).sum() / (mm * mm).sum()) if (mm * mm).sum() > 0 else np.nan

def within_slope(aa, mm, lab):
    """OLS slope of aa on mm with one free intercept per group in lab (groups of 1 drop out)."""
    u, inv = np.unique(lab, return_inverse=True)
    cnt = np.bincount(inv)
    ac = aa - np.bincount(inv, weights=aa)[inv] / cnt[inv]
    mc = mm - np.bincount(inv, weights=mm)[inv] / cnt[inv]
    return thr(mc, ac)

def boot_rows(pool, stat, seed=0, nb=NB):
    """TF-cluster bootstrap: resample TFs with replacement, stat(rows, group labels)."""
    r = np.random.default_rng(seed)
    out = []
    for _ in range(nb):
        draw = r.choice(pool, size=len(pool), replace=True)
        rows = np.concatenate([idx[t] for t in draw])
        lab = np.concatenate([np.full(len(idx[t]), i) for i, t in enumerate(draw)])
        out.append(stat(rows, lab))
    return np.array(out, float)

def ci(b, lo=2.5, hi=97.5):
    b = b[np.isfinite(b)]
    return float(np.percentile(b, lo)), float(np.percentile(b, hi))

# ------------------------------------------------------------------ reproduction of slopecheck.py
sub = m_obs < -1
a_s, m_s, tf_s = a[sub], m[sub], tf_e[sub]
multi_s = [t for t in sorted(set(tf_s)) if (tf_s == t).sum() >= 2]
in_s = np.isin(tf_s, multi_s)
repro = dict(n=int(sub.sum()), ntf=len(set(tf_s)), n176=int(in_s.sum()),
             corr=float(np.corrcoef(a_s, m_s)[0, 1]))
repro["slope"], repro["icpt"] = [float(v) for v in np.polyfit(m_s, a_s, 1)]
repro["within"] = within_slope(a_s[in_s], m_s[in_s], tf_s[in_s])
ok_repro = (repro["n"] == 279 and repro["ntf"] == 155 and repro["n176"] == 176
            and abs(repro["corr"] - 0.095) < 5e-3 and abs(repro["slope"] - 0.413) < 5e-3
            and abs(repro["icpt"] + 0.901) < 5e-3 and abs(repro["within"] - 0.079) < 5e-3)

print(f"\n{'='*84}\nT-034 / F-022  dose-resolved slope check on the UNFILTERED knockTF")
print(f"{len(adata)} experiments in the file, {adata.obs.source.nunique()} TFs; "
      f"{len(exps)} scored ({len(tfs_u)} TFs in CollecTRI at tmin=5), "
      f"{len(multi)} TFs with >= 2 scored experiments = {inm.sum()} experiments")
print(f"own logFC from the matrix vs obs.logFC: max |diff| {m_gap:.2e}, "
      f"{n_nomat} TFs not in the matrix (obs.logFC used)")
print(f"self-edges in the pruned network: {len(self_rc)} of {len(sources)} TFs")
print(f"REPRODUCTION of slopecheck.py on the logFC < -1 subset: n {repro['n']} (279), "
      f"TFs {repro['ntf']} (155), multi-exps {repro['n176']} (176), corr {repro['corr']:.3f} (0.095), "
      f"OLS slope {repro['slope']:.3f} (0.413) intercept {repro['icpt']:.3f} (-0.901), "
      f"within-TF slope {repro['within']:.3f} (0.079)  -> {'MATCH' if ok_repro else 'MISMATCH'}")

# ------------------------------------------------------------------ A1 correlation and slopes
r_all = float(np.corrcoef(a, m)[0, 1])
sl, ic = [float(v) for v in np.polyfit(m, a, 1)]
ols_b = boot_rows(tfs_u, lambda rows, lab: np.polyfit(m[rows], a[rows], 1))
sl_ci, ic_ci = ci(ols_b[:, 0]), ci(ols_b[:, 1])
w_all = within_slope(a[inm], m[inm], tf_e[inm])
w_stat = lambda rows, lab: within_slope(a[rows], m[rows], lab)
wb = boot_rows(multi, w_stat)
w_ci = ci(wb)
seed_cis = [ci(boot_rows(multi, w_stat, seed=s, nb=500)) for s in range(5)]

print(f"\n(A1) all {len(a)} scored experiments")
print(f"     corr(a, m) = {r_all:.3f}")
print(f"     OLS with intercept: slope {sl:.3f} 95% CI [{sl_ci[0]:.2f}, {sl_ci[1]:.2f}]   "
      f"intercept {ic:.3f} 95% CI [{ic_ci[0]:.2f}, {ic_ci[1]:.2f}]")
print(f"       intercept models: through-origin dose-response predicts ~0, pure knockdown offset "
      f"~{a[m_obs < -1].mean():.2f} (mean a on the deep knockdowns)")
print(f"{'within-TF slope (TF fixed intercepts)':44s} {'slope':>7s} {'95% TF-cluster CI':>20s} {'exps':>6s} {'TFs':>5s}")
for lab, pool in [("all multi-experiment TFs", multi),
                  ("TFs with within-TF range > 1 logFC", np.array([t for t in multi if rng_tf[t] > 1])),
                  ("TFs with within-TF range > 2 logFC", np.array([t for t in multi if rng_tf[t] > 2]))]:
    rows = np.concatenate([idx[t] for t in pool])
    s = within_slope(a[rows], m[rows], tf_e[rows])
    l, h = ci(boot_rows(pool, w_stat))
    print(f"{lab:44s} {s:7.3f} {f'[{l:+.2f}, {h:+.2f}]':>20s} {len(rows):6d} {len(pool):5d}"
          f"  -> {'EXCLUDES' if l > 0 or h < 0 else 'INCLUDES'} 0")
print(f"     CI spread over 5 bootstrap seeds (500 draws each, V-018 caught a seed artefact): "
      + "  ".join(f"[{l:+.2f},{h:+.2f}]" for l, h in seed_cis))

# ------------------------------------------------------------------ A2 binned dose-response
def bin_stat(rows, lab, fn=np.mean):
    return [fn(a[rows][bin_e[rows] == b]) if (bin_e[rows] == b).any() else np.nan for b in range(5)]
bmean_b = boot_rows(tfs_u, lambda rows, lab: bin_stat(rows, lab))
print(f"\n(A2) binned dose-response, self-edge-dropped ULM score of the perturbed TF")
print(f"{'own logFC bin':16s} {'n':>5s} {'TFs':>5s} {'mean a':>8s} {'95% TF-cluster CI':>20s} {'median a':>9s} {'mean m':>8s}")
for b in range(5):
    j = bin_e == b
    l, h = ci(bmean_b[:, b])
    print(f"{BIN_NAMES[b]:16s} {j.sum():5d} {len(set(tf_e[j])):5d} {a[j].mean():8.2f} "
          f"{f'[{l:+.2f}, {h:+.2f}]':>20s} {np.median(a[j]):9.2f} {m[j].mean():8.2f}")

# ------------------------------------------------------------------ A3 placebo for the within-TF slope
terc = pd.qcut(pd.Series(size), 3, labels=False).values
perm = []
for _ in range(500):
    mp = m.copy()
    for g in np.unique(terc):
        j = np.flatnonzero(terc == g); mp[j] = m[rng.permutation(j)]
    perm.append(within_slope(a[inm], mp[inm], tf_e[inm]))
perm = np.array(perm)
print(f"\n(A3) placebo: own logFC permuted across experiments within regulon-size tercile, 500 draws")
print(f"     within-TF slope  real {w_all:+.3f}   placebo mean {perm.mean():+.3f} "
      f"(5-95% {np.percentile(perm,5):+.3f} to {np.percentile(perm,95):+.3f}); "
      f"real above {np.mean(perm < w_all):.1%} of placebos")

# ------------------------------------------------------------------ A4 leave-one-experiment-out
loo = {k: np.full(len(a), np.nan) for k in
       ["null (0)", "global slope x m", "per-TF mean", "per-TF mean + pooled within-slope x dm",
        "per-TF through-origin slope x m"]}
for i in np.flatnonzero(inm):
    o = np.ones(len(a), bool); o[i] = False
    t = np.flatnonzero((tf_e == tf_e[i]) & o)
    om = o & inm
    loo["null (0)"][i] = 0.0
    loo["global slope x m"][i] = thr(m[o], a[o]) * m[i]
    loo["per-TF mean"][i] = a[t].mean()
    loo["per-TF mean + pooled within-slope x dm"][i] = (
        a[t].mean() + within_slope(a[om], m[om], code[om]) * (m[i] - m[t].mean()))
    loo["per-TF through-origin slope x m"][i] = thr(m[t], a[t]) * m[i]
err = {k: a - v for k, v in loo.items()}
print(f"\n(A4) leave-one-experiment-out on the {inm.sum()} multi-experiment experiments "
      f"(everything re-estimated inside the fold)")
print(f"{'predictor of held-out a':42s} {'RMSE':>7s} {'med|err|':>9s} {'folds':>7s}")
for k, v in err.items():                 # a per-TF through-origin slope is undefined when the
    f = np.isfinite(v) & inm             # TF's remaining experiments all have m = 0
    print(f"{k:42s} {rmse(v[f]):7.2f} {np.median(np.abs(v[f])):9.2f} {f.sum():7d}")
K1, K0 = "per-TF mean + pooled within-slope x dm", "per-TF mean"
d = rmse(err[K1][inm]) - rmse(err[K0][inm])
db = boot_rows(multi, lambda rows, lab: rmse(err[K1][rows]) - rmse(err[K0][rows]))
l, h = ci(db)
print(f"     RMSE(mean + slope) - RMSE(mean) = {d:+.3f}  95% TF-cluster CI [{l:+.3f}, {h:+.3f}]"
      f"  -> {'EXCLUDES' if l > 0 or h < 0 else 'INCLUDES'} 0"
      f"{', slope WORSE' if d > 0 else ', slope better'}")

# ------------------------------------------------------------------ A5 confounds and leverage
def eta2(g, x):                      # fraction of variance of x explained by grouping g
    s = pd.Series(x).groupby(pd.Series(g)).transform("mean").values
    return 1 - ((x - s) ** 2).sum() / ((x - x.mean()) ** 2).sum()
def dummies(v):
    u = pd.unique(v); return np.stack([(v == k).astype(float) for k in u], 1)
D_tf, D_me = dummies(tf_e[inm]), dummies(method[inm])
X_fe = np.column_stack([D_tf, m[inm]])
X_me = np.column_stack([D_tf, D_me, m[inm]])
b_fe = float(np.linalg.lstsq(X_fe, a[inm], rcond=None)[0][-1])
b_me = float(np.linalg.lstsq(X_me, a[inm], rcond=None)[0][-1])
print(f"\n(A5) confounds. variance of own logFC explained by Knock.Method {eta2(method, m):.3f}, "
      f"by Platform {eta2(platform, m):.3f} ({len(set(method))} methods, {len(set(platform))} platforms)")
print(f"     within-TF slope, TF dummies only        {b_fe:+.3f}  (lstsq check of the {w_all:+.3f} above)")
print(f"     within-TF slope, TF + Knock.Method dummies {b_me:+.3f}")
comp = pd.crosstab(pd.Series([BIN_NAMES[b] for b in bin_e], name="bin"),
                   pd.Series(method, name="method"), normalize="index").reindex(BIN_NAMES)
print("     dose-bin composition by Knock.Method (row fractions):")
print(comp.round(2).to_string().replace("\n", "\n       "))
ss = sorted(((float(((m[idx[t]] - m[idx[t]].mean()) ** 2).sum()), t) for t in multi), reverse=True)[:5]
drop = {t for _, t in ss}
keep = np.array([t for t in multi if t not in drop])
rows_k = np.concatenate([idx[t] for t in keep])
w_k = within_slope(a[rows_k], m[rows_k], tf_e[rows_k])
l, h = ci(boot_rows(keep, w_stat))
print(f"     leverage: top-5 TFs by within-TF SS(m) = "
      + ", ".join(f"{t} ({v:.0f})" for v, t in ss)
      + f"; they are {sum(v for v, _ in ss)/sum(((m[idx[t]]-m[idx[t]].mean())**2).sum() for t in multi):.0%} of total SS(m)")
print(f"     within-TF slope without them {w_k:+.3f} 95% CI [{l:+.2f}, {h:+.2f}] "
      f"({len(rows_k)} exps, {len(keep)} TFs) -> {'EXCLUDES' if l > 0 or h < 0 else 'INCLUDES'} 0")

# ------------------------------------------------------------------ CONTROL: shuffled network
def shuffle_tf_gene(adj, seed):
    """Degree- and sign-mix-preserving TF-gene null: each TF keeps its number of +1 and of -1
    targets, drawn at random from the measured genes other than its own (so no self-edge)."""
    r = np.random.default_rng(seed)
    B = np.zeros_like(adj)
    for s, c in src_i.items():
        w = adj[:, c][adj[:, c] != 0]
        if len(w) == 0: continue
        pool = np.arange(adj.shape[0])
        if s in var_i: pool = pool[pool != var_i[s]]
        B[r.choice(pool, size=len(w), replace=False), c] = w
    return B

shuf_bins, shuf_w = [], []
for s in range(10):
    es_sh = score(mat, shuffle_tf_gene(adjm_ns, s))
    a_sh = np.array([es_sh[obs_i[e], src_i[t]] for e, t in exps])
    shuf_bins.append([a_sh[bin_e == b].mean() for b in range(5)])
    shuf_w.append(within_slope(a_sh[inm], m[inm], tf_e[inm]))
shuf_bins, shuf_w = np.array(shuf_bins), np.array(shuf_w)
print(f"\n(CONTROL) degree-preserving shuffled TF-gene network, 10 seeds "
      f"(targets redrawn at random, each TF's +1/-1 counts kept)")
print(f"{'':16s} " + " ".join(f"{n:>22s}" for n in BIN_NAMES))
print(f"{'real mean a':16s} " + " ".join(f"{a[bin_e==b].mean():22.2f}" for b in range(5)))
print(f"{'shuffled mean a':16s} " + " ".join(
    f"{shuf_bins[:,b].mean():+.2f} [{shuf_bins[:,b].min():+.2f},{shuf_bins[:,b].max():+.2f}]".rjust(22)
    for b in range(5)))
print(f"     within-TF slope: real {w_all:+.3f}   shuffled mean {shuf_w.mean():+.3f} "
      f"[{shuf_w.min():+.3f}, {shuf_w.max():+.3f}] over 10 seeds")

# ------------------------------------------------------------------ B benchmark by dose bin
mat_m = mat.copy()
for e, t in exps:
    if t in var_i: mat_m[obs_i[e], var_i[t]] = 0
S_un = pd.DataFrame(score(mat, adjm), index=obs, columns=sources)
S_ma = pd.DataFrame(score(mat_m, adjm), index=obs, columns=sources)

def fracs(S):
    out = []
    for e, t in exps:
        s = adata.obs.loc[e, "type_p"] * S.loc[e]
        out.append(s.rank(ascending=False)[t] / S.shape[1])
    return np.array(out)
fr_un, fr_ma = fracs(S_un), fracs(S_ma)

def hit_stat(fr, thr_):
    return lambda rows, lab: [np.mean(fr[rows][bin_e[rows] == b] <= thr_)
                              if (bin_e[rows] == b).any() else np.nan for b in range(5)]
print(f"\n(B) benchmark by dose bin. rank of the true TF among {len(sources)} scored TFs, "
      f"type_p = -1 convention; chance = 0.10 / 0.05")
print(f"{'variant':22s} {'bin':16s} {'n':>5s} {'top10':>7s} {'hits/n':>10s} {'95% CI':>18s} {'top5':>7s} {'95% CI':>18s}")
for name, fr in [("ULM (unmasked)", fr_un), ("ULM own gene MASKED", fr_ma)]:
    b10 = boot_rows(tfs_u, hit_stat(fr, .10)); b5 = boot_rows(tfs_u, hit_stat(fr, .05))
    for b in range(5):
        j = bin_e == b
        l10, h10 = ci(b10[:, b]); l5, h5 = ci(b5[:, b])
        print(f"{name if b == 0 else '':22s} {BIN_NAMES[b]:16s} {j.sum():5d} "
              f"{np.mean(fr[j] <= .1):7.3f} {f'{int((fr[j] <= .1).sum())}/{j.sum()}':>10s} "
              f"{f'[{l10:.3f}, {h10:.3f}]':>18s} {np.mean(fr[j] <= .05):7.3f} "
              f"{f'[{l5:.3f}, {h5:.3f}]':>18s}")
    print(f"{'  all bins':22s} {'':16s} {len(fr):5d} {np.mean(fr <= .1):7.3f} "
          f"{f'{int((fr <= .1).sum())}/{len(fr)}':>10s} {'':>18s} {np.mean(fr <= .05):7.3f}")
print(f"     masked vs unmasked, paired on {len(fr_un)} experiments: "
      f"better {np.mean(fr_ma < fr_un):.2f}  worse {np.mean(fr_ma > fr_un):.2f}  "
      f"tied {np.mean(fr_ma == fr_un):.2f}")
print(f"{'='*84}")

# F-029: persist. One row per scored experiment on the UNFILTERED file: the dose m, the
# self-edge-dropped score a, the dose bin, and both benchmark rank fractions. Prints unchanged.
pd.DataFrame(dict(exp=[e for e, _ in exps], tf=tf_e, a=a, m=m, m_obs=m_obs,
                  regulon_size=size, knock_method=method, platform=platform,
                  dose_bin=[BIN_NAMES[b] for b in bin_e], in_filtered_279=m_obs < -1,
                  multi_tf=inm, ulm_unmasked_frac=fr_un, ulm_masked_frac=fr_ma)
             ).to_csv("../results/fulldose.csv", index=False)

# the reproduction is the load-bearing claim (same pipeline as slopecheck.py), and the shuffle
# must keep every TF's regulon size and sign mix while moving every target off its real gene
B0 = shuffle_tf_gene(adjm_ns, 0)
assert (ok_repro and m_gap < 1e-3
        and ((B0 != 0).sum(0) == (adjm_ns != 0).sum(0)).all()
        and ((B0 > 0).sum(0) == (adjm_ns > 0).sum(0)).all()
        and sum(B0[var_i[s], src_i[s]] != 0 for s in sources if s in var_i) == 0
        and len(a) == len(m) == len(exps) and np.isfinite([r_all, sl, ic, w_all]).all()), \
    "design broke: the logFC<-1 subset no longer reproduces slopecheck.py, the matrix logFC " \
    "disagrees with obs.logFC, or the shuffled network is not degree/sign-preserving"
