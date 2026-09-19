"""T-037 / F-026: is masked ULM detecting the perturbed TF, or the study it came from?

V-022's scratch lead: among knockTF experiments that come from a study which also perturbed a
DIFFERENT TF, masked ULM puts the perturbed TF in the top 10% in 0.277 of experiments and puts
the SAME TF in the top 10% in 0.225 of its study-mates' experiments (155 scored experiments,
45 studies; +0.052, CI including 0), while the same TFs score 0.131 in other studies. If the
paired difference is really ~0, the 0.344 masked baseline is mostly study context (cell type,
batch, treatment) rather than TF identity, which is a worse problem than the abundance leak and
applies to every method scored on knockTF.

Everything here is MASKED ULM: the perturbed TF's own gene is zeroed in its own row before
scoring (leak.py's masking, applied to all 907 rows of the unfiltered file). The abundance
channel has to be shut for this to be a question about specificity.

Pipeline copied from fulldose.py (dc.pp.extract / prune tmin=5 / adjmat / dc.mt.ulm.func on the
unfiltered file), rank convention copied from leak.py / propagate.py (type_p * score, ranked
descending, divided by the number of scored TFs). leak.py is imported so its 279-experiment
masked table prints above these, and the final assert checks this script's masked rank
fractions reproduce it exactly on the logFC < -1 subset.

PREDICTIONS (the main session's, stated blind before running):
  P1  paired difference (own vs study-mate top-10%) +0.05 to +0.10, 95% CI EXCLUDING 0 once
      clustered on STUDY rather than on TF
  P2  dropping study-mate pairs with regulon Jaccard > 0.05 INCREASES the difference
  P3  the confound is carried mostly by CELL TYPE (Biosample.Name), not by study identity:
      a TF's rank in a different study on the same cell line should also be elevated
  P4  study-controlled baseline 0.20-0.25 (above the 0.123 per-TF base rate, below 0.344)
  Devastating outcome: paired difference CI containing 0 AND study-controlled baseline <= 0.15,
  which would mean most of the masked baseline is context, not TF identity.

REFERENCE LINE is the per-TF base rate, 0.123 (V-022), never 0.10: the rate at which the same
TF lands in the top 10% of experiments where a DIFFERENT TF was perturbed, recomputed here.

CONTROLS. Masked-gene: built in everywhere, there is no unmasked variant in this script. A
degree-preserving shuffled network is owed only if some variant comes out ABOVE the 0.344
masked baseline; the script prints whether any did, and none is expected, since every variant
here is a restriction or a penalisation of that same baseline.

JUDGMENT CALLS (one line each):
  (a) the study is Profile.ID (the GEO series). Pubmed.ID has 107 '-' and 52 'N' placeholders
      and one PMID covering 228 experiments; Data.Source has two values ('GEO'/other) and
      identifies nothing. Profile.ID reproduces V-022's 45-study / 155-experiment count.
  (b) the experiment being judged must be SCORED (its perturbed TF is in CollecTRI at tmin=5),
      but a study-mate need not be: reading TF t's rank in a mate experiment works whatever was
      perturbed there, and that mate's own gene is masked in its own row too. Requiring both to
      be scored costs 16 experiments and 14 studies (139/31 instead of 155/45) and loses the
      cleanest mates, the ones whose perturbed TF is not in CollecTRI at all.
  (c) the primary set is the unfiltered file (643 scored), where the lead lives and where there
      is power; the filtered 279 that carries the 0.344 baseline is reported as a subset row.
  (d) each experiment contributes one mate rate (mean over its mates), so experiments are
      weighted equally regardless of how many mates they have.
  (e) regulon = the nonzero rows of that TF's column in the pruned adjacency, self-edge
      included; Jaccard on those target sets.
  (f) study-controlled hit = top-10% in its own experiment AND strictly better rank there than
      the same TF's rank across its study-mates' experiments; reported against the MEDIAN mate
      (lenient) and against the BEST mate (strict, "beat every study-mate").
  (g) the degree-preserving shuffled network is not owed (no variant proposes a new scoring
      rule) but is run anyway on the paired difference, 5 seeds, since it is cheap and it is
      the only thing that says whether the own-vs-mate gap is regulon-specific.

====================================================================================================
T-053 / T-054 / F-028 EXTENSION (2026-09-18). Blocks (6), (7), (8) below are APPENDED; everything
above them is byte-identical to the F-026 run and still prints under the ASYMMETRIC mask, so
F-026's numbers do not move. V-026 found, in scratch, the one comparison in this area whose CI
excludes 0 and which F-026 does not contain: on the 106 scored experiments that have BOTH a
study-mate (different TF, same Profile.ID) and a same-cell-line mate from a DIFFERENT Profile.ID,
the gap between the two mate sets is the STUDY-CONTEXT LIFT. Both arms share the judged row, which
is what makes it well identified; F-026 part (4)'s +0.116 same-cell-different-study row is not
(V-026: matching on cell line LOWERS d, 0.141 -> 0.116, so that row is the ordinary base-rate
contrast with cell-line composition held fixed).

  (6) T-053: the 106-row contrast, own / within-study mate / same-cell-other-study mate, both
      clusterings, with every n.
  (7) T-054: the SYMMETRIC mask (every perturbed TF's gene zeroed in EVERY row, not only in the
      row where it was perturbed) becomes the reporting convention. One extra ULM run. Every
      touched block is printed under both conventions, side by side, labelled.
  (8) The control V-026 did not run: fulldose.py's degree-preserving shuffler, 10 seeds, on the
      106-row lift specifically. If the lift survives a shuffle it is not regulon-specific and is
      a property of the expression matrix (series/batch structure) rather than of the network,
      which is a different and stronger claim than F-026's.

PREDICTIONS for the extension (stated before running, as T-053 requires):
  P5  all six V-026 numbers reproduce exactly: own 0.245, same-cell-other-study mate 0.163
      (d +0.082, study-CI [+0.027, +0.188]), within-study mate 0.236 (d +0.009 [-0.072, +0.071]),
      lift +0.073 with study-CI [+0.001, +0.192] and TF-CI [+0.018, +0.128], on 106 rows /
      24 studies.
  P6  the symmetric mask moves the LIFT by less than 0.02 (it moves the core d by +0.010,
      +0.056 -> +0.066 [+0.003, +0.134], and leaves the filtered-279 baseline at 0.344/0.348).
  P7  the SHUFFLED lift SHRINKS but does NOT vanish: +0.03 to +0.05, i.e. most of it survives,
      because series/batch structure lives in the expression matrix, not in the network. This is
      the opposite of what the shuffle did to the own-vs-mate gap (+0.056 -> +0.003).
  What would kill the result: the lift's CI containing 0 under the symmetric mask, or the
  shuffled lift MATCHING the real one (which would mean the contrast measures nothing
  network-specific at all and "the series raises a TF's score" needs another mechanism).

REFERENCE LINE below is the per-TF base rate, 0.122, never 0.10 (T-038).
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/studymate.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd, scipy.sparse as sps
import leak                                      # reruns propagate.py + controls.py + the masked 279 table

NB = 2000
print(f"\n{'='*96}\nT-037 / F-026  does masked ULM detect the perturbed TF or the study?\n{'='*96}")

# ------------------------------------------------------------------ load, exactly as fulldose.py
adata = ad.read_h5ad("../data/knocktf_full.h5ad")
net = pd.read_parquet("../data/collectri.parquet")
mat, obs, var = dc.pp.extract(adata)
if sps.issparse(mat): mat = mat.toarray()
pnet = dc.pp.prune(features=var, net=net, tmin=5)
sources, feats, adjm = dc.pp.adjmat(features=var, net=pnet)
src_i = {s: i for i, s in enumerate(sources)}
var_i = {v: i for i, v in enumerate(var)}
obs_i = {o: i for i, o in enumerate(obs)}
nS = len(sources)

tf_of = adata.obs.source.astype(str).values                    # perturbed TF per row
mat_m = mat.copy()                                             # own-gene masking, leak.py's move
for i, t in enumerate(tf_of):
    if t in var_i: mat_m[i, var_i[t]] = 0
S = pd.DataFrame(dc.mt.ulm.func(mat_m, adjm)[0], index=obs, columns=sources)
signed = S.values * adata.obs.type_p.values[:, None]           # type_p = -1 everywhere
R = pd.DataFrame(signed, index=obs, columns=sources).rank(axis=1, ascending=False).values / nS

scored = np.array([t in src_i for t in tf_of])
own_col = np.array([src_i.get(t, 0) for t in tf_of])
fr_own = R[np.arange(len(obs)), own_col]                       # rank fraction of the perturbed TF
logfc = adata.obs.logFC.values.astype(float)
filt = scored & (logfc < -1)                                   # the 279 that carry the 0.344 baseline

# per-TF base rate (V-022): how often TF t is top-10% where a DIFFERENT TF was perturbed
base_rate = {}
for t in sorted(set(tf_of[scored])):
    other = tf_of != t
    base_rate[t] = float(np.mean(R[other, src_i[t]] <= .10)), float(np.mean(R[other, src_i[t]] <= .05))
BASE10 = float(np.mean([base_rate[t][0] for t in tf_of[filt]]))
BASE5 = float(np.mean([base_rate[t][1] for t in tf_of[filt]]))
print(f"\nmasked ULM, {len(obs)} experiments, {scored.sum()} scored over {len(set(tf_of[scored]))} TFs, "
      f"{nS} scored TFs; filtered (logFC < -1) subset {filt.sum()}")
print(f"  masked top-10% on the filtered subset {np.mean(fr_own[filt] <= .10):.3f} "
      f"({int((fr_own[filt] <= .10).sum())}/{int(filt.sum())}), top-5% {np.mean(fr_own[filt] <= .05):.3f}"
      f"   [F-022/V-022: 0.344, 96/279, 0.280]")
print(f"  per-TF BASE RATE on the same rows {BASE10:.3f} top-10%, {BASE5:.3f} top-5% "
      f"[V-022: 0.123] -- this is the reference line everywhere below, not 0.10")
print(f"  masked top-10% on all {scored.sum()} scored {np.mean(fr_own[scored] <= .10):.3f} "
      f"[V-022: 0.299], base rate {np.mean([base_rate[t][0] for t in tf_of[scored]]):.3f}")

# ------------------------------------------------------------------ (1) the study grouping
print(f"\n{'-'*96}\n(1) STUDY GROUPING. candidate obs columns")
for col in ["Profile.ID", "Pubmed.ID", "Data.Source"]:
    v = adata.obs[col].astype(str)
    g = v[scored].groupby(v[scored], observed=True).size()
    ntf = adata.obs[scored].groupby(v[scored], observed=True).source.nunique()
    print(f"  {col:12s} {v.nunique():4d} values ({v[scored].nunique()} among scored), "
          f"{(ntf >= 2).sum():3d} with >= 2 distinct perturbed TFs covering {g[ntf[ntf>=2].index].sum():4d} "
          f"scored experiments; placeholders {(v.isin(['-','N','NA','']).sum())}")
study = adata.obs["Profile.ID"].astype(str).values
print("  -> Profile.ID (GEO series) is the study identifier; Pubmed.ID pools unrelated series "
      "behind placeholder ids, Data.Source is the repository name.")

sdf = pd.DataFrame(dict(study=study, tf=tf_of))
per_study = sdf.groupby("study").tf.nunique()
print(f"  TFs per study over all {len(obs)} experiments: "
      + ", ".join(f"{k} TF(s): {v} studies" for k, v in per_study.value_counts().sort_index().items()))
multi_studies = set(per_study[per_study >= 2].index)
eligible = scored & np.isin(study, list(multi_studies))
print(f"  multi-TF studies {len(multi_studies)}, holding {eligible.sum()} scored experiments "
      f"({eligible.sum()/scored.sum():.0%} of scored) over {len(set(tf_of[eligible]))} TFs "
      f"[V-022 lead: 155 experiments, 101 TFs, 45 studies]")
print(f"  of those, {int((eligible & filt).sum())} are inside the filtered 279 [V-022: 77]")

# ------------------------------------------------------------------ mate machinery
targets = {s: frozenset(np.flatnonzero(adjm[:, c])) for s, c in src_i.items()}
def jac(a, b):
    if a not in targets or b not in targets: return 0.0    # mate TF absent from CollecTRI: no overlap
    A, B = targets[a], targets[b]
    return len(A & B) / len(A | B) if (A or B) else 0.0

def mate_pairs(group, mask, jmax=None):
    """{judged row i -> [mate rows j]} : j shares group[i] and perturbs a different TF.
    i must satisfy mask (scored, so it has a rank of its own); j may be any of the 907 rows."""
    by = {}
    for j in range(len(obs)): by.setdefault(group[j], []).append(j)
    out = {}
    for i in np.flatnonzero(mask):
        js = [j for j in by[group[i]] if tf_of[j] != tf_of[i]
              and (jmax is None or jac(tf_of[i], tf_of[j]) <= jmax)]
        if js: out[i] = js
    return out

def stats(pairs, thr=.10, Rm=None):
    """own hit, mate hit-rate, median and best mate rank per experiment (equal weight per exp)."""
    Rm = R if Rm is None else Rm
    rows = sorted(pairs)
    own = np.array([Rm[i, src_i[tf_of[i]]] <= thr for i in rows], float)
    vals = [Rm[pairs[i], src_i[tf_of[i]]] for i in rows]
    mate = np.array([np.mean(v <= thr) for v in vals])
    med = np.array([np.median(v) for v in vals])
    best = np.array([v.min() for v in vals])
    return np.array(rows), own, mate, med, best

def boot(rows, vals, clusters, seed=0, nb=NB):
    """cluster bootstrap of mean(vals): resample cluster labels with replacement."""
    lab = clusters[rows]
    groups = {}
    for k, l in enumerate(lab): groups.setdefault(l, []).append(k)
    keys = list(groups)
    r = np.random.default_rng(seed)
    out = np.empty(nb)
    for b in range(nb):
        pick = np.concatenate([groups[keys[k]] for k in r.integers(0, len(keys), len(keys))])
        out[b] = vals[pick].mean()
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

# ------------------------------------------------------------------ (2) the core comparison
print(f"\n{'-'*96}\n(2) CORE COMPARISON, paired within study. masked ULM throughout.")
print(f"    (a) rank of TF t in its own experiment   vs   (b) rank of t in the same study's "
      f"experiments on a different TF")
print(f"{'set':26s} {'n_exp':>6s} {'studies':>8s} {'pairs':>6s} {'own10':>7s} {'mate10':>7s} "
      f"{'own5':>6s} {'mate5':>6s} {'ownMed':>7s} {'mateMed':>8s} {'d10':>7s} {'95% study-CI':>18s} "
      f"{'95% TF-CI':>18s}")
core = {}
for label, mask in [("all scored (643)", scored), ("filtered logFC<-1", scored & filt)]:
    P = mate_pairs(study, mask)
    rows, own, mate, med, best = stats(P)
    _, own5, mate5, _, _ = stats(P, .05)
    d = own - mate
    lo_s, hi_s = boot(rows, d, study)
    lo_t, hi_t = boot(rows, d, tf_of)
    allmate = np.concatenate([[R[j, src_i[tf_of[i]]] for j in P[i]] for i in rows])
    core[label] = dict(rows=rows, own=own, mate=mate, med=med, d=d, P=P,
                       ci_s=(lo_s, hi_s), ci_t=(lo_t, hi_t))
    print(f"{label:26s} {len(rows):6d} {len(set(study[rows])):8d} {sum(len(v) for v in P.values()):6d} "
          f"{own.mean():7.3f} {mate.mean():7.3f} {own5.mean():6.3f} {mate5.mean():6.3f} "
          f"{np.median(fr_own[rows]):7.3f} {np.median(allmate):8.3f} {d.mean():+7.3f} "
          f"{f'[{lo_s:+.3f}, {hi_s:+.3f}]':>18s} {f'[{lo_t:+.3f}, {hi_t:+.3f}]':>18s}")
print(f"  reference: per-TF base rate {BASE10:.3f}; the same TFs' top-10% rate in OTHER studies "
      + ", ".join(f"{lab}: {np.mean([np.mean(R[(study != study[i]) & (tf_of != tf_of[i]), src_i[tf_of[i]]] <= .10) for i in core[lab]['rows']]):.3f}"
                  for lab in core))

# ------------------------------------------------------------------ (3) regulon-overlap control
print(f"\n{'-'*96}\n(3) CONTROL V-022 ASKED FOR: drop study-mate pairs with regulon Jaccard > 0.05")
rng = np.random.default_rng(0)
pool = sorted(set(tf_of[scored]))
rand_j = np.mean([jac(a, b) for a, b in zip(rng.choice(pool, 5000), rng.choice(pool, 5000)) if a != b])
P_all = core["all scored (643)"]["P"]
js = np.array([jac(tf_of[i], tf_of[j]) for i in P_all for j in P_all[i]])
tfpairs = sorted({tuple(sorted((tf_of[i], tf_of[j]))) for i in P_all for j in P_all[i]})
jt = np.array([jac(a, b) for a, b in tfpairs])
print(f"  mean regulon Jaccard: study-mate EXPERIMENT pairs {js.mean():.3f} (n={len(js)}), "
      f"study-mate TF pairs {jt.mean():.3f} (n={len(jt)}, V-022's 0.066 is this one), "
      f"random TF pairs {rand_j:.3f} [V-022: 0.024]")
print(f"  {int((js > .05).sum())} of {len(js)} experiment pairs ({np.mean(js > .05):.0%}) and "
      f"{int((jt > .05).sum())} of {len(jt)} TF pairs ({np.mean(jt > .05):.0%}) exceed Jaccard 0.05")
big = max(set(study[np.isin(study, list(multi_studies))]), key=lambda s: (study == s).sum())
small = {tuple(sorted((tf_of[i], tf_of[j]))) for i in P_all for j in P_all[i] if study[i] != big}
print(f"  by study size: TF pairs from the largest study ({big}, {(study == big).sum()} experiments) "
      f"{np.mean([jac(*p) for p in tfpairs if p not in small]):.3f}, TF pairs from every other "
      f"multi-TF study {np.mean([jac(*p) for p in small]):.3f} (n={len(small)}) -- "
      f"the enrichment V-022 reported lives in the small deliberate studies, not in the big screen")
print(f"{'set':38s} {'n_exp':>6s} {'studies':>8s} {'pairs':>6s} {'own10':>7s} {'mate10':>7s} "
      f"{'d10':>7s} {'95% study-CI':>18s} {'95% TF-CI':>18s}")
for label, mask in [("all scored (643)", scored), ("filtered logFC<-1", scored & filt)]:
    keep = None
    for jmax, tag in [(None, "all mates"), (.05, "Jaccard <= 0.05"), (None, "same rows, all mates")]:
        P = mate_pairs(study, mask, jmax)
        if tag == "Jaccard <= 0.05": keep = set(P)              # rows that survive the filter
        if tag.startswith("same rows"): P = {i: v for i, v in P.items() if i in keep}
        rows, own, mate, med, best = stats(P)
        d = own - mate
        lo_s, hi_s = boot(rows, d, study); lo_t, hi_t = boot(rows, d, tf_of)
        print(f"{label + ', ' + tag:38s} {len(rows):6d} {len(set(study[rows])):8d} "
              f"{sum(len(v) for v in P.values()):6d} {own.mean():7.3f} {mate.mean():7.3f} "
              f"{d.mean():+7.3f} {f'[{lo_s:+.3f}, {hi_s:+.3f}]':>18s} {f'[{lo_t:+.3f}, {hi_t:+.3f}]':>18s}")
print("  the third row of each block holds the experiments fixed and puts the high-overlap mates "
      "back,\n  so it separates which experiments survive the filter from what the filter does to "
      "the mate rate")

# ------------------------------------------------------------------ (4) what does "study" carry?
print(f"\n{'-'*96}\n(4) WHAT DOES 'STUDY' CARRY? same comparison, grouping by each candidate.")
print(f"    the key contrast is the last row: same cell line, DIFFERENT study.")
print(f"{'grouping':34s} {'n_exp':>6s} {'groups':>7s} {'pairs':>7s} {'own10':>7s} {'mate10':>7s} "
      f"{'d10':>7s} {'95% group-CI':>18s}")
cell = adata.obs["Biosample.Name"].astype(str).values
method = adata.obs["Knock.Method"].astype(str).values
plat = adata.obs["Platform"].astype(str).values
groupings = [("study (Profile.ID)", study), ("cell type (Biosample.Name)", cell),
             ("knockdown (Knock.Method)", method), ("platform (Platform)", plat)]
for name, g in groupings:
    P = mate_pairs(g, scored)
    rows, own, mate, med, best = stats(P)
    d = own - mate
    lo, hi = boot(rows, d, g)
    print(f"{name:34s} {len(rows):6d} {len(set(g[rows])):7d} {sum(len(v) for v in P.values()):7d} "
          f"{own.mean():7.3f} {mate.mean():7.3f} {d.mean():+7.3f} {f'[{lo:+.3f}, {hi:+.3f}]':>18s}")
# same cell line, different study: mates restricted to a different Profile.ID
P = mate_pairs(cell, scored)
P_x = {i: [j for j in v if study[j] != study[i]] for i, v in P.items()}
P_x = {i: v for i, v in P_x.items() if v}
rows, own, mate, med, best = stats(P_x)
d = own - mate
lo, hi = boot(rows, d, cell)
print(f"{'same cell line, DIFFERENT study':34s} {len(rows):6d} {len(set(cell[rows])):7d} "
      f"{sum(len(v) for v in P_x.values()):7d} {own.mean():7.3f} {mate.mean():7.3f} {d.mean():+7.3f} "
      f"{f'[{lo:+.3f}, {hi:+.3f}]':>18s}")
# and the complement: same study, same cell line (does study add anything beyond cell line?)
P_ss = {i: [j for j in v if study[j] == study[i]] for i, v in mate_pairs(cell, scored).items()}
P_ss = {i: v for i, v in P_ss.items() if v}
rows, own, mate, med, best = stats(P_ss)
d = own - mate
lo, hi = boot(rows, d, study)
print(f"{'same cell line, SAME study':34s} {len(rows):6d} {len(set(cell[rows])):7d} "
      f"{sum(len(v) for v in P_ss.values()):7d} {own.mean():7.3f} {mate.mean():7.3f} {d.mean():+7.3f} "
      f"{f'[{lo:+.3f}, {hi:+.3f}]':>18s}")

# ------------------------------------------------------------------ (5) study-controlled baseline
print(f"\n{'-'*96}\n(5) STUDY-CONTROLLED BASELINE: hit = top-10% in its own experiment AND better "
      f"rank there\n    than the same TF's rank in its study-mates' experiments (vs the MEDIAN "
      f"mate, and vs the BEST mate).")
print(f"{'set':34s} {'n_exp':>6s} {'naive10':>8s} {'vs median':>10s} {'95% study-CI':>18s} "
      f"{'vs best':>8s} {'95% study-CI':>18s} {'naive5':>7s} {'vsMed5':>8s} {'baseline':>9s}")
sc_out = {}
for label, mask in [("multi-TF studies, all scored", scored), ("multi-TF studies, filtered 279", scored & filt)]:
    P = mate_pairs(study, mask)
    rows, own, mate, med, bst = stats(P)
    _, own5, _, _, _ = stats(P, .05)
    sc = own * (fr_own[rows] < med)
    scb = own * (fr_own[rows] < bst)
    sc5 = own5 * (fr_own[rows] < med)
    lo, hi = boot(rows, sc, study)
    lob, hib = boot(rows, scb, study)
    b = float(np.mean([base_rate[t][0] for t in tf_of[rows]]))
    sc_out[label] = (own.mean(), sc.mean(), lo, hi, scb.mean(), lob, hib, b)
    print(f"{label:34s} {len(rows):6d} {own.mean():8.3f} {sc.mean():10.3f} "
          f"{f'[{lo:.3f}, {hi:.3f}]':>18s} {scb.mean():8.3f} {f'[{lob:.3f}, {hib:.3f}]':>18s} "
          f"{own5.mean():7.3f} {sc5.mean():8.3f} {b:9.3f}")
print(f"  for comparison the whole filtered set: masked ULM {np.mean(fr_own[filt] <= .10):.3f} "
      f"(0.344), per-TF base rate {BASE10:.3f} (0.123). Only "
      f"{(eligible & filt).sum()}/{filt.sum()} = {(eligible & filt).sum()/filt.sum():.0%} of the "
      f"filtered experiments can be study-controlled at all.")

# ------------------------------------------------------------------ CONTROL: shuffled network
def shuffle_tf_gene(adj, seed):
    """copied from fulldose.py: each TF keeps its number of +1 and of -1 targets, drawn at random
    from the measured genes other than its own, so regulon size and sign mix survive and identity
    does not."""
    r = np.random.default_rng(seed)
    B = np.zeros_like(adj)
    for s, c in src_i.items():
        w = adj[:, c][adj[:, c] != 0]
        if len(w) == 0: continue
        pool = np.arange(adj.shape[0])
        if s in var_i: pool = pool[pool != var_i[s]]
        B[r.choice(pool, size=len(w), replace=False), c] = w
    return B

hi_variant = max([core[l]["own"].mean() for l in core] + [v[1] for v in sc_out.values()])
print(f"\n{'-'*96}\n(CONTROL) degree-preserving shuffled TF-gene network, 5 seeds, masked matrix, "
      f"same study pairing.\n  Not owed -- no variant here proposes a new scoring rule, they all "
      f"re-count the same masked ULM scores;\n  the one number above the 0.344 baseline "
      f"({hi_variant:.3f}) is a row subset of it, not a method. Run anyway.")
P = mate_pairs(study, scored)
print(f"{'seed':>6s} {'own10':>7s} {'mate10':>7s} {'d10':>8s} {'own5':>7s} "
      f"{'better':>7s} {'worse':>7s}   (better/worse = shuffled vs real masked rank, 643 scored)")
sh_d = []
for seed in range(5):
    S_sh = dc.mt.ulm.func(mat_m, shuffle_tf_gene(adjm, seed))[0]
    R_sh = pd.DataFrame(S_sh * adata.obs.type_p.values[:, None]).rank(axis=1, ascending=False).values / nS
    rows, own_sh, mate_sh, _, _ = stats(P, .10, R_sh)
    _, own_sh5, _, _, _ = stats(P, .05, R_sh)
    fr_sh = R_sh[np.arange(len(obs)), own_col][scored]
    sh_d.append(own_sh.mean() - mate_sh.mean())
    print(f"{seed:6d} {own_sh.mean():7.3f} {mate_sh.mean():7.3f} {sh_d[-1]:+8.3f} {own_sh5.mean():7.3f} "
          f"{np.mean(fr_sh < fr_own[scored]):7.2f} {np.mean(fr_sh > fr_own[scored]):7.2f}")
print(f"  real +{core['all scored (643)']['d'].mean():.3f}   shuffled mean {np.mean(sh_d):+.3f} "
      f"[{min(sh_d):+.3f}, {max(sh_d):+.3f}]")

# ==================================================================================================
# T-053 / T-054 / F-028.  Everything above this line is byte-identical to the F-026 run.
# ==================================================================================================
print(f"\n{'='*96}\n(6) T-053  THE STUDY-CONTEXT LIFT, on IDENTICAL rows.  [ASYMMETRIC mask, "
      f"F-026's convention]\n{'='*96}")
Pw_all = mate_pairs(study, scored)                             # within-study mate, different TF
Pc_all = mate_pairs(cell, scored)
Px_all = {i: [j for j in v if study[j] != study[i]] for i, v in Pc_all.items()}
Px_all = {i: v for i, v in Px_all.items() if v}                # same cell line, DIFFERENT study
rows106 = np.array(sorted(set(Pw_all) & set(Px_all)))          # rows where BOTH contrasts exist
Pw = {i: Pw_all[i] for i in rows106}
Px = {i: Px_all[i] for i in rows106}
BASE106 = float(np.mean([base_rate[t][0] for t in tf_of[rows106]]))
BASE106_5 = float(np.mean([base_rate[t][1] for t in tf_of[rows106]]))
print(f"  set: {len(rows106)} scored experiments in {len(set(study[rows106]))} studies over "
      f"{len(set(tf_of[rows106]))} TFs and {len(set(cell[rows106]))} cell lines; "
      f"{sum(len(v) for v in Pw.values())} within-study mate pairs, "
      f"{sum(len(v) for v in Px.values())} same-cell-different-study mate pairs; "
      f"{int((filt[rows106]).sum())} of them inside the filtered 279")
print(f"  reference line: per-TF base rate on these rows {BASE106:.3f} top-10% / {BASE106_5:.3f} "
      f"top-5% (global 0.122, T-038) -- never 0.10")

def lift(Rm, thr=.10):
    """on rows106: own hit, within-study mate rate, same-cell-other-study mate rate (one per row)."""
    own = np.array([Rm[i, src_i[tf_of[i]]] <= thr for i in rows106], float)
    mw = np.array([np.mean(Rm[Pw[i], src_i[tf_of[i]]] <= thr) for i in rows106])
    mx = np.array([np.mean(Rm[Px[i], src_i[tf_of[i]]] <= thr) for i in rows106])
    return own, mw, mx

def contrast(label, a, b, av, bv, ref=""):
    d = av - bv
    lo_s, hi_s = boot(rows106, d, study)
    lo_t, hi_t = boot(rows106, d, tf_of)
    print(f"{label:46s} {a:7.3f} {b:7.3f} {d.mean():+7.3f} {f'[{lo_s:+.3f}, {hi_s:+.3f}]':>18s} "
          f"{f'[{lo_t:+.3f}, {hi_t:+.3f}]':>18s}  {ref}")
    return d.mean(), (lo_s, hi_s), (lo_t, hi_t)

L = {}
for thr, tag in [(.10, "top-10%"), (.05, "top-5%")]:
    own, mw, mx = lift(R, thr)
    print(f"\n  {tag}, masked ULM (asymmetric), one rate per experiment, 106 rows")
    print(f"{'contrast (a vs b)':46s} {'a':>7s} {'b':>7s} {'d':>7s} {'95% study-CI':>18s} "
          f"{'95% TF-CI':>18s}")
    L[(tag, "own vs cellx")] = contrast("own  vs  same-cell OTHER-study mate", own.mean(), mx.mean(),
                                        own, mx, "V-026: 0.245 / 0.163 / +0.082 [+0.027, +0.188]")
    L[(tag, "own vs within")] = contrast("own  vs  WITHIN-study mate", own.mean(), mw.mean(),
                                         own, mw, "V-026: 0.245 / 0.236 / +0.009 [-0.072, +0.071]")
    L[(tag, "lift")] = contrast("LIFT: within-study mate  vs  other-study mate", mw.mean(), mx.mean(),
                                mw, mx, "V-026: +0.073 [+0.001, +0.192] / [+0.018, +0.128]")
own10, mw10, mx10 = lift(R, .10)
print(f"\n  reading: both arms share the judged row, so the lift is the study-context effect with "
      f"TF identity,\n  cell line and row composition held fixed. Being in the same GEO series "
      f"raises a TF's top-10% rate\n  from {mx10.mean():.3f} to {mw10.mean():.3f}; being the TF that "
      f"was actually perturbed raises it to {own10.mean():.3f}.")

# ------------------------------------------------------------------ (7) T-054 symmetric mask
print(f"\n{'='*96}\n(7) T-054  SYMMETRIC MASK as the reporting convention: every perturbed TF's gene\n"
      f"    zeroed in EVERY row, not only in the row where it was perturbed. One extra ULM run.\n"
      f"    F-026's blocks above keep the ASYMMETRIC convention so their numbers do not move;\n"
      f"    every block this touches is printed under both, labelled.\n{'='*96}")
mat_sym = mat.copy()
pert_cols = sorted({var_i[t] for t in set(tf_of) if t in var_i})
mat_sym[:, pert_cols] = 0
S_sym = pd.DataFrame(dc.mt.ulm.func(mat_sym, adjm)[0], index=obs, columns=sources)
R_sym = pd.DataFrame(S_sym.values * adata.obs.type_p.values[:, None],
                     index=obs, columns=sources).rank(axis=1, ascending=False).values / nS
fr_sym = R_sym[np.arange(len(obs)), own_col]
base_sym = {}
for t in sorted(set(tf_of[scored])):
    other = tf_of != t
    base_sym[t] = float(np.mean(R_sym[other, src_i[t]] <= .10)), float(np.mean(R_sym[other, src_i[t]] <= .05))
print(f"  {len(pert_cols)} perturbed-TF genes zeroed in all {len(obs)} rows "
      f"({int((np.abs(mat_m) > 0).sum() - (np.abs(mat_sym) > 0).sum())} extra nonzero cells removed "
      f"vs the asymmetric mask); masked/unmasked own-gene check is in the assert")
print(f"\n  BASELINE (the 0.344 number) under both conventions")
print(f"{'convention':22s} {'filt10':>7s} {'filt5':>7s} {'all10':>7s} {'base10':>7s} {'base5':>7s} "
      f"{'better':>7s} {'worse':>7s}   (better/worse = vs asymmetric masked ULM, 643 scored)")
print(f"{'asymmetric (F-026)':22s} {np.mean(fr_own[filt] <= .10):7.3f} "
      f"{np.mean(fr_own[filt] <= .05):7.3f} {np.mean(fr_own[scored] <= .10):7.3f} "
      f"{BASE10:7.3f} {BASE5:7.3f} {'--':>7s} {'--':>7s}")
BASE10s = float(np.mean([base_sym[t][0] for t in tf_of[filt]]))
BASE5s = float(np.mean([base_sym[t][1] for t in tf_of[filt]]))
print(f"{'symmetric (T-054)':22s} {np.mean(fr_sym[filt] <= .10):7.3f} "
      f"{np.mean(fr_sym[filt] <= .05):7.3f} {np.mean(fr_sym[scored] <= .10):7.3f} "
      f"{BASE10s:7.3f} {BASE5s:7.3f} {np.mean(fr_sym[scored] < fr_own[scored]):7.2f} "
      f"{np.mean(fr_sym[scored] > fr_own[scored]):7.2f}   [V-026: 0.348 vs 0.344]")

def core_row(label, Rm, mask):
    P = mate_pairs(study, mask)
    rows, own, mate, med, best = stats(P, .10, Rm)
    _, own5, mate5, _, _ = stats(P, .05, Rm)
    d = own - mate
    lo_s, hi_s = boot(rows, d, study)
    lo_t, hi_t = boot(rows, d, tf_of)
    print(f"{label:40s} {len(rows):6d} {own.mean():7.3f} {mate.mean():7.3f} {own5.mean():6.3f} "
          f"{mate5.mean():6.3f} {d.mean():+7.3f} {f'[{lo_s:+.3f}, {hi_s:+.3f}]':>18s} "
          f"{f'[{lo_t:+.3f}, {hi_t:+.3f}]':>18s}")
    return d.mean()

print(f"\n  BLOCK (2) CORE COMPARISON under both conventions")
print(f"{'set / convention':40s} {'n_exp':>6s} {'own10':>7s} {'mate10':>7s} {'own5':>6s} "
      f"{'mate5':>6s} {'d10':>7s} {'95% study-CI':>18s} {'95% TF-CI':>18s}")
core_row("all scored (643), asymmetric", R, scored)
d_sym_core = core_row("all scored (643), SYMMETRIC", R_sym, scored)
core_row("filtered logFC<-1, asymmetric", R, scored & filt)
core_row("filtered logFC<-1, SYMMETRIC", R_sym, scored & filt)
print(f"  [V-026: the asymmetry biases d DOWNWARD; symmetric core d +0.066 [+0.003, +0.134]]")

def cell_row(label, Rm, P):
    rows, own, mate, med, best = stats(P, .10, Rm)
    d = own - mate
    lo, hi = boot(rows, d, cell)
    print(f"{label:40s} {len(rows):6d} {sum(len(v) for v in P.values()):7d} {own.mean():7.3f} "
          f"{mate.mean():7.3f} {d.mean():+7.3f} {f'[{lo:+.3f}, {hi:+.3f}]':>18s}")

print(f"\n  BLOCK (4) the two cell-line rows that feed (6), under both conventions")
print(f"{'grouping / convention':40s} {'n_exp':>6s} {'pairs':>7s} {'own10':>7s} {'mate10':>7s} "
      f"{'d10':>7s} {'95% cell-CI':>18s}")
cell_row("same cell, DIFFERENT study, asymmetric", R, Px_all)
cell_row("same cell, DIFFERENT study, SYMMETRIC", R_sym, Px_all)
P_ss2 = {i: [j for j in v if study[j] == study[i]] for i, v in Pc_all.items()}
P_ss2 = {i: v for i, v in P_ss2.items() if v}
cell_row("same cell, SAME study, asymmetric", R, P_ss2)
cell_row("same cell, SAME study, SYMMETRIC", R_sym, P_ss2)
print(f"  [V-026: matching on cell line LOWERS d (0.141 -> 0.116), so the +0.116 row is the "
      f"ordinary base-rate\n   contrast with cell-line composition held fixed, NOT a study-context "
      f"estimate. (6) is the identified one.]")

print(f"\n  BLOCK (5) STUDY-CONTROLLED BASELINE under both conventions")
print(f"{'set / convention':40s} {'n_exp':>6s} {'naive10':>8s} {'vs median':>10s} "
      f"{'95% study-CI':>18s} {'vs best':>8s} {'95% study-CI':>18s} {'baseline':>9s}")
for label, mask in [("multi-TF, all scored", scored), ("multi-TF, filtered 279", scored & filt)]:
    for tag, Rm, fr_, br in [("asymmetric", R, fr_own, base_rate), ("SYMMETRIC", R_sym, fr_sym, base_sym)]:
        Ps = mate_pairs(study, mask)
        rows, own, mate, med, bst = stats(Ps, .10, Rm)
        sc = own * (fr_[rows] < med)
        scb = own * (fr_[rows] < bst)
        lo, hi = boot(rows, sc, study)
        lob, hib = boot(rows, scb, study)
        print(f"{label + ', ' + tag:40s} {len(rows):6d} {own.mean():8.3f} {sc.mean():10.3f} "
              f"{f'[{lo:.3f}, {hi:.3f}]':>18s} {scb.mean():8.3f} {f'[{lob:.3f}, {hib:.3f}]':>18s} "
              f"{float(np.mean([br[t][0] for t in tf_of[rows]])):9.3f}")

print(f"\n  BLOCK (6) THE LIFT under both conventions (the T-053 result)")
print(f"{'contrast / convention':46s} {'a':>7s} {'b':>7s} {'d':>7s} {'95% study-CI':>18s} "
      f"{'95% TF-CI':>18s}")
lift_real = {}
for tag, Rm in [("asymmetric", R), ("SYMMETRIC", R_sym)]:
    own_, mw_, mx_ = lift(Rm)
    contrast(f"own vs other-study mate, {tag}", own_.mean(), mx_.mean(), own_, mx_)
    contrast(f"own vs within-study mate, {tag}", own_.mean(), mw_.mean(), own_, mw_)
    lift_real[tag] = contrast(f"LIFT within vs other study, {tag}", mw_.mean(), mx_.mean(), mw_, mx_)
print(f"  symmetric minus asymmetric on the LIFT: "
      f"{lift_real['SYMMETRIC'][0] - lift_real['asymmetric'][0]:+.3f} "
      f"(P6 predicted a move smaller than 0.020)")
print(f"  CAUTION on 'the CI excludes 0': the study-clustered lower bound on the lift sits ON zero "
      f"under the\n  ASYMMETRIC mask ({lift_real['asymmetric'][1][0]:+.4f} here, +0.001 in V-026's "
      f"scratch -- the difference is bootstrap noise between two\n  implementations, not data). It "
      f"excludes 0 under the SYMMETRIC mask "
      f"({lift_real['SYMMETRIC'][1][0]:+.3f}) and under\n  TF clustering in both. Quote the "
      f"symmetric row; do not quote a knife-edge bound as if it were robust.")
print(f"  CONVENTION GOING FORWARD: SYMMETRIC. It is the honest one -- under the asymmetric mask a "
      f"mate row\n  still carries the judged TF's own transcript, so the mate arm keeps an "
      f"abundance channel the own arm\n  has had removed, which biases every own-vs-mate d DOWNWARD.")

# ------------------------------------------------------------------ (8) CONTROL on the 106 rows
print(f"\n{'='*96}\n(8) CONTROL T-053 OWES, which V-026 did not run: degree-preserving shuffled\n"
      f"    TF-gene network (fulldose.py's shuffler), 10 seeds, on the 106-row lift specifically.\n"
      f"    PREDICTED (P7) before running: the shuffled lift SHRINKS but SURVIVES, +0.03 to +0.05,\n"
      f"    because series/batch structure is in the expression matrix, not in the network.\n{'='*96}")
print(f"{'seed':>4s} {'own10':>7s} {'within':>7s} {'other':>7s} {'LIFT':>8s} {'own10':>7s} "
      f"{'within':>7s} {'other':>7s} {'LIFT':>8s} {'better':>7s} {'worse':>7s}")
print(f"{'':>4s} {'--- asymmetric mask ---':^32s} {'--- SYMMETRIC mask ---':^32s} "
      f"{'(sym shuffled vs sym real, 643)':>16s}")
sh_lift = {"asymmetric": [], "SYMMETRIC": []}
for seed in range(10):
    B = shuffle_tf_gene(adjm, seed)
    cells, line = {}, f"{seed:4d}"
    for tag, M in [("asymmetric", mat_m), ("SYMMETRIC", mat_sym)]:
        R_sh = pd.DataFrame(dc.mt.ulm.func(M, B)[0] * adata.obs.type_p.values[:, None]
                            ).rank(axis=1, ascending=False).values / nS
        own_, mw_, mx_ = lift(R_sh)
        sh_lift[tag].append(mw_ - mx_)
        cells[tag] = R_sh
        line += f" {own_.mean():7.3f} {mw_.mean():7.3f} {mx_.mean():7.3f} {(mw_ - mx_).mean():+8.3f}"
    fr_shs = cells["SYMMETRIC"][np.arange(len(obs)), own_col][scored]
    print(line + f" {np.mean(fr_shs < fr_sym[scored]):7.2f} {np.mean(fr_shs > fr_sym[scored]):7.2f}")
print(f"\n{'shuffled lift (mean over 10 seeds)':46s} {'':>7s} {'':>7s} {'d':>7s} "
      f"{'95% study-CI':>18s} {'95% TF-CI':>18s}")
for tag in ("asymmetric", "SYMMETRIC"):
    m = np.mean(sh_lift[tag], axis=0)                      # per-row lift averaged over seeds
    lo_s, hi_s = boot(rows106, m, study)
    lo_t, hi_t = boot(rows106, m, tf_of)
    per_seed = [float(a.mean()) for a in sh_lift[tag]]
    print(f"{'shuffled network, ' + tag:46s} {'':>7s} {'':>7s} {np.mean(per_seed):+7.3f} "
          f"{f'[{lo_s:+.3f}, {hi_s:+.3f}]':>18s} {f'[{lo_t:+.3f}, {hi_t:+.3f}]':>18s}  "
          f"per-seed [{min(per_seed):+.3f}, {max(per_seed):+.3f}]")
    print(f"{'  real (same rows, same pairing), ' + tag:46s} {'':>7s} {'':>7s} "
          f"{lift_real[tag][0]:+7.3f} {f'[{lift_real[tag][1][0]:+.3f}, {lift_real[tag][1][1]:+.3f}]':>18s} "
          f"{f'[{lift_real[tag][2][0]:+.3f}, {lift_real[tag][2][1]:+.3f}]':>18s}")
surv = np.mean([a.mean() for a in sh_lift["SYMMETRIC"]]) / lift_real["SYMMETRIC"][0]
print(f"\n  fraction of the symmetric lift surviving the shuffle: {surv:.0%}, with a CI covering 0 "
      f"and a per-seed\n  range straddling it. P7 IS CONTRADICTED: the lift does NOT survive a "
      f"degree-preserving shuffle, it dies\n  like the own-vs-mate gap did (+0.056 -> +0.003). So "
      f"the study-context effect is NOT generic series/batch\n  variance that any TF-gene map of "
      f"the right degree sequence would read out -- it runs THROUGH the real\n  regulons: inside a "
      f"GEO series a TF's actual targets move together, which is a stronger and more\n  specific "
      f"claim than the one predicted, and it is the same mechanism class the 0.344 baseline "
      f"rests on.")
print(f"\n  ROUTE (V-026's language). This is TF-perturbing (knockTF), masked throughout, and it "
      f"concerns\n  {len(rows106)} of {int(scored.sum())} scored experiments "
      f"({len(rows106)/scored.sum():.0%}), {int(filt[rows106].sum())} of them in the filtered 279. "
      f"It does NOT move the\n  0.344 baseline, which stands; it does not license calling that "
      f"baseline mostly study context, and it\n  says nothing about the {int(scored.sum()) - 155} "
      f"scored experiments with no study-mate at all. What it DOES license: on\n  the rows where "
      f"the test is possible, being in the same GEO series is worth about as much to a TF's\n"
      f"  score as being the TF that was perturbed, so a study-controlled baseline is the honest "
      f"lower bound.")

# ------------------------------------------------------------------ F-029: persist
# One row per SCORED experiment: study, cell line, both mask conventions' rank fractions, the
# per-TF base rate under each, and the mate counts that define the 155- and 106-row subsets.
_s = np.flatnonzero(scored)
pd.DataFrame(dict(exp=np.asarray(obs)[_s], tf=tf_of[_s], study=study[_s], cell=cell[_s],
                  own_logfc=logfc[_s], in_filtered_279=filt[_s],
                  masked_frac=fr_own[_s], masked_sym_frac=fr_sym[_s],
                  base_rate10=[base_rate[t][0] for t in tf_of[_s]],
                  base_rate10_sym=[base_sym[t][0] for t in tf_of[_s]],
                  n_within_study_mates=[len(Pw_all.get(i, [])) for i in _s],
                  n_same_cell_other_study_mates=[len(Px_all.get(i, [])) for i in _s],
                  in_106=np.isin(_s, rows106))).to_csv("../results/studymate.csv", index=False)

# ------------------------------------------------------------------ assert
ref = np.sort(leak.evaluate(leak.S))                          # masked ULM on the filtered file
mine = np.sort(fr_own[filt])          # max |diff| 0.004: leak.py's file drops 50 all-zero genes
rg = [(i, var_i[t]) for i, t in enumerate(tf_of) if t in var_i]
no_self = all(tf_of[j] != tf_of[i] and study[j] == study[i] for i in P for j in P[i])
assert (len(ref) == filt.sum() and np.abs(ref - mine).max() < 0.02
        and (ref <= .10).sum() == (mine <= .10).sum() and no_self
        and all(mat_m[i, g] == 0 for i, g in rg)
        and all(mat[i, g] != 0 for i, g in rg if logfc[i] < -1)
        # T-053: both mate sets exist on every one of the 106 rows, and each is what it claims
        and len(rows106) > 0 and set(rows106) <= set(np.flatnonzero(scored))
        and all(study[j] == study[i] and tf_of[j] != tf_of[i] for i in Pw for j in Pw[i])
        and all(cell[j] == cell[i] and study[j] != study[i] and tf_of[j] != tf_of[i]
                for i in Px for j in Px[i])
        # T-054: the symmetric mask zeroes every perturbed TF's gene in every row, and strictly more
        # cells than the asymmetric one, while leaving non-perturbed genes untouched
        and np.abs(mat_sym[:, pert_cols]).max() == 0
        and (np.abs(mat_sym) > 0).sum() < (np.abs(mat_m) > 0).sum()
        and np.array_equal(np.delete(mat_sym, pert_cols, axis=1), np.delete(mat, pert_cols, axis=1))), \
    "masked ranks do not reproduce leak.py's 279, or a study-mate shares the perturbed TF, or masking " \
    "failed, or the 106-row mate sets are not both-kinds-on-identical-rows, or the symmetric mask is wrong"
