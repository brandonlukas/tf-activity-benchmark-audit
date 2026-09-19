"""F-029 / T-055: the machine-readable record. One row per SCORED knockTF experiment, with every
condition's rank fraction and hit flag plus the covariates, so that most headline numbers in
findings.md are recomputable from results/master.csv alone without rerunning anything.

WHY IT IS BUILT ON THE UNFILTERED FILE. data/knocktf_full.h5ad is `dc.ds.knocktf(thr_fc=None)`:
907 experiments, 643 of them scored (perturbed TF is a CollecTRI source surviving tmin=5) over
263 TFs. The 279/155 set that carries F-001 through F-021 is the logFC < -1 subset, marked by
`in_filtered_279`. V-013 established that the denominator is the single most-misquoted thing in
this log, so both denominators live in one file and every rate is a groupby away.
fulldose.py's docstring is the licence for doing it this way: the var axis of the two files is
NOT identical after dc.pp.extract: it drops 50 all-zero genes from the filtered file, so leak.py/selfedge.py score on 21935 features and this script on 21985 (V-029) and ULM scores rows independently, so the logFC < -1 subset of the
unfiltered pipeline reproduces the filtered pipeline exactly. Propagation (S + a*S@A) is also
row-independent, so the same holds for the propagated conditions.

MACHINERY IS COPIED, NOT REIMPLEMENTED, from the scripts the Index says are citable:
  load / prune / adjmat / dc.mt.ulm.func     selfedge.py, via fulldose.py's unfiltered version
  asymmetric own-gene mask                   leak.py  (F-005/F-013: 0.344, the citable baseline)
  symmetric mask                             studymate.py (7)  (F-028/V-028's convention)
  no-self-edge / no-other-regulon            selfedge.py  (F-010)
  propagation up/sum/alpha=0.5               propagate.py + leak.py  (F-004)
  per-TF base rate                           studymate.py  (V-022: 0.123, never 0.10)
  study grouping = Profile.ID                studymate.py (1)  (Pubmed.ID has 107 '-' placeholders)
  max regulon Jaccard                        gapcheck.py  (F-002's overlap measure)

CONVENTIONS WHERE TWO SCRIPTS DISAGREE, resolved by what the Index calls citable:
  * MASKING. `ulm_masked` is leak.py's ASYMMETRIC mask: the perturbed TF's own gene is zeroed
    only in ITS OWN row. That is the convention behind the citable 0.344 / 96 of 279
    (F-005, V-013) and behind every F-010/F-011/F-016/F-017 number, so it is the default here.
    `ulm_masked_symmetric` is F-028's convention: EVERY perturbed TF's gene zeroed in EVERY row
    (766 columns unchanged, 263 gene columns blanked everywhere). V-028 calls symmetric the
    honest convention going forward for between-experiment comparisons; it is not the
    convention any pre-F-026 headline was computed under, and it is NOT asserted against one.
  * BASE RATE. V-022's construction, as coded in studymate.py: for TF t, the fraction of the
    907 rows where a DIFFERENT TF was perturbed in which t is itself top-10% under the
    ASYMMETRIC masked scores. Reported per row as base_rate10 / base_rate5. Mean over the
    filtered 279 is V-022's 0.123. studymate.py's own docstring says 0.122; the Index says
    0.123, so 0.123 is the target.
  * REGULON SIZE. Nonzero column count of the PRUNED (tmin=5) adjacency, self-edge INCLUDED.
    fulldose.py drops the self-edge first; gapcheck.py counts unpruned net targets. The pruned
    count is the one F-017's size bins and V-027's law are about, so it is the one stored.
  * type_p. -1 in every knockTF row; rank is on `type_p * score`, descending, over all 766
    scored TFs, exactly as propagate.py/leak.py/selfedge.py/studymate.py all do it.

VALIDATION IS THE POINT, NOT THE FILE. Every published number below is recomputed from the
assembled table alone and asserted. PREDICTIONS, stated before running:
  EXACTLY (integer counts and 3-decimal rates), because nothing here is seeded and every step
  is row-independent, so the unfiltered pipeline must land on the filtered pipeline's rows:
    F-001            0.452, 126/279, 155 TFs
    F-005/F-013      masked 0.344, 96/279
    F-010            no-self-edge 96/279, no-other 127/279, 197/82 split,
                     0.538/0.386 (self-edge unmasked/masked), 0.244/0.244 (no self-edge)
    F-011            176 +1 / 21 -1; +1 group 98 -> 65 hits, -1 group 8 -> 11
    F-016/F-017      the 176 split 92 default activation / 74 PMID / 10 regulon
    F-022            643 scored, 263 TFs, masked 0.299
    V-022            base rate 0.123
  TO 3 DECIMALS ONLY:
    F-004            propagated 0.566 unmasked, 0.276 masked. These are the two numbers whose
                     reproduction depends on rebuilding the TF-TF matrix A by a different route
                     (pivot_table over the unpruned net on this file's source list) rather than
                     on re-ranking the same scores, so a tie-break or a duplicate-edge
                     difference could move them by one experiment (0.0036) without anything
                     being wrong. Nothing here uses a random seed, so no number is seed-noisy.
  WHAT A FAILURE MEANS. A headline that does NOT come back is not a bug in this file to be
  tuned away. It means the number in findings.md was computed on a denominator, a mask
  convention or a file other than the one its entry names, and the entry is wrong about its
  own scope. The mismatch is then the finding, and it is reported as such.

OUTCOME (written after the first run; the PREDICTIONS above are left exactly as stated).
  28 of 29 reproduced, all of them exactly, including both F-004 numbers that were only
  predicted to 3 decimals. ONE MISSED, and it is a denominator slip of the kind V-013 warned
  about, not a computation error:
    V-022's per-TF base rate "0.123 on the filtered 279" recomputes to 0.1225, i.e. 0.122.
    0.1231 -> 0.123 IS this construction's value on all 643 SCORED rows. V-022's own table
    subtracts the same 0.123 from both rows ("FILTERED 279 0.344 - 0.123 = +0.221" and "all
    643 0.299 - 0.123 = +0.175"), so the 643 figure was carried onto the 279 row. Two things
    confirm the construction here IS V-022's and not a near-miss: the top-5% companion
    reproduces exactly (0.0666 -> 0.067, V-022's "0.280 - 0.067 = +0.213"), and studymate.py
    already prints 0.122 for this quantity, which V-028 records as "global 0.122". So the log
    contains both values for one number and the Index, the V-022 entry and CLAUDE.md all quote
    the 643 one against the 279 denominator.
    Consequence, and it is small but it runs the right way: masked ULM's margin over its base
    rate on the 279 is +0.222, not +0.221, and the multiple is 2.81x, not 2.80x. Nothing about
    "chance for perturbed TFs is about 12%, not 10%" changes. The file is NOT adjusted; the
    validation table prints the row as FAIL, and the final assert requires this mismatch and
    ONLY this mismatch, so it breaks if the slip is ever silently fixed in the file or if any
    other number starts to drift.

JUDGMENT CALLS: (a) `n_studymates` counts, per row, the OTHER rows of the 907 that share its
Profile.ID and perturb a DIFFERENT TF, which is studymate.py's mate_pairs set; (b)
`max_jaccard` is gapcheck.py's off-diagonal maximum over all 766 scored TFs, computed by a
matrix product instead of its 766^2 python loop (same numbers, seconds instead of minutes);
(c) `frac_regulon_default` is over PRUNED regulon edges, F-017's "scored regulon"; (d) rows
whose perturbed TF is not a scored source are DROPPED (they have no rank), so the file is 643
rows, not 907; (e) `self_edge_sign_decision` is read from data/collectri.parquet's
sign_decision column (V-016) and is NA when the TF has no self-edge in the pruned network.
"""
import sys, pathlib

pathlib.Path("../results").mkdir(exist_ok=True)
_log = open("../results/master.log", "w", buffering=1)
class _Tee:
    def write(self, t): sys.__stdout__.write(t); _log.write(t)
    def flush(self): sys.__stdout__.flush(); _log.flush()
sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd, scipy.sparse as sps

ALPHA, DIRECTION, NORM = 0.5, "up", "sum"          # F-004's winning propagation setting

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
tf_of = adata.obs.source.astype(str).values
type_p = adata.obs.type_p.values.astype(float)
logfc = adata.obs.logFC.values.astype(float)

def score(m, adj, rows=None):                      # selfedge.py's wrapper around decoupler's kernel
    return dc.mt.ulm.func(m if rows is None else m[rows], adj)[0]

print(f"\n{'='*96}\nF-029  master record: one row per scored knockTF experiment\n{'='*96}")
print(f"{len(obs)} experiments in data/knocktf_full.h5ad, {adata.obs.source.nunique()} perturbed TFs, "
      f"{len(var)} genes; {nS} TFs scored by CollecTRI at tmin=5")

# ------------------------------------------------------------------ the seven score conditions
mat_ma = mat.copy()                                            # leak.py: ASYMMETRIC own-gene mask
for i, t in enumerate(tf_of):
    if t in var_i: mat_ma[i, var_i[t]] = 0
mat_sym = mat.copy()                                           # studymate.py (7): SYMMETRIC mask
pert_cols = sorted({var_i[t] for t in set(tf_of) if t in var_i})
mat_sym[:, pert_cols] = 0

adjm_ns = adjm.copy()                                          # selfedge.py: every self-edge gone
self_rc = [(var_i[s], src_i[s]) for s in sources if s in var_i and adjm[var_i[s], src_i[s]] != 0]
for r, c in self_rc: adjm_ns[r, c] = 0

es_full = pd.DataFrame(score(mat, adjm), index=obs, columns=sources)
es_ma = pd.DataFrame(score(mat_ma, adjm), index=obs, columns=sources)
es_sym = pd.DataFrame(score(mat_sym, adjm), index=obs, columns=sources)
es_ns = pd.DataFrame(score(mat, adjm_ns), index=obs, columns=sources)

exps = [(e, t) for e, t in adata.obs.source.items() if t in src_i]
print(f"self-edges in the pruned network: {len(self_rc)} of {nS} TFs; "
      f"{len(pert_cols)} perturbed-TF gene columns blanked by the symmetric mask; "
      f"{len(exps)} scored experiments over {len({t for _, t in exps})} TFs")

# selfedge.py: only the perturbed TF's OWN column loses its self-edge, every other column intact
es_noself = es_full.copy()
for e, t in exps: es_noself.loc[e, t] = es_ns.loc[e, t]

# selfedge.py: gene X dropped from every OTHER TF's regulon, X->X kept; one network per TF
es_noother = es_full.copy()
for t in sorted({t for _, t in exps}):
    if t not in var_i: continue
    rows = [obs_i[e] for e, tt_ in exps if tt_ == t]
    g, keep = var_i[t], adjm[var_i[t]].copy()
    adjm[g] = 0; adjm[g, src_i[t]] = keep[src_i[t]]
    es_noother.loc[obs[rows]] = pd.DataFrame(score(mat, adjm, rows), index=obs[rows], columns=sources)
    adjm[g] = keep

# propagate.py: one-step propagation along CollecTRI TF->TF edges, autoregulation dropped
tfs = pd.Index(sources)
tt = net[net.source.isin(tfs) & net.target.isin(tfs)]
A = (tt.pivot_table(index="source", columns="target", values="weight", aggfunc="first")
       .reindex(index=tfs, columns=tfs).fillna(0))
A = A.mask(np.eye(len(tfs), dtype=bool), 0)
def propagate(S):                                   # direction "up", norm "sum" -> P = A, unnormalised
    return pd.DataFrame(S.values + ALPHA * (S.values @ A.values), index=S.index, columns=S.columns)

COND = {"ulm_unmasked": es_full, "ulm_masked": es_ma, "ulm_masked_symmetric": es_sym,
        "ulm_noselfedge": es_noself, "ulm_noother": es_noother,
        "propagated_unmasked": propagate(es_full), "propagated_masked": propagate(es_ma)}

def rankfrac(S):                                    # type_p * score, descending, over all nS TFs
    return pd.DataFrame(S.values * type_p[:, None], index=S.index,
                        columns=S.columns).rank(axis=1, ascending=False).values / nS

RF = {k: rankfrac(v) for k, v in COND.items()}

# ------------------------------------------------------------------ covariates
sd = (net.assign(_k=net.source + "|" + net.target)
         .drop_duplicates("_k").set_index("_k").sign_decision)
pn = pnet.assign(_k=pnet.source + "|" + pnet.target)
pn["sign_decision"] = pn._k.map(sd)
frac_def = pn.groupby("source").sign_decision.apply(lambda s: float((s == "default activation").mean()))
reg_size = pd.Series({s: int((adjm[:, c] != 0).sum()) for s, c in src_i.items()})

# gapcheck.py's max off-diagonal regulon Jaccard, by matrix product instead of its 766^2 loop
gt = net[net.source.isin(tfs)]
tgt = pd.Index(sorted(set(gt.target)))
B = np.zeros((len(tgt), nS), np.float32)
B[tgt.get_indexer(gt.target), tfs.get_indexer(gt.source)] = 1.0
inter = B.T @ B
sz = np.diag(inter).copy()
union = sz[:, None] + sz[None, :] - inter
J = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)
np.fill_diagonal(J, 0.0)
max_jac = pd.Series(J.max(1), index=tfs)

# studymate.py (1): the study is Profile.ID, the GEO series
study = adata.obs["Profile.ID"].astype(str).values
by_study = {}
for j, s in enumerate(study): by_study.setdefault(s, []).append(j)
n_mates = np.array([sum(tf_of[j] != tf_of[i] for j in by_study[study[i]]) for i in range(len(obs))])

# studymate.py: V-022's per-TF base rate, over the ASYMMETRIC masked ranks, all 907 rows
R_ma = RF["ulm_masked"]
base = {t: (float(np.mean(R_ma[tf_of != t, src_i[t]] <= .10)),
            float(np.mean(R_ma[tf_of != t, src_i[t]] <= .05)))
        for t in sorted({t for _, t in exps})}

# ------------------------------------------------------------------ assemble
rows = [obs_i[e] for e, _ in exps]
tf_r = np.array([t for _, t in exps])
D = pd.DataFrame(dict(
    exp=[e for e, _ in exps], tf=tf_r,
    in_filtered_279=logfc[rows] < -1,
    profile_id=study[rows],
    biosample=adata.obs["Biosample.Name"].astype(str).values[rows],
    knock_method=adata.obs["Knock.Method"].astype(str).values[rows],
    platform=adata.obs["Platform"].astype(str).values[rows],
    own_logfc=logfc[rows],
    own_gene_measured=[t in var_i for t in tf_r],
    regulon_size=reg_size.reindex(tf_r).values,
    has_self_edge=[t in var_i and adjm[var_i[t], src_i[t]] != 0 for t in tf_r],
    self_edge_weight=[float(adjm[var_i[t], src_i[t]]) if t in var_i else 0.0 for t in tf_r],
    frac_regulon_default=frac_def.reindex(tf_r).values,
    n_studymates=n_mates[rows],
    max_jaccard=max_jac.reindex(tf_r).values,
    base_rate10=[base[t][0] for t in tf_r],
    base_rate5=[base[t][1] for t in tf_r],
))
D["self_edge_weight"] = D.self_edge_weight.where(D.has_self_edge)
D["self_edge_sign_decision"] = np.where(D.has_self_edge, sd.reindex(tf_r + "|" + tf_r).values, None)
for k, R in RF.items():
    D[k + "_frac"] = R[rows, [src_i[t] for t in tf_r]]
    D[k + "_hit10"] = D[k + "_frac"] <= .10
    D[k + "_hit5"] = D[k + "_frac"] <= .05
D = D.sort_values(["tf", "exp"]).reset_index(drop=True)
D.to_csv("../results/master.csv", index=False)
print(f"\nwrote ../results/master.csv  {D.shape[0]} rows x {D.shape[1]} columns "
      f"({int(D.in_filtered_279.sum())} in the filtered 279)")
print("columns: " + ", ".join(D.columns))

# ------------------------------------------------------------------ VALIDATION, the point of this
F = D[D.in_filtered_279]                      # the 279
S_ = D[D.has_self_edge & D.in_filtered_279]   # the 197
pos = S_[S_.self_edge_weight > 0]             # the 176
neg = S_[S_.self_edge_weight < 0]             # the 21
V = []
def chk(entry, what, exp, got, tol=0.0):
    V.append(dict(entry=entry, quantity=what, published=exp, recomputed=got,
                  ok=bool(abs(float(got) - float(exp)) <= tol)))

chk("F-001", "ULM unmasked top-10%, 279", 0.452, round(F.ulm_unmasked_hit10.mean(), 3))
chk("F-001", "ULM unmasked hits / 279", 126, int(F.ulm_unmasked_hit10.sum()))
chk("F-001", "n experiments (filtered)", 279, len(F))
chk("F-001", "n TFs (filtered)", 155, F.tf.nunique())
chk("F-005/F-013", "masked top-10%, 279", 0.344, round(F.ulm_masked_hit10.mean(), 3))
chk("F-005/F-013", "masked hits / 279", 96, int(F.ulm_masked_hit10.sum()))
chk("F-010", "no-self-edge hits / 279", 96, int(F.ulm_noselfedge_hit10.sum()))
chk("F-010", "no-other-regulon hits / 279", 127, int(F.ulm_noother_hit10.sum()))
chk("F-010", "self-edge experiments", 197, len(S_))
chk("F-010", "no-self-edge experiments", 82, len(F) - len(S_))
chk("F-010", "self-edge group, unmasked top-10%", 0.538, round(S_.ulm_unmasked_hit10.mean(), 3))
chk("F-010", "self-edge group, masked top-10%", 0.386, round(S_.ulm_masked_hit10.mean(), 3))
NS = F[~F.has_self_edge]
chk("F-010", "no-self-edge group, unmasked top-10%", 0.244, round(NS.ulm_unmasked_hit10.mean(), 3))
chk("F-010", "no-self-edge group, masked top-10%", 0.244, round(NS.ulm_masked_hit10.mean(), 3))
chk("F-011", "+1 self-edge experiments", 176, len(pos))
chk("F-011", "-1 self-edge experiments", 21, len(neg))
chk("F-011", "+1 group, unmasked hits", 98, int(pos.ulm_unmasked_hit10.sum()))
chk("F-011", "+1 group, no-self-edge hits", 65, int(pos.ulm_noselfedge_hit10.sum()))
chk("F-011", "-1 group, unmasked hits", 8, int(neg.ulm_unmasked_hit10.sum()))
chk("F-011", "-1 group, no-self-edge hits", 11, int(neg.ulm_noselfedge_hit10.sum()))
vc = pos.self_edge_sign_decision.value_counts()
chk("F-016/F-017", "+1 self-edge, default activation", 92, int(vc.get("default activation", 0)))
chk("F-016/F-017", "+1 self-edge, PMID", 74, int(vc.get("PMID", 0)))
chk("F-016/F-017", "+1 self-edge, regulon", 10, int(vc.get("regulon", 0)))
chk("F-022", "scored experiments (unfiltered)", 643, len(D))
chk("F-022", "TFs scored (unfiltered)", 263, D.tf.nunique())
chk("F-022", "masked top-10%, all 643", 0.299, round(D.ulm_masked_hit10.mean(), 3))
chk("V-022", "per-TF base rate, 279", 0.123, round(F.base_rate10.mean(), 3))
chk("V-022", "per-TF base rate, 643 (where 0.123 is)", 0.123, round(D.base_rate10.mean(), 3))
chk("V-022", "per-TF base rate top-5%, 279", 0.067, round(F.base_rate5.mean(), 3))
chk("F-004", "propagated unmasked top-10%, 279", 0.566, round(F.propagated_unmasked_hit10.mean(), 3), 1e-3)
chk("F-004", "propagated masked top-10%, 279", 0.276, round(F.propagated_masked_hit10.mean(), 3), 1e-3)

VT = pd.DataFrame(V)
print(f"\n{'-'*96}\nVALIDATION: published numbers recomputed from master.csv alone\n{'-'*96}")
print(f"{'entry':14s} {'quantity':40s} {'published':>10s} {'recomputed':>11s}  ok")
for r in V:
    print(f"{r['entry']:14s} {r['quantity']:40s} {r['published']:>10} {r['recomputed']:>11}  "
          f"{'.' if r['ok'] else 'FAIL'}")
bad = VT[~VT.ok]
print(f"\n{int(VT.ok.sum())} of {len(VT)} reproduce."
      + ("" if bad.empty else "  MISMATCHES: " + "; ".join(
          f"{r.entry} {r.quantity} {r.published} -> {r.recomputed}" for r in bad.itertuples())))
VT.to_csv("../results/master_validation.csv", index=False)

# the one known, diagnosed mismatch (see OUTCOME in the docstring). Kept as a FAIL row on
# purpose: master.csv is not adjusted to it, and the assert below breaks if it goes away.
KNOWN = {("V-022", "per-TF base rate, 279")}
print(f"\nKNOWN MISMATCH, not a bug in this file: V-022's per-TF base rate on the FILTERED 279 is "
      f"{F.base_rate10.mean():.4f} = 0.122, not the 0.123 the Index, the V-022 entry and CLAUDE.md "
      f"quote.\n  0.123 is this same construction on all 643 SCORED rows ({D.base_rate10.mean():.4f}); "
      f"V-022 subtracts one figure from both rows.\n  The top-5% companion on the 279 reproduces "
      f"exactly ({F.base_rate5.mean():.4f} -> 0.067), and studymate.py already prints 0.122 here "
      f"(V-028: 'global 0.122'),\n  so the construction is V-022's and the log carries two values "
      f"for one number. Effect: masked ULM is +{F.ulm_masked_hit10.mean()-F.base_rate10.mean():.3f} "
      f"over base on the 279, {F.ulm_masked_hit10.mean()/F.base_rate10.mean():.2f}x, not +0.221 / 2.80x.")

# context the entries quote but do not pin down; printed, not asserted
print(f"\nunasserted context: symmetric-mask top-10% {F.ulm_masked_symmetric_hit10.mean():.3f} on the "
      f"279 and {D.ulm_masked_symmetric_hit10.mean():.3f} on the 643 (F-028's convention, NOT the "
      f"0.344 baseline's)\n  median frac_regulon_default over the 279 {F.frac_regulon_default.median():.2f} "
      f"(F-017: 0.53, IQR 0.15-0.62; here {F.frac_regulon_default.quantile(.25):.2f}-"
      f"{F.frac_regulon_default.quantile(.75):.2f})\n  -1 self-edge sign_decision on the 21: "
      + ", ".join(f"{k} {v}" for k, v in neg.self_edge_sign_decision.value_counts().items())
      + f" (F-016: 20 PMID, 1 regulon)\n  top-5%: unmasked {F.ulm_unmasked_hit5.mean():.3f} masked "
        f"{F.ulm_masked_hit5.mean():.3f} (F-005: 0.280); rows with >=1 study-mate "
        f"{int((D.n_studymates > 0).sum())} of {len(D)}")
print(f"{'='*96}")

assert set(zip(bad.entry, bad.quantity)) == KNOWN, (
    "the set of published numbers master.csv fails to reproduce changed. expected exactly "
    + str(KNOWN) + ", got " + str({(r.entry, r.quantity, r.published, r.recomputed)
                                   for r in bad.itertuples()}))
