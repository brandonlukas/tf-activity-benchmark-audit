"""T-028 / F-017: does the SOURCE of a CollecTRI self-edge's sign change what F-011 measured?

F-016 found data/collectri.parquet carries a `sign_decision` column ("PMID", "regulon",
"default activation"). Every -1 edge is evidence-based; 58% of +1 edges are the default fill
("TF-gene interactions without any information were assigned an activating mode of regulation").
F-011's +1 self-edge group (176 experiments, 84 TFs) is therefore a mix of curated activating
self-edges and edges with no sign information at all. This script:

  (1) commits F-016's crosstabs (all 42990 edges, 461 self-edges, 96 knockTF TFs, 279 experiments);
  (2) splits F-011's conditions (full network / self-edge removed / own gene masked) by
      sign_decision, and attributes the 33 top-10% hits ULM loses on self-edge removal (V-011:
      all 33 are +1) to default vs PMID vs regulon;
  (3) size control: median rank change on self-edge removal within regulon-size bins for
      default vs PMID +1, because V-010 says the self-edge effect is a size effect and
      default-filled signs should sit on less-studied, smaller regulons;
  (4) whole-network accounting: default-activation fraction of each perturbed TF's scored
      regulon (median, IQR), and masked-ULM hit rate for mostly-evidence-signed (>= 50%
      PMID/regulon) vs mostly-default regulons, within the same size bins.

Conditions, ULM kernel and evaluate() are selfedge.py's, imported, not reimplemented (selfsign.py
is being edited by someone else, so the sign split is recomputed here from adjm_orig the same way).
No gain over ULM is claimed, so no shuffled network: "no self-edge" and "masked" are the controls
and both are in every table.

PREDICTIONS (stated before running):
  (2) the 33 lost hits split roughly in proportion to the +1 group sizes (92 default / 74 PMID /
      10 regulon), i.e. about 17 default / 14 PMID / 2 regulon; the top-10% drop on self-edge
      removal is similar in the default and PMID groups once regulon size is matched;
  (3) provenance does NOT matter beyond size -- by the algebra a +1 is a +1 whatever labelled it,
      so within a size bin the median rank change should be the same for default and PMID;
  (4) masked ULM top-10% higher for mostly-evidence-signed regulons than for mostly-default,
      about 0.40 vs 0.28, and mostly explained by regulon size.
Point (2)/(3) is accounting, not mechanism: nothing here can distinguish "the sign is right" from
"the sign is a guess", only which bucket the benchmark's hits sit in.

Runtime ~40 s (imports selfedge.py, which reruns propagate.py, controls.py and leak.py).
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/signsource.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import numpy as np, pandas as pd
from selfedge import (adata, adjm_orig, es_full, es_ns, exps, fr, has_self, pnet, src_i,
                      tfs, var_i)   # reruns F-001/F-004/F-006/F-010's tables

net = pd.read_parquet("../data/collectri.parquet")
self_rows = net[net.source == net.target]

# ---------------- (1) F-016's crosstabs, committed ----------------
tf_exp = pd.DataFrame({"exp": [e for e, _ in exps], "tf": [t for _, t in exps]})
sd_self = self_rows.set_index("source")                       # one row per self-edge TF
knock_tf = sd_self.reindex(sorted(set(tf_exp.tf))).dropna(subset=["weight"])
knock_exp = tf_exp.join(sd_self, on="tf").dropna(subset=["weight"])

def xt(df, title):
    t = pd.crosstab(df.sign_decision, df.weight)
    print(f"\n{title} (n={len(df)})")
    print(t.to_string())

print(f"\n{'='*78}\n(5) T-028: CollecTRI sign_decision crosstabs (F-016, committed)")
xt(net, "all edges")
xt(self_rows, "self-edges (source == target)")
xt(knock_tf, f"self-edges of the {len(set(tf_exp.tf))} knockTF perturbed TFs, one row per TF")
xt(knock_exp, f"same, one row per experiment ({len(exps)} experiments)")
print(f"  {len(exps) - len(knock_exp)} experiments whose perturbed TF has no self-edge in the raw network")

# ---------------- per-experiment labels ----------------
w_self = np.array([adjm_orig[var_i[tf], src_i[tf]] if tf in var_i else 0.0 for _, tf in exps])
sd = sd_self.sign_decision.to_dict()
label = np.array([f"{sd[tf]} ({w:+.0f})" if h else "no self-edge"
                  for (_, tf), w, h in zip(exps, w_self, has_self)])
reg_size = np.array([(adjm_orig[:, src_i[tf]] != 0).sum() for _, tf in exps])
n_tfs_of = lambda m: len({t for (_, t), k in zip(exps, m) if k})

VAR = ["ulm (full net)", "no self-edge", "masked (leak.py)"]
ref = fr["ulm (full net)"]
GROUPS = ["default activation (+1)", "PMID (+1)", "regulon (+1)",
          "PMID (-1)", "regulon (-1)", "no self-edge"]

def block(title, m):
    print(f"\n{title}: n={m.sum()} experiments, {n_tfs_of(m)} TFs, "
          f"median regulon {np.median(reg_size[m]) if m.sum() else float('nan'):.0f} targets")
    if m.sum() == 0: return
    print(f"  {'variant':22s} top10  hits/n     top5  hits/n    better worse")
    for k in VAR:
        f, r = fr[k][m], ref[m]
        print(f"  {k:22s} {np.mean(f<=.1):.3f} {np.sum(f<=.1):3d}/{m.sum():<4d} "
              f"{np.mean(f<=.05):.3f} {np.sum(f<=.05):3d}/{m.sum():<4d} "
              f"{np.mean(f<r):.2f}   {np.mean(f>r):.2f}")

print(f"\n{'='*78}\n(6) F-011's conditions split by sign_decision "
      f"(better/worse = that variant vs unmasked full-network ULM, within the group)")
for g in GROUPS: block(g, label == g)

# ---------------- (2) where do the 33 lost hits sit? ----------------
hit_full, hit_ns = ref <= .1, fr["no self-edge"] <= .1
lost, gained = hit_full & ~hit_ns, ~hit_full & hit_ns
print(f"\n  top-10% hits lost / gained on self-edge removal, by sign_decision "
      f"(V-011: 33 lost, all +1; 3 gained, all -1)")
print(f"  {'group':26s} lost gained  hits full -> no self-edge")
for g in GROUPS:
    m = label == g
    print(f"  {g:26s} {np.sum(lost&m):4d} {np.sum(gained&m):6d}  "
          f"{np.sum(hit_full&m):3d} -> {np.sum(hit_ns&m):3d}")
print(f"  {'TOTAL':26s} {lost.sum():4d} {gained.sum():6d}  {hit_full.sum():3d} -> {hit_ns.sum():3d}")

# ---------------- (3) size control: default vs PMID +1 ----------------
d_rank = (fr["no self-edge"] - ref) * len(tfs)      # positive = worse rank after removal
bins = [("<= 60", reg_size <= 60), ("61-150", (reg_size > 60) & (reg_size <= 150)),
        ("> 150", reg_size > 150)]
print(f"\n  (3) median rank change of the true TF on self-edge removal (of {len(tfs)} TFs; "
      f"+ = pushed down), within regulon-size bins")
print(f"  {'size bin':10s} {'default (+1)':>28s} {'PMID (+1)':>28s}")
print(f"  {'':10s} {'n  TFs  med.size  drank  hits':>28s} {'n  TFs  med.size  drank  hits':>28s}")
for bl, bm in bins:
    cells = []
    for g in ["default activation (+1)", "PMID (+1)"]:
        m = (label == g) & bm
        cells.append(f"{m.sum():2d} {n_tfs_of(m):4d} {np.median(reg_size[m]) if m.sum() else np.nan:9.0f} "
                     f"{np.median(d_rank[m]) if m.sum() else np.nan:6.0f} "
                     f"{np.sum(hit_full&m):3d}->{np.sum(hit_ns&m):<3d}" if m.sum() else f"{0:2d}" + " " * 26)
    print(f"  {bl:10s} {cells[0]:>28s} {cells[1]:>28s}")

# the algebraic quantity behind the rank change: how far the self-edge moves the TF's own score.
# If provenance is irrelevant (a +1 is a +1), this matches at matched size and any rank-change
# difference is a starting-rank effect (V-011's finding for the -1 group).
d_score = np.array([abs(es_full.loc[e, tf] - es_ns.loc[e, tf]) for e, tf in exps])
own_lfc = np.array([abs(adata.to_df().loc[e, tf]) if tf in adata.var_names else np.nan
                    for e, tf in exps])
print(f"  same bins: median |d own score| (t-units), median |own logFC|, median starting rank "
      f"fraction under full ULM")
print(f"  {'size bin':10s} {'default (+1): |ds|  |lfc|  start':>34s} {'PMID (+1): |ds|  |lfc|  start':>34s}")
for bl, bm in bins:
    cells = []
    for g in ["default activation (+1)", "PMID (+1)"]:
        m = (label == g) & bm
        cells.append(f"{np.median(d_score[m]):.3f}  {np.median(own_lfc[m]):.2f}   "
                     f"{np.median(ref[m]):.3f}" if m.sum() else "-")
    print(f"  {bl:10s} {cells[0]:>34s} {cells[1]:>34s}")

# ---------------- (4) whole-regulon sign provenance ----------------
# dc.pp.prune drops the extra columns, so put sign_decision back by (source, target)
scored = pnet[pnet.source.isin(src_i)].merge(net[["source", "target", "sign_decision"]],
                                             on=["source", "target"], how="left")
frac_def = scored.groupby("source").sign_decision.apply(lambda s: (s == "default activation").mean())
fd = np.array([frac_def[tf] for _, tf in exps])
q1, med, q3 = np.percentile(fd, [25, 50, 75])
print(f"\n  (4) default-activation fraction of the SCORED regulon (pruned network) of the "
      f"{len(exps)} experiments' perturbed TFs:\n      median {med:.2f}, IQR {q1:.2f}-{q3:.2f} "
      f"({len(set(tf_exp.tf))} distinct TFs: median {np.median(frac_def[sorted(set(tf_exp.tf))]):.2f})")
ev, df_ = fd < .5, fd >= .5
print(f"  masked-ULM hit rate by regulon sign provenance (size is confounded; bins given)")
print(f"  {'stratum':22s} {'mostly evidence (<50% default)':>32s} {'mostly default (>=50%)':>32s}")
print(f"  {'':22s} {'n  TFs med.size top10 hits/n':>32s} {'n  TFs med.size top10 hits/n':>32s}")
for bl, bm in [("all", np.ones(len(exps), bool))] + bins:
    cells = []
    for m0 in (ev, df_):
        m = m0 & bm
        f = fr["masked (leak.py)"][m]
        cells.append(f"{m.sum():3d} {n_tfs_of(m):4d} {np.median(reg_size[m]) if m.sum() else np.nan:8.0f} "
                     f"{np.mean(f<=.1) if m.sum() else np.nan:5.3f} {np.sum(f<=.1):3d}/{m.sum():<4d}")
    print(f"  {bl:22s} {cells[0]:>32s} {cells[1]:>32s}")
print("  NOTE: accounting only. The SELF-EDGE's label tracks how well studied a TF is (median regulon"
      "\n  78 targets for default vs 186 for PMID), but the WHOLE-REGULON default fraction in this block"
      "\n  does not: TF-level Spearman(default fraction, regulon size) = -0.02, p = 0.81, and the two"
      "\n  halves' median sizes are 75 and 79 (V-017 (5)). Block 4 is therefore not confounded by size,"
      "\n  it is flat: the difference is +0.014 with TF-cluster 95% CI [-0.14, +0.17]. No mechanism is"
      "\n  claimed and 'sign quality does not limit ULM' is not established, only 'cannot tell at this n'.")

# F-029: persist. One row per experiment with its self-edge sign_decision label, the whole-regulon
# default fraction and every variant's rank fraction. No printed line changes.
pd.DataFrame(dict(exp=[e for e, _ in exps], tf=[t for _, t in exps], sign_decision_label=label,
                  self_edge_weight=w_self, has_self_edge=has_self, regulon_size=reg_size,
                  frac_regulon_default=fd, hit_full=hit_full, hit_noselfedge=hit_ns,
                  **{k.split(" (")[0].replace(" ", "_"): v for k, v in fr.items()})
             ).to_csv("../results/signsource.csv", index=False)

assert (len(net), len(self_rows), len(knock_exp) + (~has_self).sum()) == (42990, 461, 279) \
    and net.sign_decision.value_counts().reindex(
        ["PMID", "regulon", "default activation"]).tolist() == [19396, 1895, 21699] \
    and sum(np.sum(hit_full & (label == g)) for g in GROUPS) == hit_full.sum() == 126 \
    and sum(np.sum(hit_ns & (label == g)) for g in GROUPS) == hit_ns.sum() == 96, \
    "crosstab totals or the sign_decision grouping do not reproduce F-016 / F-010's hit counts"
