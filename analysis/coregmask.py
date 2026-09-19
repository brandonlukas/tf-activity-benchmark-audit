"""T-012 / F-012: the mask controls V-007 had to run from scratch, as a committed script.

randmask.py (F-007) ships two arms that V-007 showed are near no-ops: 71% of measured genes sit in
no CollecTRI regulon at all, so most "any" picks cannot change a ULM score, and the regulon-count
"matched" pick zeroes a median |logFC| of ~0.13 where the real mask zeroes a ~5 sd value that is on
the leak path by construction. This is randmask.py copied (F-007 is not edited, and its printed
table must stay as V-007 verified it) with its two arms kept for reference and two added that do sit
on the path:

  coreg   the measured gene, other than the perturbed TF, with the MOST NEGATIVE logFC among the 100
          genes whose CollecTRI upstream-regulator set is most Jaccard-similar to the perturbed TF's
          own gene  (V-007's "most negative among the top 100")
  logfc   the measured gene whose logFC is NEAREST the perturbed TF's own logFC, among genes sitting
          in at least (regulons of the TF's gene minus 3) CollecTRI regulons

Both new picks are deterministic given the TF and the experiment, so the seed loop only re-derives
them and the [min-max] ranges are zero-width by construction; the seeds are kept so the assert covers
"every seed and arm" and so the any/matched arms reproduce F-007 exactly.

Reported for every arm: median logFC and median |logFC| of the genes actually zeroed (the point of
T-012 -- an arm that zeroes ~0 is not a control), plus the regulon gap of the matched arm on the 279
SCORED experiments (V-007: 2.69) next to the 2.05 randmask.py prints over all 388.

No gain over ULM is claimed here; the controls the convention owes are both in the output anyway,
via the imports: leak.py's perturbed-TF-masked table and controls.py's 10-draw degree-preserving
shuffled network (the SHUF rows).

PREDICTION (stated before running, from V-007's scratch numbers):
  coreg   ulm 0.44,  up/sum/a=0.5 0.56,   SHUF 0.29
  logfc   ulm 0.45,  up/sum/a=0.5 0.565,  SHUF 0.30
i.e. neither arm comes near the perturbed-TF mask (0.34 ulm / 0.28 propagated): the 57->28 collapse
is the TF's own mRNA, not "any co-regulated gene with a large negative logFC".

--- T-032 / F-020 additions (V-012's owed items; appended so every block above prints
byte-identically to the run V-012 checked) ---------------------------------------------------
The masked-gene table and the shared-regulator mean above are pooled over all 388 experiments,
but only 279 are SCORED (the other 109 have no ULM score for their TF, F-013), so the numbers
in the table are not the numbers behind the hit rates. Both are reprinted on the 279. Two
variant arms then test the two judgment calls V-012 named in the picker:
  coreg-v  SKIPS the experiment (no cell zeroed) when the perturbed TF's gene has no CollecTRI
           regulator at all -- the coreg arm fills those with a gene that is not co-regulated,
           since every Jaccard is 0 and argsort takes the first -- and steps past already-zero
           cells to the next candidate, so a pick is always a real zeroing.
  logfc-v  regulon-count floor max(memb - 3, 1) instead of max(memb - 3, 0): with the 0 floor
           the constraint is vacuous for TF genes in <= 3 regulons, so genes in NO regulon are
           eligible and a pick can be a no-op for ULM. Also steps past already-zero cells, so
           the extended assert ("every variant pick is nonzero BEFORE masking") covers both.

PREDICTION for the variants (stated before running): both land within 0.005 of the arms they
copy, i.e. coreg-v ~ coreg 0.437 ulm / 0.563 propagated / 0.292 shuffled (V-007's scratch
0.441 / 0.563 / 0.292) and logfc-v ~ logfc 0.444 / 0.566 / 0.300 (V-007's 0.452 / 0.565 /
0.304). The 53 experiments with no regulator (19 of them scored) and the handful of already-zero
picks are too few to move a rate over 279 experiments; if either variant moves by more than
0.01 the original arm was being carried by its filler picks and F-012's arms need rewording.
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/coregmask.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import anndata as ad, decoupler as dc, numpy as np, pandas as pd
from propagate import A, propagate, tfs          # reruns the propagation grid
from controls import shuffle_preserving_degree, Ashuf   # reruns unmasked table (incl. 10-seed shuffle)
import leak                                      # reruns the perturbed-TF-masked table

adata = ad.read_h5ad("../data/knocktf.h5ad"); net = pd.read_parquet("../data/collectri.parquet")
X0 = adata.to_df(); genes = X0.columns
memb = net.groupby("target").size().reindex(genes).fillna(0)   # regulons each measured gene sits in

# upstream-regulator sets: M[gene, source] = 1 if that TF regulates the gene (all of CollecTRI, not
# just the TF-TF part), used for the Jaccard that defines the coreg arm
srcs = pd.Index(sorted(net.source.unique()))
sub = net[net.target.isin(genes)]
M = np.zeros((len(genes), len(srcs)), np.float32)
M[genes.get_indexer(sub.target), srcs.get_indexer(sub.source)] = 1.0
cnt = M.sum(1)
regs = net.groupby("target").source.unique()

TOPJ, SHARED = {}, {}
def top_jaccard(tf, k=100):
    """Indices of the k measured genes whose regulator set is most Jaccard-similar to gene `tf`'s."""
    if tf not in TOPJ:
        v = np.zeros(len(srcs), np.float32)
        if tf in regs.index: v[srcs.get_indexer(regs[tf])] = 1.0
        inter = M @ v
        jac = inter / np.maximum(cnt + v.sum() - inter, 1e-9)
        if tf in genes: jac[genes.get_loc(tf)] = -1.0          # never the perturbed TF itself
        TOPJ[tf] = np.argsort(-jac, kind="stable")[:k]
        SHARED[tf] = inter[TOPJ[tf]]
    return TOPJ[tf]

def evaluate(S):                                  # same statistic as propagate.py / leak.py
    fr = []
    for exp, tf in adata.obs.source.items():
        if tf not in S.columns: continue
        s = adata.obs.loc[exp, "type_p"] * S.loc[exp]; fr.append(s.rank(ascending=False)[tf] / S.shape[1])
    return np.array(fr)

def pick(exp, tf, rng, arm):
    """A measured gene that is not the perturbed TF."""
    if arm == "any":     return rng.choice(memb.index[memb.index != tf])
    if arm == "matched":                                        # randmask.py's picker, verbatim
        c = memb.drop(tf, errors="ignore").sub(memb.get(tf, 0)).abs().nsmallest(50).index
        return rng.choice(c)
    row = X0.loc[exp].values
    if arm == "coreg":
        c = top_jaccard(tf); return genes[c[np.argmin(row[c])]]
    own = row[genes.get_loc(tf)] if tf in genes else 0.0
    ok = cnt >= max(memb.get(tf, 0) - 3, 0)                     # logfc arm: regulon count >= TF's minus 3
    if tf in genes: ok[genes.get_loc(tf)] = False
    i = np.flatnonzero(ok)
    return genes[i[np.argmin(np.abs(row[i] - own))]]

scored = np.array([t in tfs for t in adata.obs.source])         # 279 of 388
rows, picks, zeroed = [], [], []
for arm, nseed in [("any", 5), ("matched", 5), ("coreg", 3), ("logfc", 3)]:
    for seed in range(nseed):
        rng = np.random.default_rng(seed)
        X = X0.copy(); ij = []
        for exp, tf in adata.obs.source.items():
            g = pick(exp, tf, rng, arm)
            picks.append(dict(arm=arm, seed=seed, exp=exp, tf=tf, gene=g, lfc=X0.loc[exp, g],
                              scored=tf in tfs))
            ij.append((X.index.get_loc(exp), genes.get_loc(g))); X.loc[exp, g] = 0
        r, c = np.array(ij).T
        t = np.array([(X.index.get_loc(e), genes.get_loc(tf)) for e, tf in adata.obs.source.items()
                      if tf in genes]).T                        # the TF's OWN cell must survive
        zeroed.append(bool((X.values[r, c] == 0).all() and (X.values[t[0], t[1]] == X0.values[t[0], t[1]]).all()))
        a = ad.AnnData(X, obs=adata.obs); dc.mt.ulm(a, net, tmin=5)
        S = a.obsm["score_ulm"].reindex(columns=tfs)
        fr0 = evaluate(S)
        P = pd.DataFrame(propagate(S, A, .5, "up", "sum"), index=S.index, columns=tfs)
        frp = evaluate(P)
        shuf = [evaluate(pd.DataFrame(propagate(S, B, .5, "up", "sum"), index=S.index, columns=tfs)) for B in Ashuf]
        rows.append(dict(variant=f"{arm} ulm", seed=seed, top10=np.mean(fr0 <= .1), top5=np.mean(fr0 <= .05),
                         better=0.0, worse=0.0))
        rows.append(dict(variant=f"{arm} up sum a=0.5", seed=seed, top10=np.mean(frp <= .1), top5=np.mean(frp <= .05),
                         better=np.mean(frp < fr0), worse=np.mean(frp > fr0)))
        rows.append(dict(variant=f"{arm} up sum a=0.5 SHUF", seed=seed,
                         top10=np.mean([np.mean(f <= .1) for f in shuf]), top5=np.mean([np.mean(f <= .05) for f in shuf]),
                         better=np.mean([np.mean(f < fr0) for f in shuf]), worse=np.mean([np.mean(f > fr0) for f in shuf])))

df = pd.DataFrame(rows)
g = df.groupby("variant", sort=False).agg(top10=("top10", "mean"), t10_lo=("top10", "min"), t10_hi=("top10", "max"),
                                          top5=("top5", "mean"), better=("better", "mean"), worse=("worse", "mean"))
print("\nfour mask arms (better/worse are paired vs the SAME-seed, SAME-arm masked ulm)")
print(g.round(3).to_string())
print("\nreference rows above: unmasked ulm 0.45 / up sum a=0.5 0.57 ; perturbed-TF-masked ulm 0.34 / 0.28")

P = pd.DataFrame(picks)
own_lfc = np.array([X0.loc[e, t] for e, t in adata.obs.source.items() if t in genes])
print(f"\nwhat each arm actually zeroes (all seeds pooled; the perturbed TF's own logFC is "
      f"median {np.median(own_lfc):.2f}, |{np.median(np.abs(own_lfc)):.2f}|, and negative in "
      f"{np.mean(own_lfc < 0):.0%} of experiments)")
print(f"{'arm':10s} {'n picks':>7s} {'median logFC':>13s} {'median |logFC|':>15s} {'frac exactly 0':>15s} "
      f"{'frac |lfc|>=1':>14s}")
for arm in ["any", "matched", "coreg", "logfc"]:
    v = P.lfc[P.arm == arm].values
    print(f"{arm:10s} {len(v):7d} {np.median(v):13.3f} {np.median(np.abs(v)):15.3f} "
          f"{np.mean(v == 0):15.2f} {np.mean(np.abs(v) >= 1):14.2f}")

m = P[(P.arm == "matched")]
gap_all = (memb.reindex(m.gene).values - memb.reindex(m.tf).fillna(0).values).__abs__()
gap_sc = (memb.reindex(m.gene[m.scored]).values - memb.reindex(m.tf[m.scored]).fillna(0).values).__abs__()
print(f"\nmatched arm regulon gap, mean |regulons(pick) - regulons(TF gene)|: {gap_all.mean():.2f} over all "
      f"{len(m)//5} experiments (what randmask.py prints), {gap_sc.mean():.2f} over the {m.scored.sum()//5} "
      f"SCORED ones (V-007: 2.69); median regulons(TF gene) {memb.reindex(m.tf).fillna(0).median():.0f} all, "
      f"{memb.reindex(m.tf[m.scored]).fillna(0).median():.0f} scored")
co = P[(P.arm == "coreg") & (P.seed == 0)]
sh = np.array([SHARED[t][list(TOPJ[t]).index(genes.get_loc(gg))] for t, gg in zip(co.tf, co.gene)])
print(f"coreg arm: the picked gene shares a mean of {sh.mean():.1f} upstream regulators with the perturbed "
      f"TF's gene (median regulators of the TF's gene {memb.reindex(co.tf).fillna(0).median():.0f})")

# ---- T-032: the same two tables on the 279 SCORED experiments only (V-012) -----------------
own_sc = np.array([X0.loc[e, t] for e, t in adata.obs.source.items() if t in genes and t in tfs])
print(f"\nwhat each arm zeroes, restricted to the {len(own_sc)} SCORED experiments (the rows the "
      f"hit rates are computed on; the table above pools all {len(adata)} experiments). The "
      f"perturbed TF's own logFC is median {np.median(own_sc):.3f} on these.")
print(f"{'arm':10s} {'n picks':>7s} {'median logFC':>13s} {'median |logFC|':>15s} {'frac exactly 0':>15s} "
      f"{'frac |lfc|>=1':>14s}")
for arm in ["any", "matched", "coreg", "logfc"]:
    v = P.lfc[(P.arm == arm) & P.scored].values
    print(f"{arm:10s} {len(v):7d} {np.median(v):13.3f} {np.median(np.abs(v)):15.3f} "
          f"{np.mean(v == 0):15.2f} {np.mean(np.abs(v) >= 1):14.2f}")
print(f"coreg arm on the scored experiments: mean {sh[co.scored.values].mean():.2f} shared upstream "
      f"regulators (all 388: {sh.mean():.2f}), median regulators of the TF's gene "
      f"{memb.reindex(co.tf[co.scored]).fillna(0).median():.0f}")
noreg = [t for t in adata.obs.source if t not in regs.index]
print(f"experiments whose perturbed TF's gene has NO CollecTRI regulator: {len(noreg)} of "
      f"{len(adata)}, {sum(t in tfs for t in noreg)} of them scored -- the coreg arm fills these "
      f"with a gene sharing 0 regulators; the coreg-v arm below skips them instead")

# ---- T-032: variant arms (skip instead of fill; nonzero picks; floor 1) --------------------
def pick_v(exp, tf, arm):
    """Variant pickers. Return None to SKIP the experiment (leave the row unmasked) rather than
    fall back on a pick that is not what the arm is supposed to be. Never returns an
    already-zero cell: zeroing a zero is a no-op, not a control."""
    row = X0.loc[exp].values
    if arm == "coreg-v":
        if tf not in regs.index: return None                 # no regulators -> every Jaccard is 0
        cand = top_jaccard(tf)
        order = cand[np.argsort(row[cand], kind="stable")]   # most negative logFC first
    else:
        own = row[genes.get_loc(tf)] if tf in genes else 0.0
        ok = cnt >= max(memb.get(tf, 0) - 3, 1)              # floor 1, not 0
        if tf in genes: ok[genes.get_loc(tf)] = False
        i = np.flatnonzero(ok)
        order = i[np.argsort(np.abs(row[i] - own), kind="stable")]
    nz = order[row[order] != 0]
    return genes[nz[0]] if len(nz) else None

vrows, vpicks, vzeroed = [], [], []
for arm in ["coreg-v", "logfc-v"]:
    X = X0.copy(); ij = []; skipped = []
    for exp, tf in adata.obs.source.items():
        gv = pick_v(exp, tf, arm)                    # not `g`: that is the summary frame above
        if gv is None:
            skipped.append(tf); continue
        vpicks.append(dict(arm=arm, exp=exp, tf=tf, gene=gv, lfc=X0.loc[exp, gv], scored=tf in tfs))
        ij.append((X.index.get_loc(exp), genes.get_loc(gv))); X.loc[exp, gv] = 0
    r, c = np.array(ij).T
    tcell = np.array([(X.index.get_loc(e), genes.get_loc(tf)) for e, tf in adata.obs.source.items()
                      if tf in genes]).T
    vzeroed.append(bool((X.values[r, c] == 0).all()
                        and (X.values[tcell[0], tcell[1]] == X0.values[tcell[0], tcell[1]]).all()))
    a = ad.AnnData(X, obs=adata.obs); dc.mt.ulm(a, net, tmin=5)
    S = a.obsm["score_ulm"].reindex(columns=tfs)
    fr0 = evaluate(S)
    Pr = pd.DataFrame(propagate(S, A, .5, "up", "sum"), index=S.index, columns=tfs)
    frp = evaluate(Pr)
    shuf = [evaluate(pd.DataFrame(propagate(S, B, .5, "up", "sum"), index=S.index, columns=tfs)) for B in Ashuf]
    vrows.append(dict(variant=f"{arm} ulm", top10=np.mean(fr0 <= .1), top5=np.mean(fr0 <= .05),
                      better=0.0, worse=0.0, masked=len(ij), skipped=len(skipped)))
    vrows.append(dict(variant=f"{arm} up sum a=0.5", top10=np.mean(frp <= .1), top5=np.mean(frp <= .05),
                      better=np.mean(frp < fr0), worse=np.mean(frp > fr0), masked=len(ij),
                      skipped=len(skipped)))
    vrows.append(dict(variant=f"{arm} up sum a=0.5 SHUF",
                      top10=np.mean([np.mean(f <= .1) for f in shuf]),
                      top5=np.mean([np.mean(f <= .05) for f in shuf]),
                      better=np.mean([np.mean(f < fr0) for f in shuf]),
                      worse=np.mean([np.mean(f > fr0) for f in shuf]), masked=len(ij),
                      skipped=len(skipped)))

V = pd.DataFrame(vpicks)
vg = pd.DataFrame(vrows).set_index("variant")
orig = g.loc[["coreg ulm", "coreg up sum a=0.5", "coreg up sum a=0.5 SHUF",
              "logfc ulm", "logfc up sum a=0.5", "logfc up sum a=0.5 SHUF"], "top10"]
vg["orig_top10"] = [orig[v.replace("-v", "")] for v in vg.index]
vg["delta"] = vg.top10 - vg.orig_top10
print("\nvariant arms: skip-instead-of-fill + never zero an already-zero cell (coreg-v), "
      "regulon floor max(memb-3, 1) (logfc-v). One pass each: picks are deterministic.")
print(vg.round(3).to_string())
print(f"{'arm':10s} {'n picks':>7s} {'median logFC':>13s} {'median |logFC|':>15s} {'frac exactly 0':>15s} "
      f"{'frac |lfc|>=1':>14s}   (scored picks only)")
for arm in ["coreg-v", "logfc-v"]:
    v = V.lfc[(V.arm == arm) & V.scored].values
    print(f"{arm:10s} {len(v):7d} {np.median(v):13.3f} {np.median(np.abs(v)):15.3f} "
          f"{np.mean(v == 0):15.2f} {np.mean(np.abs(v) >= 1):14.2f}")
print(f"largest |top10 - original arm| = {vg.delta.abs().max():.3f} -> "
      + ("both variants reproduce their arm (prediction: within 0.005)" if vg.delta.abs().max() <= .005
         else "a variant moved more than 0.005: the original arm leans on its filler picks"))
vg.to_csv("../results/coregmask.csv")                # F-029: the printed variant table
P.to_csv("../results/coregmask_picks.csv", index=False)   # which gene each arm zeroed

assert all(zeroed) and (P.tf != P.gene).all() and all(vzeroed) and (V.tf != V.gene).all() \
    and (V.lfc != 0).all(), \
    "every arm and seed must zero exactly the genes it picked, never the perturbed TF's own " \
    "cell, and every variant pick must be nonzero before it is masked"
