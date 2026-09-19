"""T-010 / F-027: does the abundance leak transfer to the PUBLISHED 2026 benchmark's own code and own metric?

Target: Zhu, Han & Wang 2026, Brief Bioinform 10.1093/bib/bbag513, repo
zhuzhe0011/Benchmarking-Methods-for-Inferring-Single-cell-Transcription-Factor-Activity.
Only public dataset: dataset/Inhibition/Dox, 617 cells x 15498 genes, TP53 perturbation
(197 perturbed vs 420 CTRL). We run THEIR methord/decoupler/dcp.py and THEIR result/auc/auc.py,
unpatched for the anchor and for the masked condition, on THEIR shipped filtered.h5ad.

Conditions (F-010's decomposition, transplanted onto their pipeline):
  anchor   original h5ad, original net_collectri_filter.csv
  (a) MASK TP53's own gene zeroed in the count matrix before their normalize_total/log1p
  (b) NOSELF TP53->TP53 CollecTRI edge deleted, gene left in the matrix
  (c) NOOTHER TP53 deleted as a TARGET from the other 108 TFs' regulons, self-edge kept
  SHUF x10 degree-preserving, sign-mix-preserving target shuffle of the non-self edges
          (self-edges left in place, so this is "random regulon + real autoregulation")

PREDICTIONS, stated before any condition was run (the anchor was already in flight; its
prediction is simply "reproduces their shipped number"):

 1. Anchor. ULM_CollecTRI AUROC for TP53__Dox reproduces their shipped 0.957240 to within
    0.01 (decoupler 2.2.0 here vs 2.1.1 in the paper); au_rank 0.0633, auprc 0.9231.
 2. (a) MASK. THE LEAK LARGELY DOES NOT TRANSFER: I predict the ULM AUROC falls by LESS
    THAN 0.02, i.e. stays in 0.94-0.96, and the bootstrap CI over cells covers the anchor.
    Mechanism, read off their inputs before running: TP53's CollecTRI regulon here is 980
    edges / 772 measured targets, and TP53's own transcript is exactly 0 in 83% of perturbed
    and 85% of control cells (raw count mean 0.198 pert vs 0.174 ctrl -- it does not even go
    down in counts). One dropout-dominated gene out of 772 cannot move a t-statistic. This is
    the opposite regime from knockTF, where the own gene is a ~5 sd bulk logFC in a median
    97-target regulon; V-010 already showed the self-edge is worth ~3% of hits for regulons
    of several hundred targets.
 3. (b) NOSELF equals (a) for ULM to within 0.01 and within bootstrap noise (univariate).
 4. (c) NOOTHER leaves TP53's OWN AUROC unchanged to < 1e-6 (ULM is univariate, column TP53
    does not contain TP53-as-target), while moving the 108 TFs that hold TP53 in their regulon.
 5. Internal null group (no self-edge AND TP53 not in their regulon, ~600 TFs): mean |dAUROC|
    < 1e-6 under (b) and (c); under (a) nonzero but < 0.005, because zeroing TP53 changes
    every cell's library size by ~0.2 of several thousand counts.
 6. Per method: every method moves by < 0.02 under (a). aucell/ora move least (rank-based over
    a 772-gene set); mlm is multivariate so it is the one that could move more.
 7. SHUF floor. Their own "random" column is 0.4939. A degree- and sign-preserving shuffle that
    KEEPS self-edges should land well above 0.5 for TP53 -- I predict 0.60-0.85 -- because a
    980-target random regulon still reads global cell state. So 0.957 is above the floor, but
    the floor is not 0.5 (T-038: do not quote nominal chance without checking).

WHAT WOULD REFUTE THE "no transfer" READING: a drop of more than 0.05 AUROC under (a), with
(b) agreeing. That would mean their headline metric sees the same leak our rank-among-TFs
metric sees on knockTF, and would upgrade F-023's "candidate abundance leak" from mechanism
argument to measurement on a published benchmark.

Run: cd analysis && uv run python bench2026.py    (SLOW: their dcp.py runs waggr/gsea with
1000 permutations each; the two full 11-method runs take ~30-60 min each on 48 cores. The
12 network-side runs use a patched copy restricted to five fast methods. Everything is cached
on the presence of the act csv, so a rerun is cheap.)
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/bench2026.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import difflib, os, shutil, subprocess, sys, time
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score, average_precision_score

REPO = os.path.abspath("../data/bench2026")
HOME = os.path.abspath("../data/bench2026_home")          # fake $HOME so their hardcoded ~/TFact resolves
BACKUP = os.path.abspath("../data/bench2026_backup")
PY = sys.executable
TYPE, TF = "Inhibition", "TP53"
# their dcp.py loops over 11 methods + consensus; mdt and udt need xgboost, which is not in this
# project's env and which nobody authorised adding, so their dcp.py prints "mdt error: xgboost is not
# installed" and skips them. 10 of their 12 output files are produced. Reported as a gap, not patched.
FULL = ["waggr", "ulm", "aucell", "viper", "gsea", "gsva", "mlm", "ora", "zscore", "consensus"]
FAST = ["ulm", "zscore", "aucell", "mlm", "ora"]
ENV = dict(os.environ, HOME=HOME)

# ---------------------------------------------------------------- 0. provenance
print("=" * 100)
head = subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
print(f"repo {REPO}\ncommit {head}")
os.makedirs(f"{REPO}/result/decoupler/{TYPE}", exist_ok=True)   # not shipped; their dcp.py assumes it
os.makedirs(f"{REPO}/result/auc/rank", exist_ok=True)           # not shipped; their auc.py line 222 assumes it
os.makedirs(HOME, exist_ok=True)
if not os.path.islink(f"{HOME}/TFact"):
    os.symlink(REPO, f"{HOME}/TFact")
print(f"$HOME for their scripts -> {HOME} (TFact -> repo). Their `~/TFact/...` paths are NOT patched.")

AUCS = ["auroc_cell", "auprc_cell", "au_rank", "au_coverages", "au_epr"]
if not os.path.isdir(BACKUP):           # their shipped result tables, kept pristine; auc.py appends to the live ones
    os.makedirs(BACKUP)
    for f in AUCS:
        shutil.copy(f"{REPO}/result/auc/{f}.csv", f"{BACKUP}/{f}.csv")
published = {f: pd.read_csv(f"{BACKUP}/{f}.csv", index_col=0) for f in AUCS}

# ---------------------------------------------------------------- 1. the one code patch, printed in full
orig = open(f"{REPO}/methord/decoupler/dcp.py").read()
patched = (orig
           .replace("""else:
    print("Net not supported!")""",
                    """else:
    netfile=os.path.expanduser(f"~/TFact/dataset/net_{net0}.csv")""")
           .replace("""for method_name in ['ulm','aucell','viper','gsea','gsva','mdt','mlm','ora','udt','zscore']:""",
                    """for method_name in os.environ['TFACT_METHODS'].split(','):""")
           .replace("""for p_name in ['waggr','ulm','aucell','viper','gsea','gsva','mdt','mlm','ora','udt','zscore']:""",
                    """for p_name in os.environ['TFACT_METHODS'].split(','):""")
           .replace("""try:
    dc.mt.waggr(adata, Net, tmin=3)
except :
    dc.mt.waggr(adata, Net, tmin=3)""", """if 'waggr' in os.environ['TFACT_METHODS']:
    dc.mt.waggr(adata, Net, tmin=3)""")
           .replace("""dc.mt.consensus(adata)""", """if len(os.environ['TFACT_METHODS'].split(',')) > 1: dc.mt.consensus(adata)""")
           .replace("""for i in ['ulm','aucell','viper','gsea','gsva','mdt','mlm','ora','udt','waggr','zscore','consensus']:""",
                    """for i in os.environ['TFACT_METHODS'].split(','):"""))
assert patched != orig
open(f"{REPO}/methord/decoupler/dcp_cond.py", "w").write(patched)
print("\n--- PATCH 1/1: methord/decoupler/dcp_cond.py (a copy; dcp.py itself is untouched and is what")
print("    the anchor and the masked condition run). It only (i) lets --net name an arbitrary")
print("    dataset/net_<name>.csv and (ii) reads the method list from $TFACT_METHODS so the 12")
print("    network-side runs can skip waggr/gsea/gsva/viper/mdt/udt. No scoring code is changed.")
for line in difflib.unified_diff(orig.splitlines(), patched.splitlines(), "dcp.py", "dcp_cond.py", lineterm="", n=1):
    print("   ", line)
print("""
--- PATCH 0/1: result/auc/auc.py is run UNPATCHED. Its `auprc` used-before-assignment bug is on
    line 86, inside `if Type == 'Activation'`; Dox is an Inhibition dataset, so that branch is
    never taken and the bug is not on this path. dataset/filter.py (the `.asType` bug) is not on
    this path either: the repo ships the filtered.h5ad.
    Noted, not patched, because nothing here depends on it: auc.py line 222 dumps
    result/auc/rank/<name>_<method>.csv AFTER the `if(True)` random block has rebound `df_rank`,
    so that per-cell file holds RANDOM ranks for every method. The headline au_rank.csv (line 129)
    is computed before the rebinding and is correct. This script therefore recomputes the per-cell
    rank fractions with a verbatim copy of their lines 97-129 instead of reading that dump.
""")

# ---------------------------------------------------------------- 2. condition inputs
net = pd.read_csv(f"{REPO}/dataset/net_collectri_filter.csv")
print(f"net_collectri_filter.csv: {len(net)} edges, {net.source.nunique()} TFs, "
      f"{(net.source == net.target).sum()} self-edges; {TF} regulon {(net.source == TF).sum()} edges; "
      f"{TF} is a target of {net[(net.target == TF) & (net.source != TF)].source.nunique()} other TFs")

def write_net(name, df):
    p = f"{REPO}/dataset/net_{name}.csv"
    if not os.path.exists(p):
        df.to_csv(p, index=False)
    return name

write_net("collectri_noself", net[~((net.source == TF) & (net.target == TF))])
write_net("collectri_noother", net[~((net.target == TF) & (net.source != TF))])

def shuffle_net(net, seed, mult=10):
    """Swap targets between non-self edges. Source and weight stay on their row, so every TF keeps
    its exact target count AND its exact activator/repressor mix; every target keeps its in-degree.
    Self-edges are held out of the swap, so autoregulation survives: this null is 'random regulon,
    real self-edge'. (controls.py shuffles an adjacency with a zero diagonal; same swap, edge-list form.)"""
    rng = np.random.default_rng(seed)
    keep = net[net.source == net.target]
    sub = net[net.source != net.target].reset_index(drop=True)
    src, tgt = sub.source.values, sub.target.values.copy()
    live = set(zip(src, tgt))
    for _ in range(mult * len(sub)):
        i, j = rng.integers(len(sub), size=2)
        a, b, c, d = src[i], tgt[i], src[j], tgt[j]
        if a == d or c == b or (a, d) in live or (c, b) in live:
            continue
        live.difference_update({(a, b), (c, d)}); live.update({(a, d), (c, b)})
        tgt[i], tgt[j] = d, b
    sub = sub.assign(target=tgt)
    return pd.concat([sub, keep], ignore_index=True)

SHUF = [f"collectri_shuf{s}" for s in range(10)]
for s, nm in enumerate(SHUF):
    if not os.path.exists(f"{REPO}/dataset/net_{nm}.csv"):
        t0 = time.time(); write_net(nm, shuffle_net(net, s)); print(f"  built net_{nm}.csv ({time.time()-t0:.0f}s)")
shuf0 = pd.read_csv(f"{REPO}/dataset/net_{SHUF[0]}.csv")
deg_ok = (shuf0.groupby("source").target.count().sort_index().equals(net.groupby("source").target.count().sort_index())
          and shuf0.groupby("source").weight.sum().sort_index().round(6).equals(net.groupby("source").weight.sum().sort_index().round(6))
          and shuf0.groupby("target").source.count().sort_index().equals(net.groupby("target").source.count().sort_index()))
surv = len(set(zip(shuf0.source, shuf0.target)) & set(zip(net.source, net.target))) / len(net)
print(f"  shuffle check: degree+sign-mix preserved = {deg_ok}; {surv:.1%} of original edges survive seed 0 "
      f"({(net.source == net.target).sum()} of them are the held-out self-edges = {(net.source==net.target).sum()/len(net):.1%})")

# masked h5ad: TP53's own gene zeroed in the counts their dcp.py reads
MASK = "DoxMASK"
os.makedirs(f"{REPO}/dataset/{TYPE}/{MASK}", exist_ok=True)
if not os.path.exists(f"{REPO}/dataset/{TYPE}/{MASK}/filtered.h5ad"):
    a = ad.read_h5ad(f"{REPO}/dataset/{TYPE}/Dox/filtered.h5ad")
    j = a.var_names.get_loc(TF)
    X = a.X.tolil() if sp.issparse(a.X) else a.X
    X[:, j] = 0
    a.X = X.tocsr() if sp.issparse(a.X) else X
    a.write_h5ad(f"{REPO}/dataset/{TYPE}/{MASK}/filtered.h5ad")
shutil.copy(f"{REPO}/result/express/{TYPE}/Dox_obs.csv", f"{REPO}/result/express/{TYPE}/{MASK}_obs.csv")

# control for the ORA anomaly found on the first run (+0.19 AUROC under masking, see the table):
# zero ONE RANDOM gene of the same detection rate instead of TP53. If ORA moves as much, the movement
# is their ORA's ordinal tie-breaking on a sparse matrix (decoupler _ora.py line 275 ranks with
# method='ordinal', so deleting any column reshuffles the top-775 cut), not the perturbed TF's mRNA.
RAND = "DoxRAND"
os.makedirs(f"{REPO}/dataset/{TYPE}/{RAND}", exist_ok=True)
if not os.path.exists(f"{REPO}/dataset/{TYPE}/{RAND}/filtered.h5ad"):
    a = ad.read_h5ad(f"{REPO}/dataset/{TYPE}/Dox/filtered.h5ad")
    det = np.asarray((a.X > 0).sum(0)).ravel()
    tgt = det[a.var_names.get_loc(TF)]
    cand = [i for i in np.where(np.abs(det - tgt) <= 0.2 * tgt)[0] if a.var_names[i] != TF]
    pick = int(np.random.default_rng(0).choice(cand))
    RANDGENE = a.var_names[pick]
    X = a.X.tolil() if sp.issparse(a.X) else a.X
    X[:, pick] = 0
    a.X = X.tocsr() if sp.issparse(a.X) else X
    a.uns["randgene"] = RANDGENE
    a.write_h5ad(f"{REPO}/dataset/{TYPE}/{RAND}/filtered.h5ad")
RANDGENE = str(ad.read_h5ad(f"{REPO}/dataset/{TYPE}/{RAND}/filtered.h5ad").uns["randgene"])
shutil.copy(f"{REPO}/result/express/{TYPE}/Dox_obs.csv", f"{REPO}/result/express/{TYPE}/{RAND}_obs.csv")

def col(a, g):
    x = a[:, g].X
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

a0 = ad.read_h5ad(f"{REPO}/dataset/{TYPE}/Dox/filtered.h5ad")
a1 = ad.read_h5ad(f"{REPO}/dataset/{TYPE}/{MASK}/filtered.h5ad")
MASK_OK = (col(a1, TF) == 0).all() and (col(a0, TF) > 0).any()
print(f"  masked h5ad: {TF} counts now all zero = {MASK_OK}; it was nonzero in "
      f"{(col(a0, TF) > 0).sum()}/{a0.n_obs} cells; {int((a0.X != a1.X).sum())} matrix entries differ")
obs = pd.read_csv(f"{REPO}/result/express/{TYPE}/Dox_obs.csv", index_col=0)
print(f"  cells: {(obs.gene == TF).sum()} {TF}, {(obs.gene == 'CTRL').sum()} CTRL")

# ---------------------------------------------------------------- 3. run their pipeline
def run_dcp(name, net0, methods, script="dcp.py"):
    want = [f"{REPO}/result/decoupler/{TYPE}/{name}_{m}_{net0}_act.csv" for m in methods]
    if all(os.path.exists(p) for p in want):
        return
    t0 = time.time()
    e = dict(ENV, TFACT_METHODS=",".join(methods))
    r = subprocess.run([PY, "-u", f"{REPO}/methord/decoupler/{script}", "--type", TYPE, "--name", name,
                        "--net", net0], env=e, cwd=REPO, capture_output=True, text=True)
    print(f"  dcp [{script} {name} {net0}] {time.time()-t0:6.0f}s rc={r.returncode} "
          f"{[l for l in r.stdout.splitlines() if 'error' in l or 'not found' in l][:3]}")

def run_auc(name, label):
    r = subprocess.run([PY, f"{REPO}/result/auc/auc.py", "--fold", "decoupler", "--method", label,
                        "--type", TYPE, "--name", name], env=ENV, cwd=REPO, capture_output=True, text=True)
    if r.returncode:
        print(f"  auc FAILED {name} {label}: {r.stderr.strip().splitlines()[-1]}")
    return r.returncode == 0

print("\nrunning their dcp.py (unpatched) for the anchor and the masked condition, and dcp_cond.py for the rest")
run_dcp("Dox", "collectri", FULL)
run_dcp(MASK, "collectri", FULL)
for net0 in ["collectri_noself", "collectri_noother"] + SHUF:
    run_dcp("Dox", net0, FAST, script="dcp_cond.py")
run_dcp(RAND, "collectri", FAST, script="dcp_cond.py")

# ---------------------------------------------------------------- 4. their metrics, via their auc.py
CONDS = [("anchor", "Dox", "collectri"), ("(a) MASK own gene", MASK, "collectri"),
         ("(b) NOSELF edge", "Dox", "collectri_noself"), ("(c) NOOTHER regulons", "Dox", "collectri_noother")] \
        + [(f"SHUF seed {s}", "Dox", nm) for s, nm in enumerate(SHUF)]

for f in AUCS:                                   # start from their shipped tables, as auc.sh does
    shutil.copy(f"{BACKUP}/{f}.csv", f"{REPO}/result/auc/{f}.csv")
for _, name, net0 in CONDS + [("rand", RAND, "collectri")]:
    ms = FULL if (net0 == "collectri" and name != RAND) else FAST
    for m in ms:
        if os.path.exists(f"{REPO}/result/decoupler/{TYPE}/{name}_{m}_{net0}_act.csv"):
            run_auc(name, f"{m}_{net0}")
res = {f: pd.read_csv(f"{REPO}/result/auc/{f}.csv", index_col=0) for f in AUCS}
for f in AUCS:                                   # restore their shipped tables
    shutil.copy(f"{BACKUP}/{f}.csv", f"{REPO}/result/auc/{f}.csv")

def their(f, name, m, net0):
    try:
        return float(res[f].loc[f"{TF}__{name}", f"{m}_{net0}"])
    except KeyError:
        return np.nan

# ---------------------------------------------------------------- 5. ANCHOR
print("\n" + "=" * 100)
print("ANCHOR: their dcp.py + their auc.py, unpatched, on their shipped Dox/filtered.h5ad")
print(f"{'metric':10s} {'ours':>10s} {'shipped (their result/auc/*.csv)':>34s}")
for f, pubcol in [("auroc_cell", "ULM_CollecTRI"), ("auprc_cell", "ULM_CollecTRI"),
                  ("au_rank", "ULM_CollecTRI"), ("au_epr", "ULM_CollecTRI")]:
    print(f"{f:10s} {their(f, 'Dox', 'ulm', 'collectri'):10.4f} {published[f].loc[f'{TF}__Dox', pubcol]:34.4f}")
print(f"{'auroc':10s} {their('auroc_cell', 'Dox', 'consensus', 'collectri'):10.4f} "
      f"{published['auroc_cell'].loc[f'{TF}__Dox', 'Consensus_CollecTRI']:34.4f}   (consensus)")
print(f"their own 'random' AUROC column for this dataset: {published['auroc_cell'].loc[f'{TF}__Dox','random']:.4f} "
      f"(rank {published['au_rank'].loc[f'{TF}__Dox','random']:.4f})")

# ---------------------------------------------------------------- 6. per-cell rank + bootstrap (their formulae)
def act_of(name, m, net0):
    p = f"{REPO}/result/decoupler/{TYPE}/{name}_{m}_{net0}_act.csv"
    if not os.path.exists(p):
        return None
    act = pd.read_csv(p, index_col=0).dropna()
    o = obs[obs.index.isin(act.columns)]
    return act.loc[:, o.index]

pert_cells = obs.index[obs.gene == TF]
ctrl_cells = obs.index[obs.gene == "CTRL"]

def auroc_all(act, align=None):
    """auc.py lines 76-93, Inhibition branch, applied to EVERY TF row (their code only loops over `pert`).
    align: reindex onto the anchor's TF list. Deleting TP53 from the other regulons can push a small
    regulon below their tmin=3, so a few TFs disappear from that condition; they become NaN here and
    are excluded pairwise (count printed), never silently compared against a different TF."""
    y = np.r_[np.ones(len(pert_cells)), np.zeros(len(ctrl_cells))]
    P = -act.loc[:, list(pert_cells) + list(ctrl_cells)].values
    s = pd.Series([roc_auc_score(y, p) for p in P], index=act.index)
    return s.reindex(align) if align is not None else s

def rankfrac(act):
    """auc.py lines 97-129, Inhibition branch: control-z-scored activity, TF ranked ascending per cell."""
    mu, sd = act[ctrl_cells].mean(axis=1), act[ctrl_cells].std(axis=1).replace(0, np.nan)
    z = act.sub(mu, axis=0).div(sd, axis=0).loc[:, pert_cells]
    out = {}
    for c in pert_cells:
        p = z[c].dropna()
        if TF in p.index:
            out[c] = rankdata(p.values, method="average")[p.index.get_loc(TF)] / act.shape[0]
    return pd.Series(out)

def boot(act, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    cells = list(pert_cells) + list(ctrl_cells)
    y = np.r_[np.ones(len(pert_cells)), np.zeros(len(ctrl_cells))]
    p = -act.loc[TF, cells].values
    out = [roc_auc_score(y[i], p[i]) for i in (rng.integers(len(cells), size=(n, len(cells))))
           if 0 < y[i].sum() < len(cells)]
    return np.percentile(out, [2.5, 97.5])

net_reg = net.groupby("source").target.apply(set)
self_tfs = set(net.loc[net.source == net.target, "source"])
holds_tf = {s for s in net_reg.index if TF in net_reg[s] and s != TF}
base_act = act_of("Dox", "ulm", "collectri")
scored = list(base_act.index)
null_tfs = [t for t in scored if t not in self_tfs and t not in holds_tf and t != TF]
print(f"\nscored TFs {len(scored)}; with a self-edge {len([t for t in scored if t in self_tfs])}; "
      f"holding {TF} in their regulon {len([t for t in scored if t in holds_tf])}; "
      f"internal null group (neither) {len(null_tfs)}")

base_auroc, base_rf = auroc_all(base_act), rankfrac(base_act)
rows = []
for lab, name, net0 in CONDS:
    act = act_of(name, "ulm", net0)
    if act is None:
        continue
    A, rf = auroc_all(act, scored), rankfrac(act)
    lo, hi = boot(act)
    hold_list = [t for t in scored if t in holds_tf]
    rows.append(dict(variant=lab, auroc=A[TF], ci=f"[{lo:.3f},{hi:.3f}]",
                     auprc=their("auprc_cell", name, "ulm", net0), rank=their("au_rank", name, "ulm", net0),
                     top10=float((rf <= .10).mean()), top5=float((rf <= .05).mean()),
                     better=float((rf < base_rf.reindex(rf.index)).mean()),
                     worse=float((rf > base_rf.reindex(rf.index)).mean()),
                     nTF=int(A.notna().sum()),
                     d_holds=float((A[hold_list] - base_auroc[hold_list]).abs().mean()),
                     d_null=float((A[null_tfs] - base_auroc[null_tfs]).abs().mean())))
R = pd.DataFrame(rows).set_index("variant")
shuf = R[R.index.str.startswith("SHUF")]
R = pd.concat([R[~R.index.str.startswith("SHUF")],
               pd.DataFrame([dict(auroc=shuf.auroc.mean(), ci=f"[{shuf.auroc.min():.3f},{shuf.auroc.max():.3f}]",
                                  auprc=shuf.auprc.mean(), rank=shuf['rank'].mean(), top10=shuf.top10.mean(),
                                  top5=shuf.top5.mean(), better=shuf.better.mean(), worse=shuf.worse.mean(),
                                  nTF=shuf.nTF.mean(), d_holds=shuf.d_holds.mean(), d_null=shuf.d_null.mean())],
                            index=["SHUF x10 mean [min-max]"])])

print("\n" + "=" * 100)
print("ULM + CollecTRI on their pipeline. auroc/auprc/rank are THEIR metrics from THEIR auc.py.")
print("ci = 2.5-97.5 percentile of 2000 bootstrap resamples of the 617 CELLS (seed 0); the unit is")
print("cells, not experiments: this is n=1 dataset, n=1 perturbed TF. top10/top5/better/worse are over")
print("the 197 perturbed cells, using their per-cell rank of TP53 among the scored TFs; better/worse")
print("are paired against the anchor (unmasked ULM) cell by cell.")
print("d_holds / d_null: mean |dAUROC| vs anchor over the TFs holding TP53 in their regulon, and over")
print("the internal null group (no self-edge, TP53 not a target).")
print(R.to_string(float_format=lambda v: f"{v:.4f}"))

# ---------------------------------------------------------------- 7. all 11 methods, unmasked vs masked
print("\n" + "=" * 100)
print("All methods their dcp.py loops over: THEIR AUROC, anchor vs (a) TP53's own gene masked")
print(f"{'method':11s} {'anchor':>8s} {'masked':>8s} {'delta':>8s} {'shipped':>8s} | {'rank':>7s} {'rank_m':>7s}")
mrows = []
for m in FULL:
    u, k = their("auroc_cell", "Dox", m, "collectri"), their("auroc_cell", MASK, m, "collectri")
    pubname = {"ulm": "ULM", "aucell": "AUCell", "viper": "VIPER", "gsea": "GSEA", "gsva": "GSVA",
               "mlm": "MLM", "ora": "ORA", "waggr": "WAGGR", "zscore": "Zscore",
               "consensus": "Consensus"}[m] + "_CollecTRI"
    sh = published["auroc_cell"].loc[f"{TF}__Dox", pubname]
    ru, rk = their("au_rank", "Dox", m, "collectri"), their("au_rank", MASK, m, "collectri")
    print(f"{m:11s} {u:8.4f} {k:8.4f} {k-u:+8.4f} {sh:8.4f} | {ru:7.4f} {rk:7.4f}")
    mrows.append((m, u, k, k - u))
M = pd.DataFrame(mrows, columns=["method", "anchor", "masked", "delta"]).set_index("method")
print(f"\nControl for the above: one RANDOM gene of the same detection rate as {TF} zeroed instead "
      f"({RANDGENE}, detected in the same number of cells +/-20%).")
print(f"{'method':11s} {'anchor':>8s} {'randmask':>9s} {'d_rand':>8s} {'d_TP53':>8s}")
for m in FAST:
    u, rd = their("auroc_cell", "Dox", m, "collectri"), their("auroc_cell", RAND, m, "collectri")
    print(f"{m:11s} {u:8.4f} {rd:9.4f} {rd-u:+8.4f} {M.loc[m,'delta']:+8.4f}")

# ---------------------------------------------------------------- 8. internal control by TF group
print("\n" + "=" * 100)
print("INTERNAL CONTROL (ULM). Mean AUROC by TF group, and mean |dAUROC| vs anchor.")
print("The null group can only move through the library-size change that masking causes; under the")
print("two network conditions it cannot move at all, which is what the assert checks.")
groups = {f"{TF} itself": [TF],
          f"holds {TF} in regulon (n={len([t for t in scored if t in holds_tf])})": [t for t in scored if t in holds_tf],
          f"internal null: no self-edge, no {TF} (n={len(null_tfs)})": null_tfs,
          f"has self-edge, no {TF} (n={len([t for t in scored if t in self_tfs and t not in holds_tf and t != TF])})":
              [t for t in scored if t in self_tfs and t not in holds_tf and t != TF]}
print(f"{'group':52s} {'anchor':>8s} " + " ".join(f"{l[:14]:>15s}" for l, _, _ in CONDS[1:4]))
for g, ts in groups.items():
    line = f"{g:52s} {base_auroc[ts].mean():8.4f} "
    for lab, name, net0 in CONDS[1:4]:
        line += f" {np.abs(auroc_all(act_of(name, 'ulm', net0), scored)[ts] - base_auroc[ts]).mean():15.3e}"
    print(line)

print("\n" + "=" * 100)
d_mask = R.loc["(a) MASK own gene", "auroc"] - R.loc["anchor", "auroc"]
d_self = R.loc["(b) NOSELF edge", "auroc"] - R.loc["anchor", "auroc"]
print(f"TRANSFER TEST: anchor {R.loc['anchor','auroc']:.4f}; masking TP53 moves their headline AUROC by "
      f"{d_mask:+.4f}; deleting only the self-edge by {d_self:+.4f}.")
print(f"Refutation threshold stated in the docstring was |delta| > 0.05 for 'the leak transfers'. "
      f"{'TRANSFERS' if abs(d_mask) > 0.05 else 'DOES NOT TRANSFER at that threshold'}.")

# F-029: persist. The per-variant AUROC table and the per-method anchor-vs-masked table.
# NOT YET GENERATED: this script takes about 43 minutes and was not rerun when F-029 was filed.
R.to_csv("../results/bench2026.csv")
M.to_csv("../results/bench2026_methods.csv")

noself_A = auroc_all(act_of("Dox", "ulm", "collectri_noself"), scored)
noother_A = auroc_all(act_of("Dox", "ulm", "collectri_noother"), scored)
assert MASK_OK and deg_ok \
    and np.nanmax(np.abs(noself_A[null_tfs] - base_auroc[null_tfs])) < 1e-9 \
    and np.nanmax(np.abs(noother_A[null_tfs] - base_auroc[null_tfs])) < 1e-9 \
    and abs(noother_A[TF] - base_auroc[TF]) < 1e-9 \
    and abs(their("auroc_cell", "Dox", "ulm", "collectri") - published["auroc_cell"].loc[f"{TF}__Dox", "ULM_CollecTRI"]) < 0.02, \
    "logic broken: masked matrix / degree-preserving shuffle / univariate null group / anchor reproduction"
print("\nassert passed: TP53 zeroed in the masked matrix, shuffle preserves degree and sign mix, the "
      "internal null group and TP53's own column are untouched by the network edits that cannot reach "
      "them, and the anchor reproduces their shipped ULM_CollecTRI AUROC to within 0.02.")
