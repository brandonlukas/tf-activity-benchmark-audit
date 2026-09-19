"""Does masked ULM on CollecTRI track TF activity that the TF's own mRNA cannot explain?

Ground truth with NO perturbed TF: SIGNOR directional regulatory phosphosites read out of the
MPXV-infected primary human fibroblast time course (Nat Commun s41467-024-51074-6, Suppl.
Dataset 2 = transcriptome log2FC, Dataset 4 = phosphoproteome log2FC WITH a protein_fold_change
column per timepoint). Because nothing is perturbed, the masked-gene control that dominates this
project is STRUCTURALLY UNNECESSARY here -- there is no own-mRNA drop to leak. It is run anyway
(and via F-010's self-edge equivalent) only to show it changes nothing.

PREDICTIONS, written before running:
  pooled Spearman(masked ULM, protein-corrected phospho activity) = +0.25, CI excluding 0
  degree-preserving shuffled network        0.00 +/- 0.10
  placebo (TF labels permuted on phospho)   0.00 +/- 0.10
  own mRNA log2FC vs phospho activity       |r| < 0.15, CI NOT excluding moderate values (small n)
WHAT WOULD KILL THE DIRECTION: masked ULM's correlation CI containing 0 at this n -- CollecTRI
regulons would carry no information about activity in a setting where abundance cannot explain it.
n is about 23 TFs x 3 timepoints, so "inconclusive at n=69" is a likely and acceptable outcome.

Retrieval analogue of the knockTF hit rate (there is no perturbed TF, so top-10%/top-5% needs a
definition): exactly leak.py's evaluate(), with sign(phospho activity) in the role of `type_p`.
For each TF-timepoint observation, rank all scored TFs by sign(phospho activity) * ULM score at
that timepoint and take the true TF's rank fraction. Chance is 0.10 / 0.05 by construction.
"""
# F-029: persist, do not only print (CLAUDE.md conventions). The tee is installed only when
# this file is the one being run, so an importing script's log stays complete and this one's
# log stays its own; it changes no printed line.
if __name__ == "__main__":
    import sys, pathlib
    pathlib.Path("../results").mkdir(exist_ok=True)
    _log = open("../results/phospho.log", "w", buffering=1)
    class _Tee:
        def write(self, t): sys.__stdout__.write(t); _log.write(t)
        def flush(self): sys.__stdout__.flush(); _log.flush()
    sys.stdout = _Tee()

import csv, re, zipfile
import anndata as ad, decoupler as dc, numpy as np, pandas as pd
from scipy.stats import rankdata
from xml.etree import ElementTree as ET

NB, NPERM, NSEED = 2000, 500, 10
TPS = ["6h", "12h", "24h"]
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

# ---------------------------------------------------------------- stdlib xlsx (openpyxl absent)
def read_sheet(path, sheet_name):
    """Minimal xlsx reader -> DataFrame of strings. No new dependency: openpyxl is not installed
    and pandas cannot open xlsx without it, so the sheet is parsed out of the zip with stdlib."""
    z = zipfile.ZipFile(path)
    rid = dict(re.findall(r'<sheet name="([^"]+)"[^>]*r:id="([^"]+)"', z.read("xl/workbook.xml").decode()))[sheet_name]
    rels = z.read("xl/_rels/workbook.xml.rels").decode()
    tgt = dict(re.findall(r'Id="([^"]+)" Type="[^"]*" Target="([^"]+)"', rels))[rid]
    shared = ["".join(t.text or "" for t in si.iter(NS + "t"))
              for si in ET.fromstring(z.read("xl/sharedStrings.xml"))]
    rows = []
    for _, row in ET.iterparse(z.open("xl/" + tgt.lstrip("/"))):
        if row.tag != NS + "row": continue
        cells = {}
        for c in row.findall(NS + "c"):
            j = 0
            for ch in re.match(r"[A-Z]+", c.get("r")).group(): j = j * 26 + ord(ch) - 64
            v = c.find(NS + "v")
            if c.get("t") == "inlineStr": cells[j - 1] = "".join(t.text or "" for t in c.iter(NS + "t"))
            elif v is not None: cells[j - 1] = shared[int(v.text)] if c.get("t") == "s" else v.text
        if cells: rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
        row.clear()
    w = max(len(r) for r in rows)
    rows = [r + [""] * (w - len(r)) for r in rows]
    return pd.DataFrame(rows[1:], columns=rows[0])

num = lambda s: pd.to_numeric(s.replace({"NA": None, "": None, "NaN": None}), errors="coerce")

# ---------------------------------------------------------------- SIGNOR: directional ACTIVITY sites
SCOLS = ["ENTITYA", "TYPEA", "IDA", "DATABASEA", "ENTITYB", "TYPEB", "IDB", "DATABASEB", "EFFECT",
         "MECHANISM", "RESIDUE", "SEQUENCE", "TAX_ID", "CELL_DATA", "TISSUE_DATA", "MODULATOR_COMPLEX",
         "TARGET_COMPLEX", "MODIFICATIONA", "MODASEQ", "MODIFICATIONB", "MODBSEQ", "PMID", "DIRECT",
         "NOTES", "ANNOTATOR", "SENTENCE", "SIGNOR_ID", "SCORE", "X"]  # the file ships without a header
AA1 = {"Ser": "S", "Thr": "T", "Tyr": "Y", "His": "H", "Asp": "D", "Cys": "C", "Lys": "K", "Arg": "R"}
sig_rows = list(csv.reader(open("../data/signor_human.tsv", encoding="utf-8", newline=""),
                           delimiter="\t", quoting=csv.QUOTE_NONE))
sig = pd.DataFrame(sig_rows, columns=SCOLS)
# V-023: never "regulates quantity" -- their showcase site YBX1 S174 is a QUANTITY annotation,
# i.e. an abundance readout, which is exactly the confound this project exists to avoid.
sig = sig[(sig.MECHANISM == "phosphorylation") & (sig.TAX_ID == "9606") & (sig.TYPEB == "protein")
          & sig.EFFECT.isin(["up-regulates activity", "down-regulates activity"])]
site_sign = {}
for gene, res, eff in zip(sig.ENTITYB, sig.RESIDUE, sig.EFFECT):
    for r in str(res).split(";"):
        m = re.fullmatch(r"\s*([A-Za-z]{3})\s*(\d+)\s*", r)
        if not m or m.group(1) not in AA1: continue
        site_sign.setdefault((gene, AA1[m.group(1)], int(m.group(2))), []).append(
            +1 if eff.startswith("up") else -1)
n_raw = len(site_sign)
conflict = [k for k, v in site_sign.items() if len(set(v)) > 1]
site_sign = {k: v[0] for k, v in site_sign.items() if len(set(v)) == 1}  # drop up/down conflicts
print(f"SIGNOR: {len(sig_rows)} rows -> {len(sig)} directional-activity phosphorylation rows -> "
      f"{n_raw} sites on {len({k[0] for k in site_sign})} proteins "
      f"({len(conflict)} sites dropped for conflicting up/down annotations)")

# ---------------------------------------------------------------- phosphoproteome
ph = read_sheet("../data/mpxv_phospho.xlsx", "A")
ph = ph[(ph.is_viral == "0") & (ph.is_contaminant == "0")]
ph["pos"] = num(ph.ptm_pos).astype("Int64")
key = list(zip(ph.gene_name, ph.ptm_AA, ph.pos))
ph["sign"] = [site_sign.get(k, 0) for k in key]
for t in TPS:
    ph["site." + t] = num(ph["fold_change_log2." + t])
    ph["prot." + t] = num(ph["protein_fold_change_log2." + t])
    ph["corr." + t] = ph["site." + t] - ph["prot." + t]      # protein-abundance corrected
dir_sites = ph[ph.sign != 0]
print(f"phosphoproteome: {len(ph)} non-viral non-contaminant site rows, "
      f"{len(dir_sites)} carry a directional SIGNOR activity annotation "
      f"({dir_sites.gene_name.nunique()} proteins)")

# ---------------------------------------------------------------- transcriptome + CollecTRI ULM
tr = read_sheet("../data/mpxv_transcriptome.xlsx", "A")
tr = tr.assign(**{t: num(tr["fold_change_log2." + t]) for t in TPS})
tr = tr.groupby("gene_name")[TPS].mean().dropna()
X = tr.T                                                       # 3 timepoints x genes
genes = np.array(X.columns)
net = pd.read_parquet("../data/collectri.parquet")
pnet = dc.pp.prune(features=genes, net=net, tmin=5)
sources, feats, adjm = dc.pp.adjmat(features=genes, net=pnet)
src_i = {s: i for i, s in enumerate(sources)}
gene_i = {g: i for i, g in enumerate(genes)}
mat = X.values.astype(np.float64)
score = lambda mm, adj: dc.mt.ulm.func(mm, adj)[0]             # decoupler's own ULM kernel

# panel = CollecTRI TFs (tmin=5) with a directional site quantified here at >=1 timepoint
panel = sorted({g for g in dir_sites.gene_name.unique() if g in src_i})
print(f"transcriptome: {len(genes)} genes, CollecTRI tmin=5 -> {len(sources)} scoreable TFs; "
      f"panel = {len(panel)} TFs with a directional site, "
      f"{len(dir_sites[dir_sites.gene_name.isin(panel)])} site rows")
print(f"   regulon sizes: " + " ".join(str((adjm[:, src_i[t]] != 0).sum()) for t in panel))

# phospho-derived activity: MEAN over the TF's directional sites of sign * site log2FC
def phospho_activity(col):
    d = dir_sites.dropna(subset=[col])
    v = (d["sign"] * d[col]).groupby(d.gene_name).mean()
    return v.reindex(panel)
PA = {"corrected": pd.DataFrame({t: phospho_activity("corr." + t) for t in TPS}),
      "uncorrected": pd.DataFrame({t: phospho_activity("site." + t) for t in TPS})}
nsite = dir_sites[dir_sites.gene_name.isin(panel)].groupby("gene_name").size().reindex(panel)

# ULM variants ------------------------------------------------------------------------------
es_ulm = pd.DataFrame(score(mat, adjm), index=TPS, columns=sources)
adjm_ns = adjm.copy()                                          # F-010: self-edge deletion == masking
for s in sources:
    if s in gene_i: adjm_ns[gene_i[s], src_i[s]] = 0
es_selfedge = pd.DataFrame(score(mat, adjm_ns), index=TPS, columns=sources)
es_masked = es_ulm.copy()                                      # leak.py-style: zero the TF's own gene
for t in panel:
    mm = mat.copy(); mm[:, gene_i[t]] = 0.0
    es_masked[t] = score(mm, adjm)[:, src_i[t]]

def shuffle_tf_gene(adj, seed):
    """fulldose.py's degree- and sign-mix-preserving null: each TF keeps its number of +1 and -1
    targets, redrawn at random from the measured genes other than its own."""
    r = np.random.default_rng(seed)
    B = np.zeros_like(adj)
    for s, c in src_i.items():
        w = adj[:, c][adj[:, c] != 0]
        if len(w) == 0: continue
        pool = np.arange(adj.shape[0])
        if s in gene_i: pool = pool[pool != gene_i[s]]
        B[r.choice(pool, size=len(w), replace=False), c] = w
    return B
es_shuf = [pd.DataFrame(score(mat, shuffle_tf_gene(adjm, s)), index=TPS, columns=sources)
           for s in range(NSEED)]

# ---------------------------------------------------------------- pooled statistics
tf_of = np.array([t for t in panel for _ in TPS])
tp_of = np.array([k for _ in panel for k in TPS])
def pair(es, pa):
    x = np.array([es.loc[k, t] for t, k in zip(tf_of, tp_of)])
    y = np.array([pa.loc[t, k] for t, k in zip(tf_of, tp_of)])
    return x, y, np.isfinite(x) & np.isfinite(y)
pear = lambda x, y: float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 and x.std() and y.std() else np.nan
spear = lambda x, y: pear(rankdata(x), rankdata(y))

def boot_ci(x, y, ok, stat, seed=0, nb=NB):
    """TF-cluster bootstrap: resample the panel TFs with replacement, recompute the pooled stat."""
    r = np.random.default_rng(seed)
    idx = {t: np.flatnonzero((tf_of == t) & ok) for t in panel}
    pool = [t for t in panel if len(idx[t])]
    out = []
    for _ in range(nb):
        rows = np.concatenate([idx[t] for t in r.choice(pool, size=len(pool), replace=True)])
        out.append(stat(x[rows], y[rows]))
    out = np.array(out, float); out = out[np.isfinite(out)]
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

# retrieval analogue (leak.py's evaluate, sign(phospho activity) in place of type_p)
def rankfrac(es, pa):
    fr, keep = [], []
    for t, k in zip(tf_of, tp_of):
        a = pa.loc[t, k]
        if not np.isfinite(a) or a == 0: fr.append(np.nan); keep.append(False); continue
        s = np.sign(a) * es.loc[k]
        fr.append(s.rank(ascending=False)[t] / es.shape[1]); keep.append(True)
    return np.array(fr), np.array(keep)

def row(name, es, pa, base_fr=None, extra=""):
    x, y, ok = pair(es, pa)
    rho, r = spear(x[ok], y[ok]), pear(x[ok], y[ok])
    lo, hi = boot_ci(x, y, ok, spear)
    fr, keep = rankfrac(es, pa)
    t10, t5 = np.mean(fr[keep] <= .10), np.mean(fr[keep] <= .05)
    b = w = np.nan
    if base_fr is not None:
        m = keep & np.isfinite(base_fr)
        b, w = np.mean(fr[m] < base_fr[m]), np.mean(fr[m] > base_fr[m])
    print(f"{name:34s} {ok.sum():4d} {rho:+.3f} [{lo:+.2f},{hi:+.2f}] {r:+.3f}  "
          f"{t10:.2f}  {t5:.2f}  {b:.2f}  {w:.2f}  {extra}")
    return rho, (lo, hi), fr

for lab, pa in PA.items():
    print(f"\n=== phospho activity: {lab} (mean over the TF's directional sites of "
          f"sign * {'site - protein log2FC' if lab == 'corrected' else 'site log2FC'})")
    print(f"{'variant':34s} {'n':>4s} {'spearman [95% TF-boot]':22s} {'pears':6s} top10  top5  bettr  worse")
    _, _, fr_ulm = row("ULM (unmasked)", es_ulm, pa)
    rho_m, ci_m, _ = row("ULM masked (own gene zeroed)", es_masked, pa, fr_ulm)
    row("ULM self-edge deleted (F-010)", es_selfedge, pa, fr_ulm)
    rr = [row(f"  SHUFFLED net seed {s}", es_shuf[s], pa, fr_ulm)[0] for s in range(NSEED)] \
         if lab == "corrected" else []
    if rr:
        print(f"{'SHUFFLED net, 10 seeds: mean':34s} {np.nanmean(rr):+.3f} "
              f"[{np.nanmin(rr):+.3f},{np.nanmax(rr):+.3f}]")
        # placebo: each TF's ULM score paired with a DIFFERENT TF's phospho trajectory
        x, y, ok = pair(es_masked, pa)
        r0 = np.random.default_rng(0); pl = []
        for _ in range(NPERM):
            while True:
                p = r0.permutation(len(panel))
                if not (p == np.arange(len(panel))).any(): break
            yp = np.array([pa.iloc[p[i], TPS.index(k)] for i, t in enumerate(panel) for k in TPS])
            m = ok & np.isfinite(yp)
            pl.append(spear(x[m], yp[m]))
        pl = np.array(pl, float)
        print(f"{'PLACEBO ground truth, 500 perms':34s} {np.nanmean(pl):+.3f} "
              f"[{np.nanpercentile(pl, 2.5):+.2f},{np.nanpercentile(pl, 97.5):+.2f}]  "
              f"(real masked {rho_m:+.3f}, one-sided p = {np.mean(pl >= rho_m):.3f})")
        PLACEBO, SHUF, RHO_M, CI_M = pl, np.array(rr, float), rho_m, ci_m

# ---------------------------------------------------------------- the direction's own number
print(f"\n=== own mRNA vs phospho-derived activity (Sugimoto's -0.02; unit = TF-timepoint)")
own = pd.DataFrame({k: [tr[k].get(t, np.nan) for t in panel] for k in TPS}, index=panel)
for lab, pa in PA.items():
    x, y, ok = pair(own.T, pa)                       # own.T is timepoints x TFs, like an es frame
    rho, r = spear(x[ok], y[ok]), pear(x[ok], y[ok])
    lo, hi = boot_ci(x, y, ok, spear); lo2, hi2 = boot_ci(x, y, ok, pear)
    print(f"  own mRNA log2FC vs {lab:12s} n={ok.sum():3d}  spearman {rho:+.3f} [{lo:+.2f},{hi:+.2f}]"
          f"   pearson {r:+.3f} [{lo2:+.2f},{hi2:+.2f}]")

# the correction subtracts protein log2FC, and protein tracks mRNA, so it can INDUCE a negative
# own-mRNA correlation by construction; quantify that before reading the number above.
prot = pd.DataFrame({k: dir_sites[dir_sites.gene_name.isin(panel)].groupby("gene_name")["prot." + k]
                     .mean().reindex(panel) for k in TPS}, index=panel)
xp, yp, okp = pair(own.T, prot)
print(f"  own mRNA log2FC vs own PROTEIN log2FC   n={okp.sum():3d}  spearman {spear(xp[okp], yp[okp]):+.3f}"
      f"   pearson {pear(xp[okp], yp[okp]):+.3f}  <- a positive value here makes the corrected "
      f"own-mRNA correlation negative by construction")

print(f"\npanel TFs (sites quantified): " + ", ".join(f"{t}:{int(nsite[t])}" for t in panel))
hw = (CI_M[1] - CI_M[0]) / 2
print(f"\nmasked-ULM Spearman CI half-width {hw:.3f} at {len(panel)} TFs; CI width scales as "
      f"1/sqrt(n_TF), so +/-0.10 needs about {int(np.ceil(len(panel) * (hw / .10) ** 2))} TFs "
      f"(~{int(np.ceil(len(panel) * (hw / .10) ** 2 / len(panel)))}x this phosphoproteome's directional panel).")

# F-029: persist. The panel's ULM scores under all three network/mask conditions, timepoint x TF,
# which is what every correlation in this script is computed from. No printed line changes.
pd.concat({"ulm": es_ulm[panel], "masked": es_masked[panel], "no_self_edge": es_selfedge[panel]},
          names=["variant", "timepoint"]).to_csv("../results/phospho.csv")

self_panel = [t for t in panel if t in gene_i and adjm[gene_i[t], src_i[t]] != 0]
print(f"panel TFs with a CollecTRI self-edge at tmin=5: {len(self_panel)} of {len(panel)} "
      f"({', '.join(self_panel) if self_panel else 'none'}); "
      f"max |masked - self-edge-deleted| over the panel = "
      f"{np.abs(es_masked[panel].values - es_selfedge[panel].values).max():.3f}")
assert (len(panel) >= 10 and np.isfinite(RHO_M) and abs(np.nanmean(PLACEBO)) < 0.10
        and np.abs(es_masked[panel].values - es_selfedge[panel].values).max() < 0.5
        and all(abs(es_masked.loc[k, t] - es_ulm.loc[k, t]) > 1e-6 for t in self_panel for k in TPS)), \
    "panel empty, placebo not centred on 0, masking did not move a self-edge TF's score, " \
    "or masking disagrees with F-010's self-edge deletion"
