"""Gap check: where does ULM+CollecTRI miss on knockTF, and is it regulon overlap or size?"""
import anndata as ad, decoupler as dc, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

adata = ad.read_h5ad("../data/knocktf.h5ad")
net = pd.read_parquet("../data/collectri.parquet")
dc.mt.ulm(adata, net, tmin=5)
S = adata.obsm["score_ulm"]                      # experiments x TFs
tfs = S.columns

# regulon size + pairwise Jaccard overlap between regulons
reg = net[net.source.isin(tfs)].groupby("source").target.apply(set)
size = reg.map(len).reindex(tfs)
M = pd.DataFrame({t: {u: len(reg[t] & reg[u]) / len(reg[t] | reg[u]) for u in tfs} for t in tfs})
M = M.mask(np.eye(len(tfs), dtype=bool), 0)

rows = []
for exp, tf in adata.obs.source.items():
    if tf not in tfs: continue
    s = S.loc[exp]
    sign = adata.obs.loc[exp, "type_p"]            # -1 knockdown, +1 overexpression
    r = (sign * s).rank(ascending=False)           # 1 = true TF is the most extreme in the expected direction
    nb = M[tf].nlargest(5)                         # 5 most overlapping TFs
    fam = (sign * s[[tf, *nb.index]]).rank(ascending=False)
    rows.append(dict(exp=exp, tf=tf, rank=r[tf], frac=r[tf] / len(tfs), fam_rank=fam[tf],
                     size=size[tf], max_jac=nb.iloc[0], top_nb=nb.index[0]))
df = pd.DataFrame(rows)
df["miss"] = df.frac > 0.10                        # true TF not in top 10%
df["fam_miss"] = df.fam_rank > 1                   # beaten by an overlap neighbour
assert df["rank"].min() >= 1 and df.frac.max() <= 1
df.to_csv("gapcheck.csv", index=False)

print(f"{len(df)} experiments, {df.tf.nunique()} TFs, top-10% hit rate {1 - df.miss.mean():.2f}, "
      f"beaten by a neighbour {df.fam_miss.mean():.2f}")
for col in ["size", "max_jac"]:
    df[col + "_bin"] = pd.qcut(df[col], 3, labels=["low", "mid", "high"], duplicates="drop")
    print(f"\nby {col}:\n", df.groupby(col + "_bin", observed=True)[["miss", "fam_miss"]].mean().round(2)
          .join(df.groupby(col + "_bin", observed=True).size().rename("n")))
print("\nmost-missed TFs:\n", df.groupby("tf").agg(n=("miss", "size"), miss=("miss", "mean"), fam_miss=("fam_miss", "mean"),
      size=("size", "first"), jac=("max_jac", "first"), nb=("top_nb", "first")).query("n>=3").sort_values("miss", ascending=False).head(15).round(2))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].scatter(df.max_jac, df.frac, c=df.fam_miss.map({True: "C3", False: "C0"}), s=12, alpha=.6)
ax[0].set(xlabel="max Jaccard with another regulon", ylabel="true-TF rank (fraction)", title="red = beaten by overlap neighbour")
ax[1].scatter(df["size"], df.frac, s=12, alpha=.6); ax[1].set(xscale="log", xlabel="regulon size", ylabel="true-TF rank (fraction)")
fig.tight_layout(); fig.savefig("gapcheck.png", dpi=130)
