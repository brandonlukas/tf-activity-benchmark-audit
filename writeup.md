# What a transcription-factor activity benchmark actually measures

Draft, 2026-09-18. Each number traces to an `F-nnn` claim in `findings.md` and its `V-nnn` verdict; where a verdict split, the wording here is the verdict's, not the claim's. Claims still lacking a verdict are marked **unverified** in place.

**Denominators, because every rate depends on them.** knockTF as `decoupler` fetches it is 388 experiments — its loader defaults to `thr_fc=-1`, silently keeping only those whose perturbed TF's own mRNA fell below −1, out of 907 (F-018/V-018). Of those 388, **279** (155 TFs) get a ULM score at `tmin=5`; the other 109 are dropped because the TF is absent from CollecTRI (90) or has fewer than five measured targets (19). Unless stated otherwise, every rate below is over those 279. It is not a rate over knockTF (F-013/V-013). Two other sets appear where named: the **643** scored experiments of the unfiltered file, and the **106** that support the study-context contrast.

## Summary

Benchmarks for transcription-factor activity inference score a method by whether it recovers the TF that was perturbed. On the standard benchmark this rewards detecting the perturbed TF's **own transcript**. We quantify the channel that carries it and give the correction.

1. The perturbed TF's own log fold change **alone** ranks it in the top 10% in **98%** of the 279 experiments, against **45%** for ULM+CollecTRI (F-001, F-013/V-013).
2. For a univariate scorer the channel is **one edge**: the autoregulatory self-edge X→X. Deleting it reproduces masking the gene at top-10% — 96/279 either way — and accounts for 24% of ULM's hits (F-010/V-010).
3. Its width obeys a law. On the 227 TFs of an independent dataset that have a self-edge, a measured own gene and a surviving regulon, `log₁₀|Δscore| ≈ −0.53 − 0.46·log₁₀(size)`, Spearman −0.868, R² = 0.80 for the joint fit including a detection term. The −0.46 exponent is the 1/√n one gene contributes to a t-statistic. The leak is a **product** of that width and the own gene's discriminative signal (F-027/V-027).
4. Corrected, ULM+CollecTRI scores **34%** against a per-TF base rate of **12%** — real signal, about 2.8× chance, and about a quarter below the 45% an unmasked run on the same 279 experiments reports (F-001/V-001, F-005, V-022).

A second, independent confound: **study context**. On the 106 experiments that support the contrast, under a symmetric mask, a TF reaches the top 10% in 22% of its own GEO series' other experiments — when a *different* TF was perturbed — against 14% in same-cell-line experiments from other series (F-028/V-028).

## 1. The benchmark rewards abundance

Adding the TF's own log fold change to ULM's score lifts the hit rate from 45% to 66%, 77% and **94%** at three weights, better than ULM on 89–91% of experiments; the own log fold change with no regulon at all reaches **98%** (F-013/V-013). On this benchmark it beats ULM+CollecTRI outright. Whether it beats other regulon methods is untested here; an own-expression-only variant did beat eleven methods on a different benchmark with a different network (Yashar et al. 2024, iScience — F-025, **unverified**).

Two structural points compound it:

- **The loader filters on the same axis.** `dc.ds.knocktf()` keeps 388 of 907 experiments by the sign and size of the TF's own mRNA change (F-018/V-018).
- **A recent single-cell benchmark's preprocessing does the same**, keeping a TF only if its mean is lower in the perturbed arm (V-027).

There is also no measured zero dose. An exact `0.0` in the knockTF matrix is a fill value for "not measured": of the 907 rows, 82 of the 140 with non-negative own logFC are missing values; among the **643 scored**, the non-negative bin is 38 unmeasured, 21 inside the experiment's noise and only 18 measurably raised (F-022/V-022).

## 2. The mechanism: one edge

ULM regresses target log fold changes on a TF's regulon column, so the perturbed TF's own transcript enters its *own* score only through a self-edge X→X.

Deleting that edge — leaving the gene in the matrix — gives 96/279 hits, the same count as masking, with the two hit sets differing on 2 of 279. Removing the gene from every *other* TF's regulon leaves the score at 0.455, unchanged. The equivalence holds **at top-10%**; at top-5% the two conditions differ slightly (0.280 masked against 0.272), and for propagation by about six experiments (F-010/V-010).

The effect splits by whether the TF has a self-edge at all:

| | n | ULM | Self-edge removed |
|---|---|---|---|
| Has a self-edge | 197 | 0.538 | 0.386 |
| No self-edge | 82 | 0.244 | 0.244 |

Masking costs the no-self-edge group nothing — it has no channel to lose (F-010/V-010).

It also splits by the edge's **sign**. In the set as fetched every knockdown lowers the TF's mRNA, so a +1 self-edge moves the score the rewarded way and a −1 edge the other. Removing the edge worsens 81 of 84 TFs with a +1 edge and improves 11 of 12 with a −1 edge (F-011/V-011).

**The sign is usually not per-edge evidence.** CollecTRI ships a `sign_decision` column: 58% of its +1 edges, and **80% of its +1 self-edges**, are a default fill for "no information" — chosen on the base rate of activation among signed edges (71–76%) *and* because it benchmarked better than defaulting to repression (F-016/V-016). Of the 176 experiments whose TF has a +1 self-edge, 92 rest on that fill (F-017/V-017). That the fill is not *per-edge* evidence does not make it uninformative: flipping all 21,699 default edges to −1 takes masked ULM from 96 hits to 62, so the default +1 is the single most consequential sign decision in the network (V-021).

### The size law

Masking each of the 613 scored TFs in turn on a single-cell dataset, and fitting on the 227 with a self-edge, a measured own gene and a regulon surviving `tmin=3`:

```
log₁₀|Δscore| ≈ −0.53 − 0.46·log₁₀(regulon size)
Spearman(size, |Δscore|) = −0.868, p = 1.7e-70 ;  R² = 0.80 for the joint fit
```

Detection is uncorrelated with the width on its own (Spearman +0.026, p = 0.69). It carries a small coefficient in the joint fit that is not separately quantified here; its real role is the second factor, the own gene's discriminative signal (V-027).

In hit-rate terms on bulk knockTF, self-edge removal costs 51% of hits in the smallest tercile of self-edge experiments (median 45 targets), 32% in the middle (median 98) and 3% in the largest (median 425) — a split within the self-edge group only (V-010). In AUROC terms on the single-cell dataset, |ΔAUROC| runs 1.1e-2 at 7 targets to 1.3e-3 above 300 (V-027).

So the leak is not a constant property of a benchmark. It is **channel width × own-gene signal**, and both must be checked before a benchmark is trusted or dismissed.

## 3. What survives

With the perturbed gene masked, ULM+CollecTRI scores **0.344** on the 279.

The right null is not 1/N: large regulons rank high in every experiment regardless of what was perturbed. We use a **per-TF base rate** — how often the same TF is top-10% where a *different* TF was perturbed — which on these 279 is **0.122** under the masking convention the 0.344 itself uses, with three independent nulls agreeing at 0.09–0.13 (the same construction over all 643 scored experiments gives 0.123; a single three-decimal figure had been quoted against both denominators). Masked ULM is 2.8× that, and on the 643 scored set it reaches 13/38 even on experiments whose own gene was never measured (V-022).

Two apparent improvements did not survive masking:

- **Propagating scores along TF→TF edges** reached 57% against 45%. Masked it falls to 28% — *below* masked ULM's 34%. Of the hits lost, 84–90% go through the upstream regulators' regulons and 10–16% through the self-edge (V-010). Masking a co-regulated gene with a comparably negative logFC instead costs propagation 0.3 points, against 29 for the TF's own (F-004, F-012/V-012).
- **A robust (Huber) loss** scored 32% against 45% unmasked. About half to two thirds of that deficit was the robust loss *clipping the leak* — it assigns the own gene a median weight of 0.09. Masked, the comparison is 30% against 34%: still worse on paired rank (p = 0.002) and top-5%, but **not distinguishable at top-10%** (12 hits, exact McNemar p = 0.10). The deficit is confined to the 197 self-edge experiments; without a self-edge Huber and ULM both score 20/82 (F-008/V-008).

The second generalises as a warning: any outlier-robust scorer will down-weight a ~5 sd own-gene value and look worse on an unmasked benchmark for reasons unrelated to activity inference.

## 4. A second confound: study context

The 34% is TF-specific only in part. On the 106 of 643 scored experiments that support the contrast, under a symmetric mask (every perturbed TF's gene zeroed in every row):

| | Top-10% rate |
|---|---|
| The TF, in its own experiment | 0.255 |
| The same TF, same series, a *different* TF perturbed | 0.218 |
| The same TF, same cell line, a different series | 0.145 |

The study-context lift is **+0.073**, 95% study-clustered CI [+0.019, +0.178] — two thirds of the perturbed-TF effect (ratio 0.68, CI [0.21, 1.60]). Under the asymmetric mask the same interval sits on zero and is retired (F-028/V-028).

**It is removable.** Series-mean centering takes the lift to +0.020 [−0.034, +0.088], at a cost of 0.036 on the own-arm baseline. "Regulon-specific" does not imply "not batch": a series-wide signed component aligned with a TF's regulon is both (V-028).

It *is* specific to the TF's real regulon — a degree-preserving shuffle leaves 6% of it, and judging a random real TF instead gives p = 0.003 — and it is **not** regulon overlap with the series' other perturbed TF, being larger where the two TFs share no targets at all.

Three limits. It exists only under per-row weighting (pair-weighted it is −0.001). It is **absent in knockTF's one systematic 78-TF screen** (−0.004) while worth +0.117 in the 23 small studies that chose two to four related TFs — a pattern favouring the benign explanation that researchers pick functionally related TFs to study together, and that a good method *should* respond to all of them. And 537 of the 643 scored experiments cannot be tested this way at all. The data cannot separate co-selection from a shared control sample, because knockTF as fetched exposes no control-sample identifier — its 14 observation columns carry no sample or control ID (V-028; verified directly against `data/knocktf_full.h5ad`).

**Practical consequence:** benchmark on knockTF with a study control or a series centering, and do not read top-10% recovery as TF-specific detection at the resolution of a GEO series.

## 5. Does it transfer?

We reran a 2026 single-cell TF-activity benchmark on its own shipped dataset (TP53, 617 cells). The anchor and the masked condition ran their **untouched** scoring and evaluation code; the network-side conditions and the shuffled band used a patched copy whose diff was audited as touching no scoring code (F-027, audit in V-027). The anchor reproduces their shipped numbers, and an independent `decoupler` call matched them to six decimals and their score matrix to 3.55e-15 — which is expected for a closed-form t-statistic on a shipped input, so it is a determinism and version-stability check rather than a remarkable agreement (V-027).

**Masking TP53's own gene moved their headline AUROC by −0.0001**, against a refutation threshold of 0.05 stated in the script docstring before the run. The leak does not transfer here, and the size law says why: TP53 has 772 measured targets and its own transcript separates the arms at AUROC 0.493 — chance. Both factors are at their minimum (F-027/V-027).

This is a negative result for our own hypothesis. Its scope: one dataset, one TF, and the worst available case for the leak. Their other 16 datasets are not public and are likely to include small-regulon TFs with stronger knockdowns, where the law predicts ~1e-2 AUROC rather than 1e-4.

The controls found two problems the main test did not:

- **A shape-matched random regulon scores 0.47 ± 0.14** on that dataset (central 95% [0.225, 0.762], 100 seeds). That band covers **34 of their 60 published columns** for this TF, including their own `random` control — on this one dataset, for this one TF; their published ranking is over 40 datasets and this licenses no statement about it (V-027).
- **`decoupler`'s ORA is unstable on sparse counts.** Rescoring an *unchanged* matrix under different feature orderings moves its AUROC by 0.17–0.19 — the size of the effect masking appeared to produce. Verified in V-027; not yet filed as a claim (T-063).

## 6. What is new

The bare observation is not. Two groups found that a TF's own expression predicts which TF was perturbed as well as or better than regulon methods: Trescher & Leser (2019, *Scientific Reports*), who read it as the network methods underperforming and who name TF self-regulation as a cause of method *failure*; and Yashar et al. (2024, *iScience*), whose own-expression-only variant of Priori beat eleven methods including ULM, and who used it constructively as a feature. Neither framed it as a validity problem; neither masked the gene (F-025, **unverified**).

A title-and-abstract screen of a 176-paper citation graph seeded on knockTF, plus full reads of two reviews on validating these methods, did not find the self-edge mechanism, the masking ablation, a degree-preserving shuffle or a per-TF base rate. That is a screen, not an exhaustive search; one of the two reviews is not in `refs/` and so is not re-checkable here; and F-025's own confidence in the negative is **medium**, with **low** confidence that the mechanism is unpublished outside that corpus (F-025, **unverified**). On present evidence the contribution is the mechanism and its controls, not the discovery that the benchmark is gameable.

## 7. The protocol

Reusable independently of any finding here:

1. **Mask the perturbed TF's own gene** before scoring. For a univariate scorer, deleting its self-edge is equivalent at top-10% and free at scoring time; check top-5% separately.
2. **Shuffle the network preserving degree and sign mix**, with at least **100 seeds**, reported as mean ± sd with a central 95% interval. Ten seeds is not a null: it produced two wrong intervals in this work before the rule was adopted.
3. **Compare against a per-TF base rate**, not 1/N. Chance here is 12%, not 10%.
4. **Study-control or centre by series** when experiments come from multi-TF studies.
5. **Declare the causal route.** *TF-perturbing* means the TF's own mRNA is on the causal path (knockTF, Perturb-seq knockdowns); *TF-sparing* means it is not (phosphoproteomics, upstream-regulator perturbations, signalling stimuli). Agreement between two TF-perturbing benchmarks is one measurement twice. Every benchmark in this paper — including §5's — is TF-perturbing.

## 8. Limitations

- **Everything is one route.** All but one claim is TF-perturbing, mostly one benchmark. The single TF-sparing measurement is inconclusive: masked ULM against phosphoproteomics-derived occupancy gives −0.018, CI [−0.33, +0.33], on 17 TFs. That route's ground truth covers 174 of CollecTRI's 1,185 TFs, so its interval caps at about ±0.11 even if every one were measurable (F-024/V-024).
- **The within-study margin is unresolved, not small.** +0.066 under the symmetric mask, with a lower bound that is positive in only 2 of 20 bootstrap seeds and must not be cited as excluding zero. Power at a true effect of +0.05 is 28%; the minimum detectable effect is +0.110 (V-026, V-028).
- **What we could not test at all** is whether these scores track activity that abundance cannot explain. That needs the same TF at several protein levels without perturbing the TF. The largest within-TF dose resource we found, scTF-seq (Liu et al. 2025), is itself TF-perturbing and reads dose off the transcript rather than the protein (`refs/README.md`); we did not find one that meets the requirement, which is weaker than a claim that none exists.
- **Pre-registration is by docstring, not by commit.** Predictions were written into each script before its run, but the repository is not under version control, so that ordering cannot be independently verified.

## Data and code

knockTF and CollecTRI as fetched by `decoupler`; the 2026 benchmark's own repository at commit `b904196`. One script per experiment in `analysis/`, each self-contained, each ending in an assertion. Claims and verdicts in `findings.md`, append-only, every claim independently rerun before it becomes citable.
