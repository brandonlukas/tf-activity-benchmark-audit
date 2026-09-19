# References

Full texts fetched by `fetch.sh` (Europe PMC XML converted to Markdown by `jats2md.py`). Publisher and bioRxiv PDF endpoints block scripted downloads, so entries marked *DOI only* need a manual save. Two citation graphs are here, not papers: `ariadne-corpus.md` (944 papers, seeded on perturbation-prediction and foundation models) is where the selection above came from, and `ariadne-corpus-perturbseq.md` (473 in canon, plus groundwork and since sections; seeded on Replogle 2022 genome-scale Perturb-seq, the UPR multiplexed screen and the genetic-interaction-manifold paper) was added 2026-09-18 for the dose-resolved calibration question. `ariadne-corpus-knocktf.md` (176 papers, seeded on knockTF and knockTF 2.0 themselves) was added the same day to ask who else benchmarks on knockTF, whether any of them mask the perturbed TF's own gene, and whether the abundance leak has been reported before.

## Directly about this project (TF activity inference and its benchmarks)

| File | Paper | Why it matters here |
|---|---|---|
| `badia-i-mompel-2022-decoupler.md` | decoupleR (Badia-i-Mompel 2022) | The ULM/MLM implementations and the original knockTF benchmark setup we reuse. |
| `muller-dott-2023-collectri.md` | CollecTRI (Müller-Dott 2023) | The network. Also the source of the sign annotations and autoregulatory edges that carry own-mRNA signal into ULM. |
| `feng-2020-knocktf.md` | knockTF (Feng 2020) | The benchmark set. Every experiment is a knockdown, so the perturbed TF's own mRNA is the strongest signal in it. |
| `garcia-alonso-2019-dorothea-benchmark.md` | DoRothEA benchmark (Garcia-Alonso 2019) | Earlier TF-activity benchmark; established the perturbation-recovery evaluation that knockTF-style tests inherit. |
| `tf-activity-benchmark-perturbation-2026.md` | Zhu, Han & Wang 2026, Brief Bioinform | Benchmarks 8 single-cell TF-activity methods on 40 Perturb-seq datasets; decoupleR+CollecTRI ranks highly. Selects perturbations by the sign of the TF's own expression change and does not mask that gene (released code, shipped data and Table S1 agree, V-009). Whether that inflates its results is untested; for ULM the channel is the self-edge only. Code: github.com/zhuzhe0011/Benchmarking-Methods-for-Inferring-Single-cell-Transcription-Factor-Activity (not runnable as published). |

## Calibration data candidates (abundance-only perturbations)

| File | Paper | Why it matters here |
|---|---|---|
| `replogle-2022-genome-scale-perturbseq.md` | Replogle 2022, Cell | Genome-wide CRISPRi in K562 and RPE1. Largest clean source of per-TF "activity change per unit mRNA change" slopes. |
| *DOI only* 10.1016/j.cell.2016.11.038 | Dixit 2016, Perturb-seq | Original method; author manuscript only on PMC. |
| *DOI only* 10.1101/2025.06.11.659105 | X-Atlas/Orion 2025 | Genome-wide Perturb-seq with dose information; would let slopes be fit across knockdown strength. |
| *DOI only* 10.1016/j.cell.2026.08.002 | CD4 T-cell genome-scale Perturb-seq 2026 | Primary cells, a context check for slopes fit in cell lines. |
| *DOI only* 10.1016/j.cell.2017.10.049 | L1000 / CMap (Subramanian 2017) | shRNA and CRISPR knockdown signatures across many cell lines; the "which perturbagen made this signature" query is the nearest existing analogue of inverse inference. |

## Perturbation models and why linear baselines matter

| File | Paper | Why it matters here |
|---|---|---|
| `ahlmann-eltze-2025-linear-baselines.md` | Ahlmann-Eltze, Huber & Anders 2025, Nat Methods | Deep perturbation models don't beat linear baselines. Support for keeping the method a weighted mean plus a subtraction. |
| `simple-controls-2025-exceed-deep-learning.md` | Bioinformatics 2025 | Same conclusion from a different group, with foundation-model embeddings tested. |
| `perturbench-2025.pdf` | PerturBench (arXiv 2408.10609) | Benchmark design for perturbation prediction; useful for how they split by context and what metrics they report. |
| *DOI only* 10.1101/2025.10.20.683304 | "Do outperform on well-calibrated metrics" 2025 | The rebuttal to the baseline papers. Read for the metric argument, which applies to any masked-benchmark design. |
| *DOI only* 10.1101/2024.12.23.630036 | Wu 2024 systematic comparison | Broad comparison of perturbation response models. |
| *DOI only* 10.1016/j.cell.2026.07.052 | State (Arc) 2026 | Current best-known perturbation model; background only. |
| *DOI only* 10.1016/j.cell.2025.06.008 | Virtual Cell Challenge 2025 | Defines the community evaluation; background only. |
| *DOI only* 10.64898/2026.02.04.703804 | Virtual cells need context, not scale 2026 | Argues context matters more than data scale; relevant to whether slopes transfer across cell types. |
| `aibar-2017-scenic.md` | SCENIC (Aibar 2017) | De novo regulons from expression; the main alternative to a prior network. |

## Prior art on the abundance leak (added 2026-09-18, F-025)

| File | Paper | Why it matters here |
|---|---|---|
| `trescher-2019-tf-activity-knockdown.md` | Trescher & Leser 2019, Sci Rep | Closest precedent to F-005's "own logFC alone wins". Different methods and data; frames it as network methods underperforming, not as benchmark inflation. Names TF self-regulation as a cause of method FAILURE, the opposite direction to F-010. |
| `yashar-2024-priori.md` | Yashar et al. 2024, iScience (Priori) | An own-expression-only variant beats 11 methods including decoupleR ULM on a 124-experiment benchmark. Used constructively as a feature, not as a critique. No masking ablation. |
| `karamveer-2024-grn-benchmarking-review.md` | Karamveer & Uzun 2024 | A review of how to benchmark GRN methods. Lists knockTF as ground truth and is silent on the perturbed gene, masking and autoregulation. The clean negative for F-025. |
| `feng-2023-knocktf2.md` | Feng et al. 2023, NAR (KnockTF 2.0) | 1468 datasets, 612 TFs + 172 TcoFs, adds mouse and plants. No knockdown-efficiency error estimate and no titration, so it does not fix the dose problem F-018/F-022 name. Download: licpathway.net/KnockTFv2/download.php |
| *DOI only* 10.1002/pmic.202200462 | Hecker et al. 2023, Proteomics | TFA-tools review with a dedicated "Limitations of TFA validation" section that never names the perturbed TF's own transcript. Fetch failed; scout read an open-access mirror. |

## Alternative routes to activity beyond abundance (added 2026-09-18)

| File | Paper | Why it matters here |
|---|---|---|
| `epiregulon-2025.md` | Epiregulon 2025, Nat Commun | Possible direction-changer: builds GRNs from paired scATAC+scRNA explicitly because expression-only methods "neglect post-transcriptional modulation of TFs". An already-published alternative to calibrating an abundance slope. Unread. |
| `liu-2025-sctf-seq.md` | Liu et al. 2025, Nat Genet (scTF-seq) | 384 mouse TFs, dox-inducible, continuous per-cell dose from barcode UMIs. The largest within-TF dose resource found, but TF-perturbing by this project's route taxonomy and dose is a transcript proxy, not protein, so it does not satisfy T-049. |

## Activity without abundance (added 2026-09-18)

| File | Paper | Why it matters here |
|---|---|---|
| `sugimoto-2026-tfactprofiler.md` | Sugimoto et al. 2026, Nucleic Acids Research (TFActProfiler) | The closest existing work to this project's question, and a possible direction-changer. Benchmarks mRNA-derived TF activity against PHOSPHOPROTEOMICS-derived activity (SIGNOR-annotated activating/inhibitory phosphosites) in an MPXV time course in human fibroblasts, i.e. ground truth for activity that abundance cannot explain. Reports TF transcript abundance vs TF phosphorylation correlates -0.02 (95% CI -0.32 to 0.29), and that CollecTRI+ULM tracks the phospho trajectories significantly worse than their own resource. Caveats not yet checked: 26 phosphosite profiles, one dataset, authors' own method winning. |
