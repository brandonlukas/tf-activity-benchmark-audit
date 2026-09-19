# tf-activity-benchmark-audit

Exploratory research on transcription-factor activity inference with the CollecTRI network. The question is whether anything beats ULM+CollecTRI in a way that survives controls, and what it would have to model to be worth building. Current direction: activity-vs-abundance decoupling. See README.md for findings so far.

## Goals

- Keep every method interpretable: it must reduce to "weighted mean of target statistics" in the simple case.
- Treat knockTF-style benchmarks with suspicion. They reward detecting the perturbed TF's own mRNA, not its activity. The masked-gene ULM baseline is 34% top-10%, not 45%; chance for the perturbed TFs is about 12%, not 10% (per-TF base rate, V-022).
- Prefer negative results written down over positive results not controlled.

## Standing rules

- When a task finishes, reread the new findings.md entries and add a follow-up task to `tasks.md` for every unresolved question, contested or unverifiable claim, and new lead. Then pick up the next open task that does not need Brandon.
- Do not stop to ask about judgment calls. Make a reasonable decision, log it under Decisions in `tasks.md` with one line of reasoning, and keep going.
- Only stop if blocked on something only Brandon can provide: credentials, a missing resource, or a change in research direction. State exactly what is needed in one sentence.
- This does not cover destructive or outward-facing actions (deleting data, pushing, rewriting findings.md entries). Those still get asked.

## Layout

- `tasks.md` open follow-up tasks and the decision log (see Standing rules). Check off tasks, do not delete them.

- `analysis/` one script per experiment, run from inside that directory: `cd analysis && uv run python <script>.py`. Data paths are `../data/`. Later scripts import earlier ones and rerun them; that is deliberate, keep scripts self-contained rather than factoring a shared module.
- `data/` knockTF and CollecTRI as fetched by decoupler. knockTF as fetched keeps only experiments whose TF's own logFC is below -1 (all 388; V-018), so it is selected on the abundance drop by construction and has no dose range. That is decoupler's default `thr_fc=-1`; `data/knocktf_full.h5ad` is `dc.ds.knocktf(thr_fc=None)`, 907 experiments, 456 TFs, own logFC from -8.5 to +9.9. In both files an exact 0.0 is a fill value for "gene not measured", not a measurement (V-022). Gitignored; regenerate with `dc.ds.knocktf()` and `dc.op.collectri()`. CollecTRI as fetched is 185 edges short of the published network (decoupler drops rows with no `references`, including HOXA5's self-edge) and carries a `sign_decision` column: PMID, regulon, or default activation (V-016).
- `refs/` full texts (Markdown from Europe PMC JATS via `jats2md.py`), `fetch.sh` to add more, and the annotated index in `refs/README.md`. Add a row there for every paper fetched, with one line on why it matters here.
- `findings.md` append-only log of claims and verdicts (see below).

## Conventions

- Name the causal route on every claim about how well a method infers activity. **TF-perturbing**: the TF's own mRNA is on the causal path (knockTF, Perturb-seq knockdowns, anything that perturbs the TF). **TF-sparing**: it is not (phosphoproteomics, upstream-regulator perturbations, signalling stimuli). A claim is general only when both routes agree; otherwise it names the route it was tested on and says the other is untested. Everything from F-001 to F-022 is TF-perturbing, one benchmark design.
- Corroboration counts only across routes. knockTF, the 2026 Perturb-seq benchmark and Sugimoto 2026 agree with each other and share one confound (F-005, F-009, F-023); agreement between sources that select on the TF's own mRNA is not evidence, it is the same measurement repeated.
- Ten seeds is not a null. Any shuffled or permuted control quoted as a mean or a range needs at least 100 draws, reported as mean +/- sd with a central 95% interval, never a min-max over seeds. Two entries here have already been corrected for this (V-024, V-027; the second had 42 of 100 seeds below its quoted minimum).
- Compare intervals, not point estimates, when sources differ in n by an order of magnitude. 279 experiments and 19 TFs do not get an equal vote.
- Any claimed benchmark gain ships with two controls in the same script: the perturbed TF's own gene masked before scoring, and a degree-preserving shuffled network. Both have already caught false positives here.
- State the predicted number before running an experiment. Report top-10% and top-5% hit rate of the true TF, plus paired better/worse against ULM.
- Every script ends with one `assert` that fails if its logic breaks. No test framework.
- Scripts print their result table AND persist it. Every script writes its per-row results to `results/<script>.csv` and tees its stdout to `results/<script>.log`. Printing alone loses the data: sixteen of the first seventeen scripts here kept nothing, so every number survived only as prose in findings.md and no figure could be drawn without rerunning everything.
- `results/master.csv` is one row per scored experiment with every condition's rank fraction and hit flag, plus the covariates (regulon size, self-edge presence and sign, sign_decision, own logFC, study, cell line). Most headline numbers in the log are recomputable from it alone; regenerate with `cd analysis && uv run python master.py`.
- The numbers that matter still go into README.md or findings.md, not into commit messages or chat only.
- `findings.md` entries: claims are `F-nnn`, verdicts are `V-nnn` referencing a claim. Nobody edits another entry. A claim is citable once a confirmed verdict exists. Fields: Claim, Evidence, Source, Confidence with reason, Depends on.
- Python via `uv`. No new dependencies for what numpy or pandas already do.
- Never modify installed packages (nothing under `.venv/` or site-packages). If a library needs changes, copy the relevant code into the project and change it there, with a comment naming the package, version and file it came from.
- Do not add, remove, or upgrade a dependency without logging the reason under Decisions in `tasks.md` first. `pyproject.toml` and `uv.lock` change only through `uv add` / `uv remove` / `uv lock --upgrade-package`.
- To inspect a package, read its files directly or use `uv pip show <pkg>`. Do not run ad-hoc code to poke at it.
- Not a place for a "senior researcher" persona. Roles in `.claude/agents/` are narrow with structured outputs.
