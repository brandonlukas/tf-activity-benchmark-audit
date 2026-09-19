#!/usr/bin/env bash
# Fetch full text for each reference. Europe PMC XML -> .md (via jats2md.py); bioRxiv source XML for preprints;
# arXiv PDF where given. Publisher PDF endpoints block scripted fetches, so PDFs are not attempted. Re-runnable.
UA="Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"
while IFS='|' read -r slug doi extra; do
  [ -z "$slug" ] && continue
  ls "$slug".* >/dev/null 2>&1 && { echo "have  $slug"; continue; }
  pmc=$(curl -s -A "$UA" "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22$doi%22&format=json&resultType=lite" | grep -o '"pmcid":"PMC[0-9]*"' | head -1 | grep -o 'PMC[0-9]*')
  if [ -n "$pmc" ]; then
    curl -s -A "$UA" -o "$slug.xml" "https://www.ebi.ac.uk/europepmc/webservices/rest/$pmc/fullTextXML"
    grep -q "<article" "$slug.xml" && { python3 jats2md.py "$slug.xml" > "$slug.md" && rm "$slug.xml"; echo "pmc   $slug ($pmc)"; continue; }
    rm -f "$slug.xml"
  fi
  case "$doi" in 10.1101/*|10.64898/*)
    src=$(curl -s -A "$UA" "https://api.biorxiv.org/details/biorxiv/$doi" | grep -o '"jatsxml":"[^"]*"' | tail -1 | cut -d'"' -f4 | sed 's|\\/|/|g')
    [ -n "$src" ] && curl -s -A "$UA" -o "$slug.xml" "$src" && grep -q "<article" "$slug.xml" && { python3 jats2md.py "$slug.xml" > "$slug.md" && rm "$slug.xml"; echo "biorx $slug"; continue; }
    rm -f "$slug.xml";;
  esac
  if [ -n "$extra" ]; then
    curl -s -L -A "$UA" -o "$slug.pdf" "$extra"; file "$slug.pdf" | grep -q PDF && { echo "pdf   $slug"; continue; }; rm -f "$slug.pdf"
  fi
  echo "MISS  $slug  https://doi.org/$doi"
done <<'LIST'
ahlmann-eltze-2025-linear-baselines|10.1038/s41592-025-02772-6|
ahlmann-eltze-2024-linear-baselines-preprint|10.1101/2024.09.16.613342|
replogle-2022-genome-scale-perturbseq|10.1016/j.cell.2022.05.013|
dixit-2016-perturbseq|10.1016/j.cell.2016.11.038|
aibar-2017-scenic|10.1038/nmeth.4463|
state-2026-perturbation-across-contexts|10.1016/j.cell.2026.07.052|
state-2025-preprint|10.1101/2025.06.26.661135|
wu-2024-systematic-comparison-perturbation-models|10.1101/2024.12.23.630036|
subramanian-2017-l1000-cmap|10.1016/j.cell.2017.10.049|
xatlas-orion-2025-genome-wide-perturbseq|10.1101/2025.06.11.659105|
simple-controls-2025-exceed-deep-learning|10.1093/bioinformatics/btaf317|
well-calibrated-metrics-2025|10.1101/2025.10.20.683304|
virtual-cell-challenge-2025|10.1016/j.cell.2025.06.008|
virtual-cells-need-context-2026|10.64898/2026.02.04.703804|
perturbench-2025|10.52202/085713-3225|https://arxiv.org/pdf/2408.10609
cd4-perturbseq-2026|10.1016/j.cell.2026.08.002|
tf-activity-benchmark-perturbation-2026|10.1093/bib/bbag513|
sugimoto-2026-tfactprofiler|10.1093/nar/gkag897|
trescher-2019-tf-activity-knockdown|10.1038/s41598-019-46053-7|
yashar-2024-priori|10.1016/j.isci.2024.109124|
feng-2023-knocktf2|10.1093/nar/gkad1016|
hecker-2023-tfa-tools-review|10.1002/pmic.202200462|
karamveer-2024-grn-benchmarking-review|10.1177/11779322241287120|https://arxiv.org/pdf/2307.08463
epiregulon-2025|10.1038/s41467-025-62252-5|
liu-2025-sctf-seq|10.1038/s41588-025-02343-7|
sugimoto-2025-tfactprofiler-preprint|10.1101/2025.10.05.680506|
muller-dott-2023-collectri|10.1093/nar/gkad841|
badia-i-mompel-2022-decoupler|10.1093/bioadv/vbac016|
feng-2020-knocktf|10.1093/nar/gkz881|
garcia-alonso-2019-dorothea-benchmark|10.1101/gr.240663.118|
LIST
