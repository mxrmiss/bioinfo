#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

TAG="foot_only"

mkdir -p "results/08_enrich/inputs/${TAG}"

awk -F'\t' 'NR>1{for(i=2;i<=NF;i++) if(($i+0)>0){print $1;break}}' results/05_matrix/tpms/gene_tpm.tsv | sort -u > /tmp/detectable.genes

{ tail -n +2 results/07_annot/gene2go.tsv | cut -f1; tail -n +2 results/07_annot/gene2pathway.tsv | cut -f1; } | sort -u > /tmp/annot.genes

comm -12 /tmp/detectable.genes /tmp/annot.genes > "results/08_enrich/inputs/${TAG}/background.list"

echo "DONE: background gene list generated at: results/08_enrich/inputs/${TAG}/background.list"
