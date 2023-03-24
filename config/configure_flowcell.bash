#!/usr/bin/env bash
if [[ "$#" -ne 1 ]] ; then
    echo "usage: ./configure_flowcell.bash [flowcell-workflow-path]"
    exit 1
fi

FLOWCELL_PATH="$1"

find "${FLOWCELL_PATH}/results/bqsr" -name "*.bam" -exec readlink -f {} \; | sed -r 's|^(.*)/(RU[0-9]+)/(SQ[0-9]+).bam|\2\t\3\t\1/\2/\3.bam|' >> manifest_bam.tsv
find "${FLOWCELL_PATH}/results/deepvariant" -name "*.g.vcf.gz" -exec readlink -f {} \; | sed -r 's|^(.*)/(RU[0-9]+)/(SQ[0-9]+).sorted.g.vcf.gz|\2\t\3\t\1/\2/\3.sorted.g.vcf.gz|' >> manifest_gvcf.tsv