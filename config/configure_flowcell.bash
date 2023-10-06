#!/usr/bin/env bash
if [[ "$#" -ne 1 ]] ; then
    echo "usage: ./configure_flowcell.bash [flowcell-export-path]"
    exit 1
fi

FLOWCELL_PATH="$1"
## this command is now modified to expect the output structure of wgs-pipeline remote data export, in anticipation of switching to that as
## the primary data source.
find "${FLOWCELL_PATH}" -name "*cram" -exec readlink -f {} \; | sed -r 's|^(.*)/([^/_]+)(.*).cram|\2\t\1/\2\3.cram|' >> manifest_cram.tsv
find "${FLOWCELL_PATH}" -name "*snv.g.vcf.gz" -exec readlink -f {} \; | sed -r 's|^(.*)/([^/_]+)(.*).snv.g.vcf.gz|\2\t\1/\2\3.snv.vcf.gz\t\1/\2\3.snv.g.vcf.gz|' >> manifest_gvcf.tsv
