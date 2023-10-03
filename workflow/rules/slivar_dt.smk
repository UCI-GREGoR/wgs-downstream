use rule slivar_dv_filter_dnm_impactful as slivar_dt_filter_dnm_impactful with:
    input:
        vcf="results/glnexus/{family_cluster}/merged_callset.filtered.regions.csq.vcf.gz",
        js="reference_data/slivar/functions.js",
        gnomad="reference_data/slivar/{}/gnomad.zip".format(reference_build),
        topmed="reference_data/slivar/{}/topmed.zip".format(reference_build),
        bed="reference_data/slivar/{}/low.complexity.bed.gz".format(reference_build),
        ped="results/deeptrio/{family_cluster}.ped",
    output:
        vcf="results/slivar/dt/{family_cluster}/putative_dnm_impactful.vcf.gz",
    benchmark:
        "results/performance_benchmarks/slivar_dv_filter_trios/{family_cluster}/putative_dnm.tsv"


use rule slivar_dv_filter_dnm_all as slivar_dt_filter_dnm_all with:
    input:
        vcf="results/glnexus/{family_cluster}/merged_callset.filtered.regions.csq.vcf.gz",
        js="reference_data/slivar/functions.js",
        bed="reference_data/slivar/{}/low.complexity.bed.gz".format(reference_build),
        ped="results/deeptrio/{family_cluster}.ped",
    output:
        vcf="results/slivar/dt/{family_cluster}/putative_dnm_all.vcf.gz",
    benchmark:
        "results/performance_benchmarks/slivar_dt_filter_trios/{family_cluster}/putative_dnm_all.tsv"


rule slivar_dt_compound_hets:
    """
    Run slivar to compute compound heterozygotes
    """
    input:
        vcf="results/slivar/{family_cluster}/putative_dnm.vcf.gz",
        ped="results/deeptrio/{family_cluster}.ped",
    output:
        vcf="results/slivar/{family_cluster}/putative_ch.vcf.gz",
    benchmark:
        "results/performance_benchmarks/slivar_compound_hets/{family_cluster}/putative_ch.tsv"
    conda:
        "../envs/slivar.yaml" if not use_containers else None
    shell:
        "slivar compound-hets -v {input.vcf} --sample-field comphet_side --sample-field denovo -p {input.ped} -o {output.vcf}"
