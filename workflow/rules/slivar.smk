rule slivar_filter_dnm:
    """
    Use slivar logic to perform trio filtering adjusted for various brentp-assorted criteria
    """
    input:
        vcf="results/glnexus/{family_cluster}/merged_callset.filtered.regions.csq.vcf.gz",
        js="reference_data/slivar/functions.js",
        gnomad="reference_data/slivar/{}/gnomad.zip".format(reference_build),
        topmed="reference_data/slivar/{}/topmed.zip".format(reference_build),
        bed="reference_data/slivar/{}/low.complexity.bed.gz".format(reference_build),
        ped="results/deeptrio/{family_cluster}.ped",
    output:
        vcf="results/slivar/{family_cluster}/putative_dnm.vcf.gz",
    params:
        dp_min=12,
        ab_het_min=0.25,
        ab_het_max=0.75,
        ab_homref_max=0.02,
        gnomad_popmax_af=0.00001,
        gnomad_nhomalt_max=10,
        topmed_g=" -g reference_data/slivar/{}/topmed.zip ".format(reference_build)
        if "topmed-zip" in config["slivar"][reference_build]
        else "",
        topmed_filter=" && INFO.topmed_af < 0.05 "
        if "topmed-zip" in config["slivar"][reference_build]
        else "",
    benchmark:
        "results/performance_benchmarks/slivar_filter_trios/{family_cluster}/putative_dnm.tsv"
    conda:
        "../envs/slivar.yaml" if not use_containers else None
    shell:
        "slivar expr --js {input.js} -g {input.gnomad} {params.topmed_g} --vcf {input.vcf} --ped {input.ped} -x {input.bed} "
        "--info \"INFO.impactful && INFO.gnomad_popmax_af < {params.gnomad_popmax_af} {params.topmed_filter} && variant.ALT[0] != '*' \" "
        '--family-expr "denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < {params.gnomad_popmax_af} {params.topmed_filter} " '
        "--family-expr \"x_denovo:variant.CHROM == 'chrX' && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < {params.gnomad_popmax_af} {params.topmed_filter} \" "
        '--trio "comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < {params.gnomad_nhomalt_max} && '
        "kid.het && mom.hom_ref && dad.hom_ref && "
        "kid.DP > {params.dp_min} && mom.DP > {params.dp_min} && dad.DP > {params.dp_min} && "
        "(mom.AD[1] + dad.AD[1]) == 0 && "
        "kid.AB > {params.ab_het_min} && kid.AB < {params.ab_het_max} && "
        'mom.AB < {params.ab_homref_max} && dad.AB < {params.ab_homref_max}" '
        "--pass-only -o {output.vcf} --skip-non-variable"


rule slivar_compound_hets:
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
