localrules:
    create_bcftools_split_groups,


rule create_bcftools_split_groups:
    """
    Generate an annotation file telling bcftools +split how to
    create split family vcf files, to avoid multiple processing.
    """
    output:
        groups="results/slivar/group_split_file.tsv",
    params:
        samples=tc.get_probands_with_structure(gvcf_manifest),
    run:
        with open(output.groups, "w") as f:
            f.writelines(
                [
                    "{}\t-\tfamily_{}_dv_joint_calls\n".format(x, x.split("-")[2])
                    for x in params.samples
                ]
            )


rule slice_joint_callset:
    """
    Extract trio members from the joint callset.
    """
    input:
        vcf="results/glnexus/all/merged_callset.filtered.regions.csq.vcf.gz",
        groups="results/slivar/group_split_file.tsv",
    output:
        vcf=expand(
            "results/slivar/family_{family_cluster}_dv_joint_calls.vcf.gz",
            family_cluster=[
                x.split("-")[2] for x in tc.get_probands_with_structure(gvcf_manifest)
            ],
        ),
    params:
        outdir="results/slivar",
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    shell:
        "bcftools +split -Oz -o {params.outdir} -G {input.groups} {input.vcf}"


rule slivar_filter_dnm:
    """
    Use slivar logic to perform trio filtering adjusted for various brentp-assorted criteria
    """
    input:
        vcf="results/slivar/family_{family_cluster}_dv_joint_calls.vcf.gz",
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
        gnomad_popmax_af=0.01,
        gnomad_nhomalt_max=10,
        topmed_g=" -g reference_data/slivar/{}/topmed.zip ".format(reference_build)
        if "topmed-zip" in config["slivar"][reference_build]
        else "",
        topmed_filter=" " if "topmed-zip" in config["slivar"][reference_build] else "",
    benchmark:
        "results/performance_benchmarks/slivar_filter_trios/{family_cluster}/putative_dnm.tsv"
    conda:
        "../envs/slivar.yaml" if not use_containers else None
    shell:
        "slivar expr --js {input.js} -g {input.gnomad} {params.topmed_g} --vcf {input.vcf} --ped {input.ped} -x {input.bed} "
        '--info \'INFO.impactful && INFO.gnomad_popmax_af < {params.gnomad_popmax_af} {params.topmed_filter} && variant.filter == "PASS" && variant.ALT[0] != "*" \' '
        "--family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < {params.gnomad_popmax_af} {params.topmed_filter} ' "
        "--family-expr 'x_denovo:variant.CHROM == \"chrX\" && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < {params.gnomad_popmax_af} {params.topmed_filter} ' "
        "--trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < {params.gnomad_nhomalt_max} && "
        "kid.het && mom.hom_ref && dad.hom_ref && "
        "kid.DP > {params.dp_min} && mom.DP > {params.dp_min} && dad.DP > {params.dp_min} && "
        "(mom.AD[1] + dad.AD[1]) == 0 && "
        "kid.AB > {params.ab_het_min} && kid.AB < {params.ab_het_max} && "
        "mom.AB < {params.ab_homref_max} && dad.AB < {params.ab_homref_max}' "
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
