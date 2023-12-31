localrules:
    create_bcftools_split_groups,


rule create_bcftools_split_groups:
    """
    Generate an annotation file telling bcftools +split how to
    create split family vcf files, to avoid multiple processing.
    """
    output:
        groups="results/slivar/dv/group_split_file.tsv",
    params:
        samples=tc.get_probands_with_structure(gvcf_manifest, True),
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

    Don't think this is actually thread-aware, but worth a shot. Without threading
    and with ~600 subjects, this takes 10 hours :(
    """
    input:
        vcf="results/glnexus/all/merged_callset.filtered.regions.csq.vcf.gz",
        groups="results/slivar/dv/group_split_file.tsv",
    output:
        vcf=expand(
            "results/slivar/dv/family_{family_cluster}_dv_joint_calls.vcf.gz",
            family_cluster=[
                x.split("-")[2]
                for x in tc.get_probands_with_structure(gvcf_manifest, True)
            ],
        ),
    params:
        outdir="results/slivar/dv",
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bcftools +split --threads {threads} -Oz -o {params.outdir} -G {input.groups} {input.vcf}"


rule slivar_dv_filter_dnm_impactful:
    """
    Use slivar logic to perform trio filtering adjusted for various brentp-assorted criteria
    """
    input:
        vcf="results/slivar/dv/family_{family_cluster}_dv_joint_calls.vcf.gz",
        js="reference_data/slivar/functions.js",
        gnomad="reference_data/slivar/{}/gnomad.zip".format(reference_build),
        topmed="reference_data/slivar/{}/topmed.zip".format(reference_build),
        bed="reference_data/slivar/{}/low.complexity.bed.gz".format(reference_build),
        ped="results/deeptrio/{family_cluster}.ped",
    output:
        vcf="results/slivar/dv/{family_cluster}/putative_dnm_impactful.vcf.gz",
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
        "results/performance_benchmarks/slivar_dv_filter_trios/{family_cluster}/putative_dnm_impactful.tsv"
    conda:
        "../envs/slivar.yaml" if not use_containers else None
    threads: config_resources["slivar"]["threads"]
    resources:
        mem_mb=config_resources["slivar"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["slivar"]["queue"], config_resources["queues"]
        ),
    shell:
        "slivar expr --js {input.js} -g {input.gnomad} {params.topmed_g} --vcf {input.vcf} --ped {input.ped} -x {input.bed} "
        "--pass-only -o {output.vcf} --skip-non-variable "
        "--info 'INFO.impactful && INFO.gnomad_popmax_af < {params.gnomad_popmax_af} {params.topmed_filter} && variant.ALT[0] != \"*\" ' "
        '--trio "denovo:kid.het && mom.hom_ref && dad.hom_ref && kid.DP > 12 && mom.DP > 12 && dad.DP > 12 && (mom.AD[1] + dad.AD[1]) == 0 '
        ' && kid.GQ > 20 && mom.GQ > 20 && dad.GQ > 20" '
        '--trio "informative:kid.GQ > 20 && dad.GQ > 20 && mom.GQ > 20 && kid.alts == 1 && ((mom.alts == 1 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 1))" '
        '--trio "recessive:trio_autosomal_recessive(kid, mom, dad)" '


rule slivar_dv_filter_dnm_all:
    """
    Use slivar logic to perform trio filtering adjusted for various brentp-assorted criteria
    """
    input:
        vcf="results/slivar/dv/family_{family_cluster}_dv_joint_calls.vcf.gz",
        js="reference_data/slivar/functions.js",
        bed="reference_data/slivar/{}/low.complexity.bed.gz".format(reference_build),
        ped="results/deeptrio/{family_cluster}.ped",
    output:
        vcf="results/slivar/dv/{family_cluster}/putative_dnm_all.vcf.gz",
    params:
        dp_min=12,
        ab_het_min=0.25,
        ab_het_max=0.75,
        ab_homref_max=0.02,
    benchmark:
        "results/performance_benchmarks/slivar_dv_filter_trios/{family_cluster}/putative_dnm_all.tsv"
    conda:
        "../envs/slivar.yaml" if not use_containers else None
    threads: config_resources["slivar"]["threads"]
    resources:
        mem_mb=config_resources["slivar"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["slivar"]["queue"], config_resources["queues"]
        ),
    shell:
        "slivar expr --js {input.js} --vcf {input.vcf} --ped {input.ped} -x {input.bed} "
        "--pass-only -o {output.vcf} --skip-non-variable "
        "--info 'variant.ALT[0] != \"*\" ' "
        '--trio "denovo:kid.het && mom.hom_ref && dad.hom_ref && kid.DP > 12 && mom.DP > 12 && dad.DP > 12 && (mom.AD[1] + dad.AD[1]) == 0 '
        ' && kid.GQ > 20 && mom.GQ > 20 && dad.GQ > 20" '


rule slivar_dv_compound_hets:
    """
    Run slivar to compute compound heterozygotes
    """
    input:
        vcf="results/slivar/dv/{family_cluster}/putative_dnm.vcf.gz",
        ped="results/deeptrio/{family_cluster}.ped",
    output:
        vcf="results/slivar/dv/{family_cluster}/putative_ch.vcf.gz",
    benchmark:
        "results/performance_benchmarks/slivar_dv_compound_hets/{family_cluster}/putative_ch.tsv"
    conda:
        "../envs/slivar.yaml" if not use_containers else None
    threads: config_resources["slivar"]["threads"]
    resources:
        mem_mb=config_resources["slivar"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["slivar"]["queue"], config_resources["queues"]
        ),
    shell:
        "slivar compound-hets -v {input.vcf} --sample-field comphet_side --sample-field denovo -p {input.ped} -o {output.vcf}"


rule slivar_summarize_dnm_counts:
    """
    For each proband, determine how many de novos were detected before and after impactfulness/frequency filtering.
    """
    input:
        dnm_all="results/slivar/{gvcf_type}/{family_cluster}/putative_dnm_all.vcf.gz",
        dnm_impactful="results/slivar/{gvcf_type}/{family_cluster}/putative_dnm_impactful.vcf.gz",
    output:
        tsv="results/slivar/{gvcf_type}/{family_cluster}/dnm_summary.tsv",
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "echo -e \"{wildcards.family_cluster}\\t$(gunzip -c {input.dnm_all} | grep denovo | awk '! /^#/' | wc -l)\\t$(gunzip -c {input.dnm_impactful} | grep denovo | awk '! /^#/' | wc -l)\" > {output.tsv}"


localrules:
    slivar_combine_dnm_count_summary,


rule slivar_combine_dnm_count_summary:
    """
    Combine the de novo counts from all probands into a single tsv.
    """
    input:
        tsv=expand(
            "results/slivar/{{gvcf_type}}/{family_cluster}/dnm_summary.tsv",
            family_cluster=[
                x.split("-")[2]
                for x in tc.get_probands_with_structure(gvcf_manifest, True)
            ],
        ),
    output:
        tsv="results/slivar/{gvcf_type}/dnm_count_summary.tsv",
    run:
        res = ["proband_id\ttotal_dnm_count\timpactful_and_rare_dnms\n"]
        for fn in input.tsv:
            with open(fn, "r") as f:
                res.extend(f.readlines())
        with open(output.tsv, "w") as f:
            f.writelines(res)
