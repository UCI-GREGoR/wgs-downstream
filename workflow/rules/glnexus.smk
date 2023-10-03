localrules:
    glnexus_create_gvcf_list,


rule glnexus_create_gvcf_list:
    """
    Emit target set of gvcf filenames to file in preparation
    for glnexus
    """
    output:
        tsv="results/glnexus/{subset}/gvcf_list.tsv",
    params:
        gvcfs=lambda wildcards: tc.get_valid_subjectids(
            wildcards,
            gvcf_manifest["sampleid"].to_list(),
            "results/gvcfs/" if wildcards.subset == "all" else "results/deeptrio/",
            ".g.vcf.gz" if wildcards.subset == "all" else ".sorted.g.vcf.gz",
            wildcards.subset != "all",
        ),
    run:
        with open(output.tsv, "w") as f:
            valid_targets = []
            for x in params.gvcfs:
                if (
                    wildcards.subset == "all"
                    or (
                        "-{}-".format(wildcards.subset) in x
                        and re.search("-[012].sorted.g.vcf.gz$", x)
                    )
                    or "/{}.g.vcf.gz".format(wildcards.subset) in x
                ):
                    valid_targets.append(x)
            f.writelines([x + "\n" for x in valid_targets])


rule glnexus_joint_calling:
    """
    Given gvcfs, create a joint called dataset.
    """
    input:
        gvcfs=lambda wildcards: tc.get_valid_subjectids(
            wildcards,
            gvcf_manifest["sampleid"].to_list(),
            "results/gvcfs/" if wildcards.subset == "all" else "results/deeptrio/",
            ".g.vcf.gz" if wildcards.subset == "all" else ".sorted.g.vcf.gz",
            wildcards.subset != "all",
        ),
        tsv="results/glnexus/{subset}/gvcf_list.tsv",
        calling_ranges=lambda wildcards: tc.get_calling_range_by_chrom(
            wildcards, config["glnexus"]["calling-ranges"]
        ),
    output:
        bcf=temp("results/glnexus/{subset}/merged_callset_{chrom}.bcf"),
    params:
        gvcf_manifest=gvcf_manifest,
        memlimit=int(
            config_resources["glnexus"]["memory"]
            / config_resources["glnexus"]["threads"]
            / 1000
        ),
        glnexus_config=config["glnexus"]["config"],
        tmp="{}/results/glnexus/{{subset}}/tmp_{{chrom}}".format(tempDir),
    benchmark:
        "results/performance_benchmarks/glnexus_joint_calling/{subset}/{chrom}.tsv"
    conda:
        "../envs/glnexus.yaml"
    threads: config_resources["glnexus"]["threads"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt
        * config_resources["glnexus"]["memory"],
        qname=rc.select_queue(
            config_resources["glnexus"]["queue"], config_resources["queues"]
        ),
        tmpdir=lambda wildcards: "{}/results/glnexus/{}".format(
            tempDir, wildcards.subset
        ),
    shell:
        "mkdir -p $(dirname {params.tmp}) && "
        "cat {input.tsv} | "
        "LD_PRELOAD=$(jemalloc-config --libdir)/libjemalloc.so.$(jemalloc-config --revision) "
        "xargs glnexus_cli --config {params.glnexus_config} "
        "--bed {input.calling_ranges} -a -m {params.memlimit} "
        "-t {threads} --dir {params.tmp} > {output.bcf}"


rule prepare_joint_calling_output:
    """
    glnexus emits unconditional bcfs; use bcftools to create gzipped vcf
    """
    input:
        bcf=lambda wildcards: expand(
            "results/glnexus/{{subset}}/merged_callset_{chrom}.bcf",
            chrom=tc.get_all_calling_ranges(config["glnexus"]["calling-ranges"]),
        ),
    output:
        vcf=temp("results/glnexus/{subset}/merged_callset.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/prepare_joint_calling_output/{subset}/merged_callset.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb=4000,
        qname="small",
    shell:
        "bcftools concat {input.bcf} --threads {threads} -O z -o {output.vcf}"


rule filter_joint_calling_output:
    """
    once joint calling is complete, apply standard group filters.
    the sed conversion from 1/0 to 0/1 is a nonsensical thing that was
    required to provide single sample vcfs to Moon. don't know if that's
    still pertinent, but since development is going to stop, now's the time.
    """
    input:
        "results/glnexus/{subset}/merged_callset.vcf.gz",
    output:
        "results/glnexus/{subset}/merged_callset.filtered.vcf.gz",
    params:
        pipeline_version=pipeline_version,
        reference_build=tc.format_reference_build(reference_build),
    benchmark:
        "results/performance_benchmarks/filtered_joint_calling_output/{subset}/merged_callset.filtered.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb=4000,
        qname="small",
    shell:
        'bcftools annotate -h <(echo -e "##wgs-downstreamVersion={params.pipeline_version}\\n##reference={params.reference_build}") -O u {input} | '
        'bcftools view -i \'(FILTER = "PASS" | FILTER = ".")\' -O u | '
        "bcftools norm -m -both -O u | "
        "bcftools view -i 'FORMAT/DP >= 10 & FORMAT/GQ >= 20' -O v | "
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


rule remove_snv_region_exclusions:
    """
    Once SNV output data have had hard filters applied, further remove configurable exclusion regions.
    These are intended to be pulled from https://github.com/Boyle-Lab/Blacklist
    """
    input:
        vcf="results/glnexus/{subset}/merged_callset.filtered.vcf.gz",
        bed="reference_data/references/{}/ref.exclusion.regions.bed".format(
            reference_build
        ),
    output:
        vcf="results/glnexus/{subset}/merged_callset.filtered.regions.vcf.gz",
    benchmark:
        "results/performance_benchmarks/remove_snv_region_exclusions/{subset}/results.tsv"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        qname="small",
    shell:
        "bedtools intersect -a {input.vcf} -b {input.bed} -wa -v -header | bgzip -c > {output}"
