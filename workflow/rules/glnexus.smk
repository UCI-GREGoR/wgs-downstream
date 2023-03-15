localrules:
    glnexus_create_gvcf_list,


rule glnexus_create_gvcf_list:
    """
    Emit target set of gvcf filenames to file in preparation
    for glnexus
    """
    input:
        gvcfs=gvcf_manifest["gvcf"],
    output:
        tsv="results/glnexus/gvcf_list.tsv",
    run:
        with open(output.tsv, "w") as f:
            f.writelines([x + "\n" for x in input.gvcfs])


rule glnexus_joint_calling:
    """
    Given gvcfs, create a joint called dataset.
    """
    input:
        tsv="results/glnexus/gvcf_list.tsv",
        calling_ranges=config["glnexus"]["calling-ranges"],
    output:
        bcf=temp("results/glnexus/merged_callset.bcf"),
    params:
        gvcf_manifest=gvcf_manifest,
        memlimit="16000",
    conda:
        "../envs/glnexus.yaml"
    threads: 4
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "LD_PRELOAD=$(jemalloc-config --libdir)/libjemalloc.so.$(jemalloc-config --revision) "
        "glnexus_cli --config DeepVariant "
        "--bed {input.calling_ranges} -a -m {params.memlimit} "
        "-t {threads} -l {input.tsv} > {output} && "
        "rm -Rf GLnexus.DB"


rule prepare_joint_calling_output:
    """
    glnexus emits unconditional bcfs; use bcftools to create gzipped vcfs
    """
    input:
        bcf="results/glnexus/merged_callset.bcf",
    output:
        vcf="results/glnexus/merged_callset.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "bcftools view {input.bcf} --threads {threads} -O z -o {output.vcf}"
