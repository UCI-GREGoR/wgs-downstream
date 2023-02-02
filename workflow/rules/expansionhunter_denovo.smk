rule expansionhunter_denovo_profile_subject:
    """
    Run single-sample STR profiling with expansionhunter_denovo
    """
    input:
        "",
    output:
        "",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb="",
        qname="",
    shell:
        ""


rule expansionhunter_denovo_merge_profiles:
    """
    Combine single-sample STR profiling from expansionhunter_denovo
    into a multi-sample STR profile
    """
    input:
        "",
    output:
        "",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb="",
        qname="",
    shell:
        ""


rule expansionhunter_denovo_locus_outliers:
    """
    Use an expansionhunter_denovo multi-sample STR profile
    to compute locus-specific outliers, in the absence of
    dataset-specific case/control data.
    """
    input:
        "",
    output:
        "",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb="",
        qname="",
    shell:
        ""


rule expansionhunter_denovo_motif_outliers:
    """
    Use an expansionhunter_denovo multi-sample STR profile
    to compute motif-specific outliers, in the absence of
    dataset-specific case/control data.
    """
    input:
        "",
    output:
        "",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb="",
        qname="",
    shell:
        ""
