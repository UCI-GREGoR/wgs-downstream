localrules:
    expansionhunter_denovo_clone_repo,


rule expansionhunter_denovo_clone_repo:
    """
    Clone a local copy of expansionhunterdenovo because
    the conda package does not contain required auxiliary scripts.
    """
    output:
        directory("results/expansionhunter_denovo/repo"),
    params:
        github="https://github.com/Illumina/ExpansionHunterDenovo.git",
        version=config["expansionhunter_denovo"]["version"],
    shell:
        "git clone {params.github} {output} && "
        "cd {output} && "
        "git checkout {params.version}"


checkpoint expansionhunter_denovo_create_manifest:
    """
    Parse manifest and data model to get a set of subjects
    with affected status and bams, and create an expansionhunter
    denovo manifest containing just those subjects
    """
    input:
        affected_status=config["expansionhunter_denovo"]["affected-status"],
        linker="results/linker.tsv",
    output:
        tsv="results/expansionhunter_denovo/manifest.tsv",
    params:
        projectids=bam_manifest["projectid"],
        sampleids=bam_manifest["sampleid"],
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        qname="small",
    script:
        "../scripts/create_expansionhunter_denovo_manifest.R"


rule expansionhunter_denovo_profile_subject:
    """
    Run single-sample STR profiling with expansionhunter_denovo
    """
    input:
        bam="results/bams/{projectid}/{sampleid}.bam",
        bai="results/bams/{projectid}/{sampleid}.bai",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
    output:
        json="results/expansionhunter_denovo/profiles/{projectid}/{sampleid}.str_profile.json",
        motif="results/expansionhunter_denovo/profiles/{projectid}/{sampleid}.motif.tsv",
        locus="results/expansionhunter_denovo/profiles/{projectid}/{sampleid}.locus.tsv",
    params:
        outprefix="results/expansionhunter_denovo/profiles/{projectid}/{sampleid}",
        min_anchor_mapq=50,
        max_irr_mapq=40,
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        qname="small",
    shell:
        "ExpansionHunterDenovo profile --reads {input.bam} "
        "--reference {input.fasta} "
        "--output-prefix {params.outprefix} "
        "--min-anchor-mapq {params.min_anchor_mapq} "
        "--max-irr-mapq {params.max_irr_mapq}"


rule expansionhunter_denovo_merge_profiles:
    """
    Combine single-sample STR profiling from expansionhunter_denovo
    into a multi-sample STR profile
    """
    input:
        manifest="results/expansionhunter_denovo/manifest.tsv",
        jsons=lambda wildcards: tc.select_expansionhunter_denovo_subjects(
            wildcards,
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/expansionhunter_denovo/profiles",
            "str_profile.json",
        ),
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
    output:
        "results/expansionhunter_denovo/profiles/combined_profiles.multisample_profile.json",
    params:
        outprefix="results/expansionhunter_denovo/profiles/combined_profiles",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        qname="small",
    shell:
        "ExpansionHunterDenovo merge --reference {input.fasta} "
        "--manifest {input.manifest} "
        "--output-prefix {params.outprefix}"


rule expansionhunter_denovo_locus_outliers:
    """
    Use an expansionhunter_denovo multi-sample STR profile
    to compute locus-specific outliers, in the absence of
    dataset-specific case/control data.
    """
    input:
        manifest="results/expansionhunter_denovo/manifest.tsv",
        jsons=lambda wildcards: tc.select_expansionhunter_denovo_subjects(
            wildcards,
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/expansionhunter_denovo/profiles",
            "str_profile.json",
        ),
        combined_json="results/expansionhunter_denovo/profiles/combined_profiles.multisample_profile.json",
        repo="results/expansionhunter_denovo/repo",
    output:
        "results/expansionhunter_denovo/outlier_analysis/results.outlier_locus.tsv",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        qname="small",
    shell:
        "{input.repo}/scripts/outlier.py locus --manifest {input.manifest} --multisample-profile {input.combined_json} --output {output}"


rule expansionhunter_denovo_motif_outliers:
    """
    Use an expansionhunter_denovo multi-sample STR profile
    to compute motif-specific outliers, in the absence of
    dataset-specific case/control data.
    """
    input:
        manifest="results/expansionhunter_denovo/manifest.tsv",
        jsons=lambda wildcards: tc.select_expansionhunter_denovo_subjects(
            wildcards,
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/expansionhunter_denovo/profiles",
            "str_profile.json",
        ),
        combined_json="results/expansionhunter_denovo/profiles/combined_profiles.multisample_profile.json",
        repo="results/expansionhunter_denovo/repo",
    output:
        "results/expansionhunter_denovo/outlier_analysis/results.outlier_motif.tsv",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        qname="small",
    shell:
        "{input.repo}/scripts/outlier.py motif --manifest {input.manifest} --multisample-profile {input.combined_json} --output {output}"
