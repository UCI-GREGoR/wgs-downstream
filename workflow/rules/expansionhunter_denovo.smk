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
        sex_manifest=config["sample-sex"],
        exclusions=config["expansionhunter_denovo"]["excluded-samples"],
    output:
        tsv="results/expansionhunter_denovo/manifest.tsv",
    params:
        sampleids=reads_manifest["sampleid"],
    conda:
        "../envs/r.yaml"
    threads: config_resources["r"]["threads"]
    resources:
        mem_mb=config_resources["r"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["r"]["queue"], config_resources["queues"]
        ),
    script:
        "../scripts/create_expansionhunter_denovo_manifest.R"


rule expansionhunter_denovo_profile_subject:
    """
    Run single-sample STR profiling with expansionhunter_denovo
    """
    input:
        cram="results/crams/{sampleid}.cram",
        crai="results/crams/{sampleid}.crai",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
    output:
        json="results/expansionhunter_denovo/profiles/{sampleid}.str_profile.json",
        motif="results/expansionhunter_denovo/profiles/{sampleid}.motif.tsv",
        locus="results/expansionhunter_denovo/profiles/{sampleid}.locus.tsv",
    params:
        outprefix="results/expansionhunter_denovo/profiles/{sampleid}",
        min_anchor_mapq=50,
        max_irr_mapq=40,
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: config_resources["expansionhunter_denovo"]["threads"]
    resources:
        mem_mb=config_resources["expansionhunter_denovo"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["expansionhunter_denovo"]["queue"],
            config_resources["queues"],
        ),
    shell:
        "ExpansionHunterDenovo profile --reads {input.cram} "
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
            checkpoints,
            reads_manifest["sampleid"],
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
    threads: config_resources["expansionhunter_denovo"]["threads"]
    resources:
        mem_mb=config_resources["expansionhunter_denovo"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["expansionhunter_denovo"]["queue"],
            config_resources["queues"],
        ),
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
            checkpoints,
            reads_manifest["sampleid"],
            "results/expansionhunter_denovo/profiles",
            "str_profile.json",
        ),
        combined_json="results/expansionhunter_denovo/profiles/combined_profiles.multisample_profile.json",
        repo="results/expansionhunter_denovo/repo",
    output:
        "results/expansionhunter_denovo/outlier_analysis/results.outlier_locus.tsv",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: config_resources["expansionhunter_denovo"]["threads"]
    resources:
        mem_mb=config_resources["expansionhunter_denovo"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["expansionhunter_denovo"]["queue"],
            config_resources["queues"],
        ),
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
            checkpoints,
            reads_manifest["sampleid"],
            "results/expansionhunter_denovo/profiles",
            "str_profile.json",
        ),
        combined_json="results/expansionhunter_denovo/profiles/combined_profiles.multisample_profile.json",
        repo="results/expansionhunter_denovo/repo",
    output:
        "results/expansionhunter_denovo/outlier_analysis/results.outlier_motif.tsv",
    conda:
        "../envs/expansionhunter_denovo.yaml"
    threads: config_resources["expansionhunter_denovo"]["threads"]
    resources:
        mem_mb=config_resources["expansionhunter_denovo"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["expansionhunter_denovo"]["queue"],
            config_resources["queues"],
        ),
    shell:
        "{input.repo}/scripts/outlier.py motif --manifest {input.manifest} --multisample-profile {input.combined_json} --output {output}"


rule annovar_extract_tarball:
    """
    Take a user-specified copy of the annovar release tarball
    and extract it into the repo.

    Note that due to licensing restrictions on the use of annovar,
    this step requires the user to get the tarball.
    """
    output:
        directory("results/annovar"),
    params:
        tarball=lambda wildcards: config["expansionhunter_denovo"]["annovar-tarball"]
        if "annovar-tarball" in config["expansionhunter_denovo"]
        else "",
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "tar -C results -xvzf {params.tarball}"


rule annovar_download_references:
    """
    Get a copy of the annovar grch38 reference data.
    """
    input:
        "results/annovar",
    output:
        directory("results/annovar_humandb"),
    conda:
        "../envs/annovar.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "perl {input}/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene {output}"


rule expansionhunter_denovo_annotate:
    """
    Using an installation of annovar, annotate the locus
    results of expansionhunter_denovo using their shell script
    for this purpose.
    """
    input:
        annovar="results/annovar",
        annovar_humandb="results/annovar_humandb",
        expansionhunter_denovo="results/expansionhunter_denovo/repo",
        locus_results="results/expansionhunter_denovo/outlier_analysis/results.outlier_locus.tsv",
    output:
        tsv="results/expansionhunter_denovo/outlier_analysis/results.outlier_locus.annotated.tsv",
    conda:
        "../envs/annovar.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "bash {input.expansionhunter_denovo}/scripts/annotate_ehdn.sh "
        "--ehdn-results {input.locus_results} "
        "--ehdn-annotated-results {output.tsv} "
        "--annovar-annotate-variation {input.annovar}/annotate_variation.pl "
        "--annovar-humandb {input.annovar_humandb} "
        "--annovar-buildver hg38"


rule create_expansionhunter_denovo_report:
    """
    create report describing repeat distribution found in expansionhunter denovo runs
    """
    input:
        annotated_tsv="results/expansionhunter_denovo/outlier_analysis/results.outlier_locus.annotated.tsv",
    output:
        "results/expansionhunter_denovo/expansionhunter_denovo.html",
    params:
        locus_count=config["expansionhunter_denovo"]["plot-locus-count"],
        min_nonzero_subjects=lambda wildcards: len(reads_manifest) * 0.75,
    conda:
        "../envs/r.yaml"
    threads: config_resources["r"]["threads"]
    resources:
        mem_mb=config_resources["r"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["r"]["queue"], config_resources["queues"]
        ),
    script:
        "../scripts/expansionhunter_denovo.Rmd"
