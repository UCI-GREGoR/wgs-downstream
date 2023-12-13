rule download_reference_data:
    """
    Conditionally acquire reference data files from either remote
    source (e.g. S3) or somewhere else on local filesystem.

    All the interesting stuff happens in the input mapping function,
    which uses the path contents between 'reference_data' and the actual
    filename to determine which configured reference file to pull.

    In some instances, the remote source will be gzipped but the local version
    won't be; in that instance, decompress midflight.

    I'm giving up on remotes, because the FTP one is unusable with conda env creation.
    """
    output:
        "reference_data/{reference_file}",
    params:
        lambda wildcards: tc.map_reference_file(wildcards, config),
    benchmark:
        "results/performance_benchmarks/download_reference_data/{reference_file}.tsv"
    conda:
        "../envs/awscli.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        'if [[ "{params}" == "s3://"* ]] ; then aws s3 cp {params} {output}.staging ; '
        'elif [[ "{params}" == "http://"* ]] || [[ "{params}" == "https://"* ]] || [[ "{params}" == "ftp://"* ]] ; then wget -O {output}.staging {params} ; '
        "else cp {params} {output}.staging ; fi ; "
        'if [[ "{params}" = *".gz" ]] && [[ "{output}" != *".gz" ]] ; then cat {output}.staging | gunzip -c > {output} && rm {output}.staging ; '
        "else mv {output}.staging {output} ; fi"


rule index_vcf:
    """
    Use tabix to generate tbi file for vcf.gz input
    """
    input:
        "{prefix}.vcf.gz",
    output:
        temp("{prefix}.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/index_vcf/{prefix}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "tabix -p vcf {input}"


localrules:
    adjust_fasta_formatting,


rule adjust_fasta_formatting:
    """
    Once upon a time, this was supposed to adjust formatting in the input fasta
    to force compatibility between certain versions of tiddit and a custom (inherited, janky)
    build of GRCh37. However, since GRCh37 is very much not supported, and if you want it
    you should use the public version anyway, this is getting replaced.
    """
    input:
        "reference_data/references/{}/ref.fasta".format(reference_build),
    output:
        "reference_data/{{aligner}}/{}/ref.fasta".format(reference_build),
    shell:
        "ln -s $(readlink -m {input}) {output}"


rule samtools_index_fasta:
    """
    Given a fasta, generate its fai index file.
    """
    input:
        "{prefix}fasta",
    output:
        "{prefix}fasta.fai",
    benchmark:
        "results/performance_benchmarks/samtools_index_fasta/{prefix}fasta.fai.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb=config_resources["samtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["samtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "samtools faidx {input}"
