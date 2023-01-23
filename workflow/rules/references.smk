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
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
        tmpdir="temp/",
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
        "{prefix}.vcf.gz.tbi",
    benchmark:
        "results/performance_benchmarks/index_vcf/{prefix}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "tabix -p vcf {input}"


rule adjust_fasta_formatting:
    """
    exclusively because of idiosyncrasies in tiddit>=3, the fasta description lines can only
    contain the ">" character as the first character of the description line. any other instances
    of that character will cause tiddit>=3 to crash with a very cryptic error about list indices
    """
    input:
        "reference_data/references/{}/ref.fasta".format(reference_build),
    output:
        "reference_data/{{aligner}}/{}/ref.fasta".format(reference_build),
    threads: 1
    resources:
        mem_mb="1000",
        qname="small",
    shell:
        "sed 's/>/_/g' {input} | sed 's/^_/>/' > {output}"
