localrules:
    caller_scatter_tasks,


rule caller_scatter_tasks:
    """
    For scatter-gather: take input set of range annotations
    and prepare for individual caller runs
    """
    input:
        lambda wildcards: config["deeptrio"][reference_build]["calling-ranges"],
    output:
        "results/{caller}/split_ranges/{splitnum}.ssv",
    benchmark:
        "results/performance_benchmarks/{caller}/split_ranges/{splitnum}.tsv"
    shell:
        "cat $(awk 'NR == {wildcards.splitnum}' {input}) | tr '\\n' ' ' > {output}"


rule deeptrio_make_examples:
    """
    Run deeptrio make_examples in a hybrid
    embarrassingly parallel fashion.
    """
    input:
        bam="results/bams/{projectid}/{sampleid}.bam",
        bai="results/bams/{projectid}/{sampleid}.bai",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
        intervals="results/deeptrio/split_ranges/{splitnum}.ssv",
    output:
        temp(
            expand(
                "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
                suffix=["gz", "gz.example_info.json"],
            )
        ),
        temp(
            expand(
                "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_make_examples/{projectid}/{sampleid}.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        gvcf_string=expand(
            "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        tmpdir="/tmp",
    container:
        "docker://google/deepvariant:{}".format(config["deeptrio"]["docker-version"])
    threads: config_resources["deeptrio"]["threads"]
    resources:
        mem_mb=config_resources["deeptrio"]["make_examples_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
        tmpdir="/tmp",
    shell:
        "mkdir -p {params.tmpdir} && "
        "seq 0 $(({threads}-1)) | parallel -j{threads} --tmpdir {params.tmpdir} "
        "make_examples --mode calling "
        '--ref {input.fasta} --reads {input.bam} --regions "$(cat {input.intervals})" '
        "--examples {params.shard_string} --channels insert_size "
        "--gvcf {params.gvcf_string} "
        "--task {{}}"


rule deeptrio_call_variants:
    """
    Run deeptrio call_variants in an
    embarrassingly parallel fashion.
    """
    input:
        expand(
            "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.{suffix}",
            shardnum=[
                str(i).rjust(5, "0")
                for i in range(config_resources["deeptrio"]["threads"])
            ],
            shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            suffix=["gz", "gz.example_info.json"],
        ),
    output:
        gz=temp(
            "results/deeptrio/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz"
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_call_variants/{projectid}/{sampleid}.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        docker_model="/opt/models/wgs/model.ckpt",
    container:
        "docker://google/deepvariant:{}".format(config["deeptrio"]["docker-version"])
    threads: config_resources["deeptrio"]["threads"]
    resources:
        mem_mb=config_resources["deeptrio"]["call_variants_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
    shell:
        "call_variants "
        "--outfile {output.gz} "
        "--examples {params.shard_string} "
        '--checkpoint "{params.docker_model}"'


rule deeptrio_postprocess_variants:
    """
    Run deeptrio postprocess_variants in an
    embarrassingly parallel fashion.
    """
    input:
        gz="results/deeptrio/{projectid}/call_variants/{sampleid}.{splitnum}.tfrecord.gz",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        gvcf=expand(
            "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
            shardnum=[
                str(i).rjust(5, "0")
                for i in range(config_resources["deeptrio"]["threads"])
            ],
            shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
        ),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
    output:
        vcf=temp(
            "results/deeptrio/{projectid}/postprocess_variants/{sampleid}.{splitnum}.vcf.gz"
        ),
        gvcf=temp(
            "results/deeptrio/{projectid}/postprocess_variants/{sampleid}.{splitnum}.g.vcf.gz"
        ),
        tbi=temp(
            "results/deeptrio/{projectid}/postprocess_variants/{sampleid}.{splitnum}.vcf.gz.tbi"
        ),
        html="results/deeptrio/{projectid}/postprocess_variants/{sampleid}.{splitnum}.visual_report.html",
    params:
        gvcf_string=expand(
            "results/deeptrio/{{projectid}}/make_examples/{{sampleid}}.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_postprocess_variants/{projectid}/{sampleid}.{splitnum}.tsv"
    container:
        "docker://google/deepvariant:{}".format(config["deeptrio"]["docker-version"])
    threads: 1
    resources:
        mem_mb=config_resources["deeptrio"]["postprocess_variants_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
    shell:
        "postprocess_variants "
        "--ref {input.fasta} "
        "--infile {input.gz} "
        "--nonvariant_site_tfrecord_path {params.gvcf_string} "
        "--gvcf_outfile {output.gvcf} "
        "--outfile {output.vcf}"


rule deeptrio_combine_regions:
    """
    Combine per-region deeptrio vcfs
    into a single mega vcf.
    """
    input:
        expand(
            "results/deeptrio/{{projectid}}/postprocess_variants/{{sampleid}}.{splitnum}.vcf.gz",
            splitnum=[i + 1 for i in range(caller_num_intervals)],
        ),
    output:
        "results/deeptrio/{projectid}/{sampleid}.sorted.vcf.gz",
    benchmark:
        "results/performance_benchmarks/deeptrio_combine_regions/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bcftools concat --threads {threads} -O z -o {output} {input}"


use rule deeptrio_combine_regions as deeptrio_combine_gvcfs with:
    input:
        expand(
            "results/deeptrio/{{projectid}}/postprocess_variants/{{sampleid}}.{splitnum}.g.vcf.gz",
            splitnum=[i + 1 for i in range(caller_num_intervals)],
        ),
    output:
        "results/deeptrio/{projectid}/{sampleid}.sorted.g.vcf.gz",
    benchmark:
        "results/performance_benchmarks/deeptrio_combine_gvcfs/{projectid}/{sampleid}.tsv"


rule rtg_create_sdf:
    """
    Convert a fasta to an sdf format *folder* for rtg tools' particularities
    """
    input:
        fasta="reference_data/bwa/{genome}/ref.fasta",
    output:
        directory("results/{genome}/ref.fasta.sdf"),
    benchmark:
        "results/performance_benchmarks/create_sdf/{genome}.tsv"
    conda:
        "../envs/vcfeval.yaml"
    threads: 1
    resources:
        qname="small",
        mem_mb=16000,
    shell:
        "rtg RTG_MEM=12G format -f fasta -o {output} {input}"


use rule somalier_build_pedfile as rtg_create_cluster_pedigree with:
    output:
        ped="results/deeptrio/{subset}.ped",
        problems=temp("results/deeptrio/{subset}.discordant_annotations.tsv"),
    benchmark:
        "results/performance_benchmarks/rtg_create_cluster_pedigree/{subset}.tsv"


rule rtg_annotate_vcf:
    input:
        vcf="results/glnexus/{family_cluster}/merged_callset.filtered.regions.vcf.gz",
        ped="results/deeptrio/{family_cluster}.ped",
        sdf="results/{}/ref.fasta.sdf".format(reference_build),
    output:
        vcf="results/deeptrio/{family_cluster}.annotated.vcf.gz",
    conda:
        "../envs/vcfeval.yaml"
    threads: 1
    resources:
        qname="small",
        mem_mb=16000,
    shell:
        "rtg RTG_MEM=12G mendelian -i {input.vcf} -o {output.vcf} --pedigree {input.ped} -t {input.sdf}"


rule aggregate_deeptrio_output:
    """
    Dispatch "deeptrio" runs for any proband who also has
    at least one parent present
    """
    input:
        lambda wildcards: tc.get_probands_with_structure(checkpoints),
    output:
        temp("results/deeptrio/.deeptrio_calls_split"),
    shell:
        "touch {output}"
