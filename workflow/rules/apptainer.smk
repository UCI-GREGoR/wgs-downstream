rule apptainer_deeptrio:
    """
    Because snakemake's internal container support doesn't
    mesh with the deeptrio image, pull the image and create
    a sif image with apptainer.
    """
    output:
        "results/apptainer_images/deepvariant_{}.sif".format(
            config["deeptrio"]["docker-version"]
        ),
    params:
        outdir="results/apptainer_images",
        image_version=config["deeptrio"]["docker-version"],
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "apptainer pull --dir {params.outdir} docker://google/deepvariant:{params.image_version}"
