# Snakemake workflow: WGS Downstream Analyses

This workflow is intended to be the R&D space for PMGRC WGS downstream (post-alignment or calling) analyses. Target analyses are as follows:

- Cross-flowcell QC with [somalier](https://github.com/brentp/somalier)
- [ExpansionHunter Denovo](https://github.com/Illumina/ExpansionHunterDenovo)
- Joint calling with [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- others TBD

New global targets should be added in `workflow/Snakefile`. Content in `workflow/Snakefile` and the snakefiles in `workflow/rules` should be specifically _rules_; python infrastructure should be composed as subroutines under `lib/` and constructed in such a manner as to be testable with [pytest](https://docs.pytest.org/en/7.2.x/). Rules can call embedded scripts (in python or R/Rmd) from `workflow/scripts`; again, these should be constructed to be testable with pytest or [testthat](https://testthat.r-lib.org/).

## Authors

* Lightning Auriga (@lightning.auriga)

## Usage

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone git@gitlab.com:lightning.auriga1/wgs-pipeline.git
```

Note that this requires local git ssh key configuration; see [here](https://docs.gitlab.com/ee/user/ssh.html) for instructions as required.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `cram_manifest.tsv` or `gvcf_manifest.tsv` to specify your sample setup.

The following settings are recognized in `config/config.yaml`.

|Configuration Setting|Description|
|---|---|
|`cram_manifest`|Relative path to aligned read file manifest|
|`gvcf_manifest`|Relative path to g.vcf file manifest|
|`sample-sex`|Relative path to plaintext linker between sample ID and self-reported sex (Female, Male, Unknown)|
|`data-model`|Local Excel spreadsheet clone of data model information from Google docs. This file is generated as part of the file upload process to AnVIL, and is used for affected status annotations with ExpansionHunterDenovo. This is a local cloned file for the same reasons as listed above|
|`multiqc-config`|Relative path to configuration settings for cross-flowcell multiQC report|
|`genome-build`|Requested genome reference build to use for this analysis run. This should match the tags used in the reference data blocks below|
|`behaviors`|User-configurable modifiers to how the pipeline will run|
||- `symlink-crams`: whether to copy (no) or symlink (yes) input crams into workspace. symlinking is faster and more memory-efficient, but less reproducible, as the upstream files may vanish leaving no way to regenerate your analysis from scratch.
|`references`|Human genome reference data applicable to multiple tools. An arbitrary number of groupings can be specified under this key, meant to be named after the relevant human reference build. `grch38` is minimally required|
||- `fasta`: human sequence fasta file|
||- `exclusion-regions-bed`: set of bed regions to be excluded from output SNV data. expected to be from https://doi.org/10.1038/s41598-019-45839-z|
|`somalier`|Reference data specific to somalier; split by genome build, as above|
||- `sites-vcf-gz`: for use with somalier extract: sites targets in appropriate genome build|
||- `kg-labels-tsv`: for use with somalier ancestry estimation: ID/ancestry linker for 1KG samples|
||- `kg-reference-data-tar-gz`: for use with somalier ancestry estimation: extract output for 1KG samples|
|`expansionhunter_denovo`|Reference data specific to ExpansionHunterDenovo|
||`affected-status`: two column tab-delimited file mapping sample IDs to affected status|
||`repo`: relative or absolute path to local clone of [ExpansionHunterDenovo GitHub repository](https://github.com/Illumina/ExpansionHunterDenovo). This is required due to idiosyncrasies in the `expansionhunterdenovo` bioconda package. This behavior may be modified at a later date|
|`glnexus`|Reference data and settings specific to glnexus|
||- `version`: version string of glnexus docker image to use for analysis. example: `1.4.1`|
||- `config`: glnexus configuration preset name. examples: `DeepVariant`, `DeepVariant_unfiltered`|
||- `calling-ranges`: path to file containing list of filenames of calling range bedfiles|
|`deeptrio`|Settings specific to DeepTrio|
||- `docker-version`: which docker tag to use when pulling the official DeepTrio docker image|
||- `number-shards`: how many shards to break calling into. needs to be at most the number of available threads in the submission queue|
||Reference data specific to DeepTrio. An arbitrary number of groupings can be specified under genome-specific keys, of the form `grch.*`. `grch38` is minimally required|
||- `calling-ranges`: path to file containing list of filenames of calling range bedfiles|

The following columns are expected in the aligned read (cram) manifest, by default at `config/cram_manifest.tsv`:

|Column Header|Description|
|---|---|
|`sampleid`|External ID for sample|
|`cram`|Absolute or relative path to sample alignment (post-BQSR) cram file|

The following columns are expected in the g.vcf manifest, by default at `config/gvcf_manifest.tsv`:

|Column Header|Description|
|---|---|
|`sampleid`|External ID for sample|
|`gvcf`|Absolute or relative path to sample gvcf file|

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --profile sge-profile --cluster-config config/cluster.yaml --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### Cluster profiles

Snakemake interfaces with job schedulers via _cluster profiles_. For running jobs on SGE, you can use
the cookiecutter template [here](https://github.com/Snakemake-Profiles/sge).



#### Runtime Behaviors

This pipeline contains an assortment of utilities. Execution via `snakemake -jN` with no other arguments will
unconditionally run somalier, ExpansionHunterDenovo, and glnexus
for a single combined joint call across all input gvcfs. I may at some point
add the option to disable some of the tools in userspace, but there's no need currently. If you want this behavior,
you can run the phony rule named after each tool (e.g. `somalier` or `expansionhunter_denovo`).
The actual files considered by the tools are pulled from the input manifest with certain restrictions:

- somalier considers all crams in the user manifest (e.g. `config/cram_manifest.tsv`).
- for ExpansionHunterDenovo, which expects some amount of case/control annotation, only samples that are present in both the
  user manifest (e.g. `config/cram_manifest.tsv`) _and_ the data model spreadsheet `participant` tab are considered. this means,
  among other things, that the flowcell control low-depth NA24385 samples are always excluded from analysis.
- glnexus pulls all gvcfs from the relevant user manifest (by default, `config/gvcf_manifest.tsv`).
- DeepTrio is still in initial development, but the draft behavior is that the workflow will detect probands with at
  least one parent present in the inputs, and use DeepTrio to regenerate gvcfs for those individuals. Those gvcfs will
  then be joint called together, excluding all other subjects, with glnexus, and the resulting vcf will be annotated
  for _de novo_ variation. how exactly this will interact with glnexus as the primary analysis target is TBD, largely
  based on how well it turns out DeepTrio performs.

#### Triggering Reruns

This workflow is expected to be run periodically from the same installation. At some point, if there's enough iteration and
interest, I may update this workflow to handle perfect updating automatically. In lieu of that, here's how you can make
sure the pipeline runs the tools again correctly when new input samples are added to the cram manifest:

- you can always delete the `results/` directory entirely. this will purge all analysis to date. this works perfectly,
  but obviously involves redundant computation
- for somalier, reruns are gated on the contents of the linker file the workflow generates from the input. to force the
  pipeline to reevaluate its inputs and rerun if necessary, delete the file `results/linker.tsv` and relaunch
- for ExpansionHunterDenovo, reruns are gated on the contents of the linker file the workflow generates as well as
  the manifest the workflow generates from the combined information from the user cram manifest and the data model. as such,
  to force the pipeline to reevaluate its inputs and rerun if necessary, delete the files `results/linker.tsv` and
  `results/expansionhunter_denovo/manifest.tsv` and relaunch
- both DeepTrio and glnexus joint calling are tracked via single tracking files installed under `results/`. as always,
  just removing the respective results subdirectories will do the job, but the compute can be rather costly if you don't
  actually require a full rerun. you can individually remove whichever targets you want rerun, plus the trackers
  `results/split_joint_calls/.joint_calls_split` or `results/deeptrio/.deeptrio_calls_split`, to force the DAG to be
  reevaluated. note that this behavior is extremely WIP and likely to change to something more sensible.

#### A Discussion about DeepTrio

Per the tool's documentation, notably [here](https://github.com/google/deepvariant/blob/r1.5/docs/deeptrio-details.md),
DeepTrio uses the following logic for `make_examples`, the first step of the DeepVariant processing chain:

- for female proband, all available parents are use in `make_examples` for all chromosomes (excluding mitochondrion)
- for male proband:
    - for autosomes, all available parents are used
    - for chrX and chrY pseudoautosomal regions, all available parents are used
    - for the chrX non-pseudoautosomal region, only the maternal sample is used
    - for the chrY non-pseudoautosomal region, only the paternal sample is used

The default configuration of calling regions for GRCh38 is compatible with this calling logic. If you want to create
a different set of calling regions, e.g. for a different genome build, the above DAG logic is controlled by detecting
`chrX`, `chrY`, and `non-PAR` in the calling region filenames as appropriate.

Note that the above logic does not currently have a solution for the situation of a male proband, non-PAR regions,
one present and one absent parent in dataset. This is flagged for update.

### Step 5: Investigate results

#### somalier

Somalier's results are processed by MultiQC and emitted as a report `results/multiqc/multiqc.cross-flowcell.html`.

#### ExpansionHunterDenovo

WIP

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@gitlab.com:lightning.auriga1/wgs-pipeline.git` or `git remote add -f upstream https://github.com/snakemake-workflows/wgs-pipeline.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/default workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/default config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository. This project follows git flow; feature branches off of dev are welcome.

1. [Clone](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html) the fork to your local system, to a different place than where you ran your analysis.
2. Check out a branch off of dev:
```
git fetch
git checkout dev
git checkout -b your-new-branch
```
3. Make whatever changes best please you to your feature branch.
4. Commit and push your changes to your branch.
5. Create a [merge request](https://docs.gitlab.com/ee/user/project/merge_requests/) against dev.

## Testing

Testing infrastructure for embedded python and R scripts is installed under `lib/` and `workflow/scripts/`. Additional testing
coverage for the Snakemake infrastructure itself should be added once the workflow is more mature ([see here](https://github.com/lightning-auriga/snakemake-unit-tests)).

### Python testing with `pytest`
The testing under `lib/` is currently functional. Partial testing exists for the builtin scripts under `workflow/scripts`: the new utilities
for this implementation are tested, but some code inherited from the legacy pipeline(s) is not yet covered. To run the tests, do the following (from top level):

```bash
mamba install pytest-cov
pytest --cov=lib --cov=workflow/scripts lib workflow/scripts
```


### R testing with `testthat`
The testing under `workflow/scripts` is currently functional. The tests can be run with the utility script `run_tests.R`:

```bash
Rscript ./run_tests.R
```

To execute the above command, the environment must have an instance of R with appropriate libraries:

```bash
mamba install -c bioconda -c conda-forge "r-base>=4" r-testthat r-covr r-r.utils r-desctools
```
