
**Finding miRNA targets**.

[![GitHub Actions CI Status](https://github.com/payamemami/mirnatargetdetector/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/mirnatargetdetector/actions)
[![GitHub Actions Linting Status](https://github.com/payamemami/mirnatargetdetector/workflows/nf-core%20linting/badge.svg)](https://github.com/payamemami/mirnatargetdetector/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/payamemami/mirnatargetdetector.svg)](https://hub.docker.com/r/nfcore/mirnatargetdetector)

## Introduction

**payamemami/mirnatargetdetector** is a bioinformatics best-practise analysis pipeline for

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run payamemami/mirnatargetdetector -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run payamemami/mirnatargetdetector --input 'path to miRNA fasta file' --input_utr 'path to UTR fasta file' -profile docker
    ```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:

* miranda target detection
* RNAHybrid target detection

## Documentation

The nf-core/mirnatargetdetector pipeline comes with documentation about the pipeline: [usage](docs/usage.md) and [output](docs/output.md).

## Credits

nf-core/mirnatargetdetector was originally written by Payam Emami.

We thank the following people for their extensive assistance in the development
of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).
