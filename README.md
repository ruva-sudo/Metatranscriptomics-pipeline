## Metatranscriptomics-pipeline

This repository hosts a pipeline for analysing metatranscriptomic data. This is a comprehensive pipeline launched using Nextflow and incorporates a Singularity
container. It performs de novo assembly and differential expression analysis allowing for better identification of novel genes or transcripts.

# Running the pipeline
The pipeline is intended to be used on a cluster or server. The simple command for running the pipeline is as follows:
`nextflow run trinity.nf -c config -profile cluster`
