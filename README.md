## Metatranscriptomics-pipeline

This repository hosts a pipeline for analysing metatranscriptomic data. This is a comprehensive pipeline launched using Nextflow and incorporates a Singularity
container. It performs de novo assembly and differential expression analysis allowing for better identification of novel genes or transcripts.

# Running the pipeline
The pipeline is intended to be used on a cluster or server. The simple command for running the pipeline is as follows:

`nextflow run trinity.nf -c config -profile cluster`

Additional options for the command are:

`-c` [REQUIRED] path to the configuration file
`-profile` [REQUIRED] defines the infrastructure being used for execution, cluster for running on cluster and server for running on local servers
`--reads` [REQUIRED] full path to the 4 column tab delimited specimen file
`--seqtype` [REQUIRED] defines the type of reads, can either be <fa> or <fq>, default is fq
`--libtype` [REQUIRED] defines the strand-specific RNA-seq read orientation. For pair-ended reads it can either be FR or RF, for single reads it can either
be F or R
`--mem` [OPTIONAL] maximum memory available for the algorthm to run. Input should be a string in the format "10G". Default value within the pipeline is 50G,
this value will depend on available resources on execution infrastructure and the size of the dataset.
`--cpus` [REQUIRED] number of cpus to be utilized at each stage of the assembly process. The default value is 8.
`--dispersion` [REQUIRED] for differential expression analysis. Default is 0.1.
`--bfly_HeapSpaceMax` [OPTIONAL] parameter for Butterlfy execution, default is 30G
`--bfly_HeapSpaceInit` [OPTIONAL] parameter for Butterfly execution, default is 5G
`--bfly_GCThreads` [OPTIONAL] specifies threads to be made availble for Butterfly execution, default is 10
`--bfly_CPU` [OPTIONAL] specifies maximum CPUs for Butterly algorithm, default is 6
`--bowtie2_thr` [OPTIONAL] defines the threads to be made available for bowtie2, default is 8
