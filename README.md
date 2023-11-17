# Metatranscriptomics-pipeline

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
this value will depend on available resources on execution infrastructure and the size of the dataset

`--cpus` [REQUIRED] number of cpus to be utilized at each stage of the assembly process. The default value is 8

`--dispersion` [REQUIRED] for differential expression analysis. Default is 0.1

`--bfly_HeapSpaceMax` [OPTIONAL] parameter for Butterlfy execution, default is 30G

`--bfly_HeapSpaceInit` [OPTIONAL] parameter for Butterfly execution, default is 5G

`--bfly_GCThreads` [OPTIONAL] specifies threads to be made availble for Butterfly execution, default is 10

`--bfly_CPU` [OPTIONAL] specifies maximum CPUs for Butterly algorithm, default is 6

`--bowtie2_thr` [OPTIONAL] defines the threads to be made available for bowtie2, default is 8

# Requirements for running pipeline

Nextflow 21+

Singularity 3.5+

# Input

The pipeline processes raw paired-end or single-end fastq reads. The pipeline takes a 4 column tab-delimited file (for paired-end reads) and a 3 column
tab-delimited file (for single-end reads). The tab-delimited specimen file specifies the sample ID, group name and absolute paths of raw reads.

# Output

A directory named <results> is saved in the current working directory. The following output is saved:
1. butterfly.txt - log file for last stage of assembly where the Butterfly algorithm does analysis.
  
2. chrysalis.txt - log file for assembly stage where Chrysalis algorithm does analysis.

3. DE - directory consisting of output from differential expression analysis. MA plots and Volcano plots comparing groups are saved here. The matrices used
to generate the plots are saved here.

4. inchworm.txt - log file to record Inchworm algorithm during assembly.
 
5. jellyfish.txt - log file for the first step of assembly.
 
6. trinity_matrix - directory consisting of the count matrix.

7. trinity_out_dir - directory consisting of assembly intermediate files and trimmed reads.

8. trinity_out_dir.Trinity.fasta - assembly output with all assembled transcript in fasta format.

9. trinity_out_dir.Trinity.fasta.gene_trans_map - gene mapping file.
 
10. trinity_quantification - directory consisting of isoform and gene counts.

# Software available in container

The container run Trinity (version 2.14.0). Development of the container and dependency software are specified in the definition file.
