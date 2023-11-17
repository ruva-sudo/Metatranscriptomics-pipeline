nextflow.enable.dsl = 2

params.reads = "/home/rmadzime/rna_preprocessing/import_files.txt"
params.outdir = "results"
params.mem = "50G"
params.cpus = 8
params.seqtype = 'fq'
params.libtype = 'FR'
params.dispersion = 0.1
params.bfly_HeapSpaceMax = "30G"
params.bfly_HeapSpaceInit = "5G"
params.bfly_GCThreads = 10
params.bfly_CPU = 6
params.bowtie2_thr = 8

process JELLYFISH{

	publishDir "${params.outdir}", mode: 'copy'

	input:
	path reads

	output:
	path "trinity_out_dir",				emit: jellyfish_dir
	path "trinity_out_dir/both.fa",			emit: kmer_counts
	path "jellyfish.txt",				emit: log
	path "trinity_out_dir/*_R1.fq.gz.P.qtrim.gz",	emit: forward_reads
	path "trinity_out_dir/*_R2.fq.gz.P.qtrim.gz",	emit: reverse_reads
	path "trinity_out_dir/*_R1.fq.gz.U.qtrim.gz",	emit: unpaired_1
	path "trinity_out_dir/*_R2.fq.gz.U.qtrim.gz",	emit: unpaired_2

	script:
	"""
	Trinity --seqType ${params.seqtype} \
		--max_memory ${params.mem} \
		--samples_file ${reads} \
		--trimmomatic \
		--SS_lib_type ${params.libtype} \
		--CPU ${task.cpus} \
		--no_run_inchworm > jellyfish.txt 2>&1
	"""
}

process INCHWORM{

	publishDir "${params.outdir}", mode: 'copy'

	input:
	path kmers
	path reads
	path trimmed

	output:
	path "trinity_out_dir", 			  emit: inch_dir
	path "trinity_out_dir/inchworm.fa", 		  emit: isocontigs
	path "trinity_out_dir/jellyfish.kmers.25.asm.fa", emit: kmers
	path "inchworm.txt",				  emit: log

	script:
	"""
	Trinity --seqType ${params.seqtype} \
		--max_memory ${params.mem} \
		--samples_file ${reads} \
		--trimmomatic \
		--SS_lib_type ${params.libtype} \
		--CPU ${task.cpus} \
		--no_run_chrysalis > inchworm.txt 2>&1
	"""
}

process CHRYSALIS{

	publishDir "${params.outdir}", mode: 'copy'

	input:
	path kmers
	path isocontigs
	path reads
	path trimmed

	output:
	path "trinity_out_dir", 			emit: chrysalis_dir
	path "trinity_out_dir/chrysalis/*", 		emit: sub_directory
	path "trinity_out_dir/recursive_trinity.cmds",  emit: graph
	path "trinity_out_dir/scaffolding_entries.sam", emit: scaffolds
	path "chrysalis.txt",				emit: log

	script:
	"""
	Trinity --seqType ${params.seqtype} \
		--max_memory ${params.mem} \
		--samples_file ${reads} \
		--trimmomatic \
		--SS_lib_type ${params.libtype} \
		--CPU ${task.cpus} \
		--no_distributed_trinity_exec > chrysalis.txt 2>&1
	"""
}

process BUTTERFLY{

	publishDir "${params.outdir}", mode: 'copy'

	input:
	path isoforms
	path scaffolds
	path cluster_and_graphs
	path reads
	path trimmed

	output:
	path "trinity_out_dir", 	       emit: full_directory
	path "*.Trinity.fasta", 	       emit: transcripts
	path "*.Trinity.fasta.gene_trans_map", emit: gene_map
	path "butterfly.txt",		       emit: log_file

	script:
	"""
	Trinity --seqType ${params.seqtype} \
		--max_memory ${params.mem} \
		--samples_file ${reads} \
		--trimmomatic \
		--SS_lib_type ${params.libtype} \
		--bflyHeapSpaceMax ${params.bfly_HeapSpaceMax} \
		--bflyHeapSpaceInit ${params.bfly_HeapSpaceInit} \
		--bflyGCThreads ${params.bfly_GCThreads} \
		--bflyCPU ${params.bfly_CPU} \
		--CPU ${task.cpus} \
		--no_version_check > butterfly.txt 2>&1
	"""
}

process BUILD_INDEX {

	publishDir "${params.outdir}/bowtie2_alignment",mode: 'copy'

	input:
	path transcripts

	output:
	path "trinity_out_dir*",	emit: index

	script:
	"""
	bowtie2-build ${transcripts} trinity_out_dir
	"""
}

process BOWTIE_ALIGN {

	publishDir "${params.outdir}/bowtie2_quality", mode: 'copy'

	input:
	path index
	path forward_reads
	path reverse_reads
	path unpaired_1
	path unpaired_2

	output:
	path "trinity_out_dir_bt2.bam",	emit: log_file
	path "bowtie2.out",		emit: output_file

	script:
	"""
	bowtie2 -q -1 ${forward_reads} -2 ${reverse_reads} \
	-U ${unpaired_1} -U ${unpaired_2} \
	-x trinity_out_dir -k 10 -p ${params.bowtie2_thr} \
	--no-unal -S bowtie2.out 2>> trinity_out_dir_bt2.bam
	"""
}
	
process QUANT{

	publishDir "${params.outdir}/trinity_quantification", mode: 'copy'

	input:
	path transcripts
	path reads

	output:
	path "*/*.isoforms.results",	emit: isoforms
	path "*/RSEM.genes.results",	emit: gene_counts

	script:
	"""
	align_and_estimate_abundance.pl --transcripts ${transcripts} \
					--samples_file ${reads} \
					--est_method RSEM \
					--seqType ${params.seqtype} \
					--SS_lib_type ${params.libtype} \
					--aln_method bowtie2 \
					--trinity_mode \
					--prep_reference

	prefix=\$(ls -l | grep -E ^d | awk '{print \$9}')
	for i in \${prefix[@]}
	do
		mv \${i}/RSEM.isoforms.results \${i}/\${i}.isoforms.results
	done
	"""
}


process MATRIX {

	publishDir "${params.outdir}/trinity_matrix", mode: 'copy'

	input:
	path isoform_counts
	path map

	output:
	path "*.isoform.counts.matrix",	emit: isoform_matrix
	path "*.gene.counts.matrix",	emit: gene_matrix

	script:
	"""
	abundance_estimates_to_matrix.pl --est_method RSEM \
					 --gene_trans_map ${map} \
					 --out_prefix RSEM \
					 ${isoform_counts}
	"""
}

process DE_ANALYSIS {

	publishDir "${params.outdir}/DE", mode: 'copy'

	input:
	path matrix_counts

	output:
	path "edgeR/*.DE_results",			emit: DE_results
	path "edgeR/*.DE_results.MA_n_Volcano.pdf",	emit: DE_plots
	path "edgeR/*.count_matrix",			emit: count_matrix
	path "edgeR/*.EdgeR.Rscript",			emit: Rscript

	script:
	"""
	run_DE_analysis.pl --matrix ${matrix_counts} \
			   --output edgeR \
			   --method edgeR \
			   --dispersion ${params.dispersion}
	"""
}
 

workflow{

	reads_ch = Channel.fromPath(params.reads)
	//reads_ch.view()

	JELLYFISH(reads_ch)
	INCHWORM(JELLYFISH.out.kmer_counts, reads_ch, JELLYFISH.out.jellyfish_dir)
	CHRYSALIS(INCHWORM.out.isocontigs, INCHWORM.out.kmers, reads_ch, INCHWORM.out.inch_dir)
	BUTTERFLY(CHRYSALIS.out.graph, CHRYSALIS.out.scaffolds, CHRYSALIS.out.sub_directory, reads_ch, CHRYSALIS.out.chrysalis_dir)
	BUILD_INDEX(BUTTERFLY.out.transcripts)
	forward_ch = JELLYFISH.out.forward_reads.collect()
	reverse_ch = JELLYFISH.out.reverse_reads.collect()
	U1_ch = JELLYFISH.out.unpaired_1.collect()
	U2_ch = JELLYFISH.out.unpaired_2.collect()
	BOWTIE_ALIGN(BUILD_INDEX.out.index, forward_ch, reverse_ch, U1_ch, U2_ch)
	QUANT(BUTTERFLY.out.transcripts, reads_ch)
	isoform_ch = QUANT.out.isoforms.collect()
	//isoform_ch.view()
	MATRIX(isoform_ch, BUTTERFLY.out.gene_map)
	DE_ANALYSIS(MATRIX.out.isoform_matrix)
}
