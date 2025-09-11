#!/usr/bin/env nextflow

/* Elephant Analysis Pipeline version 0.3.0
Michael G. Campana, 2023-2025
Smithsonian\'s National Zoo and Conservation Biology Institute

The software is made available under the Smithsonian Institution terms of use (https://www.si.edu/termsofuse). */

gatk = 'gatk --java-options "' + params.java_options + '" ' // Simplify gatk command line
picard = "picard " + params.java_options
circulargenerator = "circulargenerator " + params.java_options
realignsamfile = "realignsamfile " + params.java_options

process prepareRef {
	
	// Prepare nuclear reference sequence
	
	input:
	path refseq
	
	output:
	path "${refseq.baseName}*.{amb,ann,bwt,pac,sa,fai,dict}"

	"""
	bwa index ${refseq}
	samtools faidx ${refseq}
	samtools dict ${refseq} > ${refseq.baseName}.dict
	"""

}

process prepareMitoRef {
	
	// Prepare mitochondrial reference sequence for BWA/CircularMapper
	
	input:
	path mtDNA
	val mtDNA_ID
	
	output:
	path "${mtDNA.baseName}*{_500.fasta,amb,ann,bwt,pac,sa,fai,dict,elongated}"
	
	"""
	$circulargenerator -e 500 -i ${mtDNA} -s ${mtDNA_ID}
	bwa index ${mtDNA.baseName}_500.fasta
	samtools faidx ${mtDNA}
	samtools dict ${mtDNA} > ${mtDNA.baseName}.dict
	"""

}

process trimReads {

	// Trim reads using AdapterRemoval v2
	
	input:
	tuple val(sample), val(library), path(reads1), path(reads2), val(rg), val(adapter1), val(adapter2)
	val(trimparams)
	
	output:
	tuple val(sample), val(library), path("${library}_R1.trunc.fastq.gz"), path("${library}_R2.trunc.fastq.gz"), val(rg)
	
	"""
	AdapterRemoval --file1 $reads1 --file2 $reads2 --basename $library --adapter1 $adapter1 --adapter2 $adapter2 --gzip $trimparams
	mv ${library}.pair1.truncated.gz ${library}_R1.trunc.fastq.gz
	mv ${library}.pair2.truncated.gz ${library}_R2.trunc.fastq.gz
	"""

}

process alignSeqs {

	// Align sequences using BWA
	
	input:
	tuple val(sample), val(library), path(reads1), path(reads2), val(rg)
	path refseq
	path "*"
	
	output:
	path("${library}_vs_genome.bam"), emit: bam
	val(sample), emit: sample
	val(library), emit: library
	val(rg), emit: rg
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa mem -t ${task.cpus} -R '${rg}' ${refseq} ${reads1} ${reads2} | samtools fixmate -@ ${samtools_extra_threads} -m - - | samtools sort -@ ${samtools_extra_threads} -o ${library}_vs_genome.bam -
	"""
	
}

process alignMitoSeqs {

	// Convert unmapped reads to FASTQ for alignment to mtDNA and align mitochondrial sequences using BWA
	
	input:
	path(bam)
	val(sample)
	val(library)
	val(rg)
	path mtDNA
	path "*"
	
	output:
	path("${library}_vs_mt.bam"), emit: bam
	val(sample), emit: sample
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools view -@ ${samtools_extra_threads} -b -f 4 ${bam} | samtools collate -@ ${samtools_extra_threads} -u -O - | samtools fastq -@ ${samtools_extra_threads} -1 ${library}.1.unmapped.fastq.gz -2 ${library}.2.unmapped.fastq.gz -0 /dev/null -s /dev/null
	bwa aln -t ${task.cpus} ${mtDNA.baseName}_500.fasta ${library}.1.unmapped.fastq.gz > ${library}.1.circ.sai
	bwa aln -t ${task.cpus} ${mtDNA.baseName}_500.fasta ${library}.2.unmapped.fastq.gz > ${library}.2.circ.sai
	bwa sampe -r '${rg}' ${mtDNA.baseName}_500.fasta ${library}.1.circ.sai ${library}.2.circ.sai ${library}.1.unmapped.fastq.gz ${library}.2.unmapped.fastq.gz > ${library}.circ.sam
	$realignsamfile -e 500 -i ${library}.circ.sam -r ${mtDNA}
	samtools fixmate -@ ${samtools_extra_threads} -m ${library}.circ_realigned.bam - | samtools sort -@ ${samtools_extra_threads} -o ${library}_vs_mt.bam -
	
	"""
}

process leftAlignIndels {

	// Left Align Indels using GATK4 LeftAlignIndels
	
	input:
	path(rg_bam)
	val(sample)
	path refseq
	path "*"
	
	output:
	tuple path("${rg_bam.simpleName}.realn.bam"), val(sample)
	
	script:
	if ( params.csi )
		"""
		samtools index -c ${rg_bam}
		$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam --create-output-bam-index false
		"""
	else
		"""
		$picard BuildBamIndex I=${rg_bam}
		$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam
		"""

}

process markDuplicates {

	// Initial marking of duplicates for unmerged library files
	
	publishDir "$params.outdir/01_LibraryBAMs", mode: 'copy'
	
	input:
	tuple path(lalnbam), val(sample)
	
	output:
	tuple path("${lalnbam.simpleName}.markdup.bam"), val(sample)
	
	script:
	samtools_extra_threads = task.cpus - 1
	if ( params.markDuplicates == "sambamba" )
		"""
		sambamba markdup ${lalnbam} ${lalnbam.simpleName}.markdup.bam
		"""
	else if ( params.markDuplicates == "samtools" )
		"""
		samtools markdup -@ ${samtools_extra_threads} ${lalnbam} ${lalnbam.simpleName}.markdup.bam
		"""
	else
		"""
		$picard MarkDuplicates I=${lalnbam} O=${lalnbam.simpleName}.markdup.bam M=${lalnbam.simpleName}.markdup.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		"""
	
}

process flagStats {

	// Calculate alignment statistics for unmerged library files using SAMtools flagstat
	
	publishDir "$params.outdir/02_LibraryFlagStats", mode: 'copy', pattern: '*.stats.txt'
	
	input:
	tuple path(mrkdupbam), val(sample)
	val(minmapped)
	
	output:
	path("${mrkdupbam.simpleName}.stats.txt")
	tuple path("${mrkdupbam.simpleName}.ok.bam"), val(sample), optional: true, emit: bam
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	#!/usr/bin/env bash
	samtools flagstat -@ ${samtools_extra_threads} ${mrkdupbam} > ${mrkdupbam.simpleName}.stats.txt
	primary=`sed -n \'2p\' ${mrkdupbam.simpleName}.stats.txt | cut -f 1 -d \" \"` # Primary alignments
	dup=`sed -n \'6p\' ${mrkdupbam.simpleName}.stats.txt | cut -f 1 -d \" \"` # Primary duplicates
	let total=\$primary-\$dup
	if [[ \$total -ge $minmapped ]]; then ln -s $mrkdupbam ${mrkdupbam.simpleName}.ok.bam; fi
	"""
	
}

process mergeSampleBAM {

	// Merge libraries by their sample IDs using SAMtools merge
	
	publishDir "$params.outdir/03_FinalBAMs", mode: 'copy', pattern: "*_merged_vs_*.markdup.bam"
	
	input:
	tuple path(bam), val(sample)
	val(marker)
	
	output:
	path "${sample}_vs_${marker}_merged.bam", emit: bam
	val(sample), emit: sample
	
	script:
	samtools_extra_threads = task.cpus - 1
	bams = 0
	bamlist = ""
	for (i in bam) {
		bams++
		bamlist = bamlist + " " + i
	}
	if (bams == 1) // Skip merging single libraries
		"""
		ln -s $bamlist ${sample}_vs_${marker}_merged.bam
		"""
	else
		"""
		samtools merge -@ ${samtools_extra_threads} -o ${sample}_vs_${marker}_merged.bam $bamlist
		"""
} 

process mergedMarkDup {

	// Mark duplicates for merged libraries after merging using Picard MarkDuplicates
	
	publishDir "$params.outdir/03_FinalBAMs", mode: 'copy'
	
	input:
	tuple path(lalnbam), val(sample)
	
	output:
	path "${sample}_vs_genome_merged.markdup.bam", optional: true, emit: genome
	path "${sample}_vs_mt_merged.markdup.bam", optional: true, emit: mt
	
	script:
	samtools_extra_threads = task.cpus - 1
	if ( params.markDuplicates == "sambamba" )
		"""
		sambamba markdup ${lalnbam} ${lalnbam.simpleName}.markdup.bam
		"""
	else if ( params.markDuplicates == "samtools" )
		"""
		samtools markdup -@ ${samtools_extra_threads} ${lalnbam} ${lalnbam.simpleName}.markdup.bam
		"""
	else
		"""
		$picard MarkDuplicates I=${lalnbam} O=${lalnbam.simpleName}.markdup.bam M=${lalnbam.simpleName}.markdup.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		"""
	
}

process mergedFlagStats {

	// Calculate alignment statistics using SAMtools flagstat
	
	publishDir "$params.outdir/04_FinalFlagStats", mode: 'copy'
	
	input:
	path mrkdupbam
	
	output:
	path "${mrkdupbam.simpleName}.stats.txt"
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools flagstat -@ ${samtools_extra_threads} ${mrkdupbam} > ${mrkdupbam.simpleName}.stats.txt
	"""

}

process callMtVariants {

	// Call mtDNA variants using GATK HaplotypeCaller
	
	publishDir "$params.outdir/05_IndividualgVCFs/mt", mode: 'copy'
	
	input:
	path final_bam
	path mtDNA
	path "*"
	
	output:
	path "${final_bam.simpleName}.vcf.gz"
	path "${final_bam.simpleName}.vcf.gz.tbi"
	
	"""
	samtools index $final_bam
	$gatk HaplotypeCaller -R ${mtDNA} -ploidy 1 -I $final_bam -O ${final_bam.simpleName}.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
	"""

}

process callGenomeVariants {

	// Call nuclear genome variants using GATK HaplotypeCaller
	
	publishDir "$params.outdir/05_IndividualgVCFs/genome", mode: 'copy'
	
	input:
	path final_bam
	path genome
	path "*"
	
	output:
	path "${final_bam.simpleName}.vcf.gz"
	path "${final_bam.simpleName}.vcf.gz.tbi"
	
	"""
	samtools index $final_bam
	$gatk HaplotypeCaller -R ${genome} -I $final_bam -O ${final_bam.simpleName}.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
	"""

}

process runPSMC {

	// Generate consensus sequence and run PSMC
	
	publishDir "$params.outdir/06_PSMC", mode: 'copy'
	
	input:
	path final_bam
	path genome
	path "*"
	val psmc_mpileup_opts
	val psmc_vcfutils_opts
	val psmc_psmcfa_opts
	val psmc_opts
	val psmc_bootstrap
	val psmc_plot_opts
	
	output:
	path "${final_bam.simpleName}.fq.gz"
	path "${final_bam.simpleName}*.psmcfa"
	path "${final_bam.simpleName}_bootstrapped.psmc"
	path "${final_bam.simpleName}_bootstrapped.eps"
	
	when:
	params.psmc == true
	
	"""
	#!/usr/bin/env bash
	bcftools mpileup $psmc_mpileup_opts -Ou --ignore-RG -f $genome $final_bam | bcftools call -c | vcfutils.pl vcf2fq $psmc_vcfutils_opts | gzip > ${final_bam.simpleName}.fq.gz
	fq2psmcfa $psmc_psmcfa_opts ${final_bam.simpleName}.fq.gz > ${final_bam.simpleName}.psmcfa
	splitfa ${final_bam.simpleName}.psmcfa > ${final_bam.simpleName}_split.psmcfa
	psmc $psmc_opts -o ${final_bam.simpleName}.psmc ${final_bam.simpleName}.psmcfa
	for i in {1..$params.psmc_bootstrap}; do psmc $psmc_opts -b -o ${final_bam.simpleName}_split\$i.psmc ${final_bam.simpleName}_split.psmcfa; done
	cat ${final_bam.simpleName}.psmc ${final_bam.simpleName}_split*.psmc > ${final_bam.simpleName}_bootstrapped.psmc
	psmc_plot.pl $psmc_plot_opts ${final_bam.simpleName}_bootstrapped ${final_bam.simpleName}_bootstrapped.psmc
	"""

}


workflow.onComplete {
	if (workflow.success) {
		println "Elephant pipeline completed successfully at $workflow.complete!"
		if (params.email != "NULL") {
			sendMail(to: params.email, subject: 'Elephant pipeline successful completion', body: "Elephant pipeline completed successfully at $workflow.complete!")
		}
	} else {
		println "Elephant pipeline terminated with errors at $workflow.complete.\nError message: $workflow.errorMessage"
		if (params.email != "NULL") {
			sendMail(to: params.email, subject: 'Elephant pipeline terminated with errors', body: "Elephant pipeline terminated with errors at $workflow.complete.\nError message: $workflow.errorMessage")
		}
	}
}

workflow merge_samples {
	// Left-align indels of merged data
	take:
		samples
		marker
		refseq
		refseq_files
	main:
		mergeSampleBAM(samples, marker)
		leftAlignIndels(mergeSampleBAM.out.bam, mergeSampleBAM.out.sample, refseq, refseq_files) | mergedMarkDup
	emit:
		mt = mergedMarkDup.out.mt
		genome = mergedMarkDup.out.genome
}

workflow mtDNA_processing {
	// Left-align indels, merge and mark duplicates for mtDNA BAMs.
	take:
		alignments
		sample
		library
		rg
	main:
		prepareMitoRef(params.mtDNA, params.mtDNA_ID)
		alignMitoSeqs(alignments, sample, library, rg, params.mtDNA, prepareMitoRef.out)
		leftAlignIndels(alignMitoSeqs.out.bam, alignMitoSeqs.out.sample, params.mtDNA, prepareMitoRef.out) | markDuplicates
		flagStats(markDuplicates.out, params.min_uniq_mapped)
		merge_samples(flagStats.out.bam.groupTuple(by: 1), "mt", params.mtDNA, prepareMitoRef.out)
		mergedFlagStats(merge_samples.out.mt)
		if (params.gatk) { callMtVariants(merge_samples.out.mt, params.mtDNA, prepareMitoRef.out) }
	emit:
		final_bams = final_mt_bams
		
}

workflow {
	main:
		prepareRef(params.refseq)
		if (params.read_trimming) {
			read_data = Channel.fromPath(params.samples).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.reads + row.Read1), file(params.reads + row.Read2), '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA', row.Adapter1, row.Adapter2)}
			trimReads(read_data, params.trimparams)
			alignSeqs(trimReads.out, params.refseq, prepareRef.out)
		} else {
			read_data = Channel.fromPath(params.samples).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.reads + row.Read1), file(params.reads + row.Read2), '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA')}
			alignSeqs(read_data, params.refseq, prepareRef.out)
		}
		if (params.circular_mtDNA) { mtDNA_processing(alignSeqs.out.bam, alignSeqs.out.sample, alignSeqs.out.library, alignSeqs.out.rg) } 
		leftAlignIndels(alignSeqs.out.bam, alignSeqs.out.sample, params.refseq, prepareRef.out) | markDuplicates
		flagStats(markDuplicates.out, params.min_uniq_mapped)
		merge_samples(flagStats.out.bam.groupTuple(by: 1), "genome", params.refseq, prepareRef.out)
		mergedFlagStats(merge_samples.out.genome)
		if (params.gatk) { callGenomeVariants(merge_samples.out.genome, params.refseq, prepareRef.out) }
		if (params.psmc) { runPSMC(merge_samples.out.genome, params.refseq, prepareRef.out, params.psmc_mpileup_opts, params.psmc_vcfutils_opts, params.psmc_psmcfa_opts, params.psmc_opts, params.psmc_bootstrap, params.psmc_plot_opts) }
}
	