#!/usr/bin/env nextflow

// Forked from Kauai pipeline on 20 Jan 2023
// Optimizations taken from Drep_Phylogenomics pipeline dated 20 Jan 2023

gatk = 'gatk --java-options "' + params.java_options + '" ' // Simplify gatk command line

Channel
	.fromPath(params.samples)
	.splitCsv(header:true)
	.map { row -> tuple(row.Sample, row.Library, file(params.reads + row.Read1), file(params.reads + row.Read2), '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA') }
	.set { readpairs_ch }

process buildRef {
	
	// Prepare nuclear reference sequence
	
	input:
	path refseq from params.refseq
	
	output:
	path "${refseq.baseName}.*" into ref_build_ch // Channel to find these files again
	path "${refseq.baseName}*.{fai,dict}" into fai_refseq_laln_ch // Channel for fai file for LeftAlignIndels

	"""
	bwa index ${refseq}
	samtools faidx ${refseq}
	samtools dict ${refseq} > ${refseq.baseName}.dict
	"""

}

process buildMitoRef {
	
	// Prepare mitochondrial reference sequence for BWA
	
	input:
	path mtDNA from params.mtDNA
	
	output:
	path "${mtDNA.baseName}.*" into ref_mtDNA_ch // Channel to find alignment files again
	path "${mtDNA.baseName}*.{fai,dict}" into fai_mtDNA_laln_ch // Channel for fai file for LeftAlignIndels
	
	"""
	bwa index ${mtDNA}
	samtools faidx ${mtDNA}
	samtools dict ${mtDNA} > ${mtDNA.baseName}.dict
	"""

}

process alignSeqs {

	// Align sequences using BWA and convert unmapped reads to FASTQ for alignment to mtDNA
	
	input:
	path refseq from params.refseq
	path "*" from ref_build_ch
	tuple val(sample), val(library), path(reads1), path(reads2), val(rg) from readpairs_ch
	
	output:
	tuple file("${library}_vs_genome.bam"), val(sample) into raw_bam_ch
	tuple val(library), path("${library}.1.unmapped.fastq.gz"), path("${library}.2.unmapped.fastq.gz"), val(sample), val(rg) into mtDNA_fastq_ch

	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa mem -t ${task.cpus} -R '${rg}' ${refseq} ${reads1} ${reads2} | samtools view -@ ${samtools_extra_threads} -b -o ${library}.bam -
	samtools view -@ ${samtools_extra_threads} -b -F 4 ${library}.bam | samtools fixmate -@ ${samtools_extra_threads} -r -m - - | samtools sort -@ ${samtools_extra_threads} -o ${library}_vs_genome.bam -
	samtools view -@ ${samtools_extra_threads} -b -f 4 ${library}.bam | samtools sort -@ ${samtools_extra_threads} -o ${library}.unmapped.bam -
	samtools collate -@ ${samtools_extra_threads} -u -O ${library}.unmapped.bam | \\
	samtools fastq -@ ${samtools_extra_threads} -1 ${library}.1.unmapped.fastq -2 ${library}.2.unmapped.fastq -0 /dev/null -s /dev/null
	gzip ${library}.*.unmapped.fastq
	"""
	
}

process alignMitoSeqs {

	// Align mitochondrial sequences using BWA
	
	input:
	tuple val(library), path(mtfastq1), path(mtfastq2), val(sample), val(rg) from mtDNA_fastq_ch
	file "*" from ref_mtDNA_ch
	path mtDNA from params.mtDNA
	
	output:
	tuple path("${library}_vs_mt.bam"), val(sample) into raw_mito_bam_ch
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa mem -t ${task.cpus} -R '${rg}' ${mtDNA} ${mtfastq1} ${mtfastq2} | samtools view -@ ${samtools_extra_threads} -b -F 4 - | samtools fixmate -@ ${samtools_extra_threads} -r -m - - | samtools sort -@ ${samtools_extra_threads} -o ${library}_vs_mt.bam -
	"""
}

raw_bam_ch2 = raw_bam_ch.mix(raw_mito_bam_ch)

process leftAlignIndels {

	// Left Align Indels using GATK4 LeftAlignIndels
	
	input:
	tuple path(rgbam), val(sample) from raw_bam_ch2
	path mtDNA from params.mtDNA
	path mtDNA_fai from fai_mtDNA_laln_ch
	path genome from params.refseq
	path genome_fai from fai_refseq_laln_ch
	
	output:
	tuple path("${rgbam.simpleName}.laln.bam"), val(sample) into laln_bam_ch
	
	script:
	 // Need to identify appropriate reference sequence for alignment
	refid = rgbam.simpleName.split('_vs_')[1]
	switch(refid) {
		case 'mt':
			laln_reference = mtDNA;
			break;
		case 'genome':
			laln_reference = genome;
			break;
	}
	"""
	$gatk LeftAlignIndels -I ${rgbam} -O ${rgbam.simpleName}.laln.bam -R ${laln_reference}
	"""
	
}	

process markDup {

	// Initial marking of duplicates for unmerged library files
	
	publishDir "$params.outdir/01_LibraryBAMs", mode: 'copy'
	
	input:
	tuple path(lalnbam), val(sample) from laln_bam_ch
	val java_options from params.java_options
	
	output:
	tuple path("${lalnbam.simpleName}.mrkdup.bam"), val(sample) into mrkdup_bam_ch
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools markdup -@ ${samtools_extra_threads} ${lalnbam} ${lalnbam.simpleName}.mrkdup.bam
	"""
	
}

process flagStats {

	// Calculate alignment statistics for unmerged library files using SAMtools flagstat
	
	publishDir "$params.outdir/02_LibraryFlagStats", mode: 'copy', pattern: '*.stats.txt'
	
	input:
	tuple path(mrkdupbam), val(sample) from mrkdup_bam_ch
	val(minmapped) from params.min_uniq_mapped
	
	output:
	path("${mrkdupbam.simpleName}.stats.txt")
	tuple path("${mrkdupbam.simpleName}.trim.bam"), val(sample) optional true into modern_bam_ch
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools flagstat -@ ${samtools_extra_threads} ${mrkdupbam} > ${mrkdupbam.simpleName}.stats.txt
	primary=`sed -n \'2p\' ${mrkdupbam.simpleName}.stats.txt | cut -f 1 -d \" \"` # Primary alignments
	dup=`sed -n \'6p\' ${mrkdupbam.simpleName}.stats.txt | cut -f 1 -d \" \"` # Primary duplicates
	let total=\$primary-\$dup
	if [[ \$total -ge $minmapped ]]; then ln $mrkdupbam ${mrkdupbam.simpleName}.trim.bam; fi
	"""
	
}

sample_bam_ch = modern_bam_ch.groupTuple(by: 2) // Get a channel of unique samples matched with their file paths

process mergeSampleBAM {

	// Merge libraries by their sample IDs using SAMtools merge
	
	input:
	tuple path(bam), val(sample) from sample_bam_ch
	
	output:
	path "${sample}_merged*.bam" // Make sure there is some output
	path "${sample}_merged*.merged.bam" optional true into merged_bam_ch // Send samples that need merging to merging processes
	path "${sample}_merged_vs_genome.mrkdup.bam" optional true into final_bam_skip_ch // Skip unnecessary merging steps
	path "${sample}_merged_vs_mt.mrkdup.bam" optional true into final_mt_skip_ch // Skip unnecessary merging steps
	
	script:
	samtools_extra_threads = task.cpus -1 
	// First make sure that an input file exists for each type of alignment since could have been removed earlier. If no alignments exist, the whole sample should have been removed previously.
	mtbamlist = ""
	genomebamlist = ""
	mtbams = 0 // Count of mt bams
	genomebams = 0 // Count of genome bams
	for (i in bam) {
		category = i.simpleName.split("_vs_")[1]
		if (category == "mt") {
			mtbams++
			mtbamlist = mtbamlist + " " + i
		} else {
			genomebams++
			genomebamlist = genomebamlist + " " + i
		}
	}
	if (genomebams == 0 && mtbams == 1)
		"""
		ln -s $mtbamlist ${sample}_merged_vs_mt.mrkdup.bam
		"""
	 else if (genomebams == 0 && mtbams > 1)
		"""
		samtools merge -@ ${samtools_extra_threads} ${sample}_merged_vs_mt.merged.bam $mtbamlist
		"""
	else if (mtbams == 0 && genomebams == 1)
		"""
		ln -s $genomebamlist ${sample}_merged_vs_genome.mrkdup.bam
		"""
	else if (mtbams == 0 && genomebams > 1)
		"""	
		samtools merge -@ ${samtools_extra_threads} ${sample}_merged_vs_genome.merged.bam $genomebamlist
		"""
	else if (mtbams == 1 && genomebams > 1) 
		"""
		ln -s $mtbamlist ${sample}_merged_vs_mt.mrkdup.bam
		samtools merge -@ ${samtools_extra_threads} ${sample}_merged_vs_genome.merged.bam $genomebamlist
		"""
	else if (mtbams > 1 && genomebams == 1)
		"""
		ln -s $genomebamlist ${sample}_merged_vs_genome.mrkdup.bam
		samtools merge -@ ${samtools_extra_threads} ${sample}_merged_vs_mt.merged.bam $mtbamlist
		"""
	else if (mtbams == 1 && genomebams == 1)
		"""
		ln -s $mtbamlist ${sample}_merged_vs_mt.mrkdup.bam
		ln -s $genomebamlist ${sample}_merged_vs_genome.mrkdup.bam
		"""
	else
		"""
		samtools merge -@ ${samtools_extra_threads} ${sample}_merged_vs_mt.merged.bam $mtbamlist
		samtools merge -@ ${samtools_extra_threads} ${sample}_merged_vs_genome.merged.bam $genomebamlist
		"""
} 

// Flatten array if both mt and genome samples need merging
merged_bam_ch2 = merged_bam_ch.flatten()


process mergedLeftAlignIndels {

	// Left align indels for merged libraries
	
	input:
	path mrgbam from merged_bam_ch2
	path mtDNA from params.mtDNA
	path mtDNA_fai from fai_mtDNA_laln_ch
	path genome from params.refseq
	path genome_fai from fai_refseq_laln_ch
	
	output:
	file "${mrgbam.simpleName}.laln.bam" into laln_merged_bam_ch
	
	script:
	 // Need to identify appropriate reference sequence for alignment
	refid = mrgbam.simpleName.split('_vs_')[1]
	switch(refid) {
		case 'mt':
			laln_reference = mtDNA;
			break;
		case 'genome':
			laln_reference = genome;
			break;
	}
	"""
	$gatk LeftAlignIndels -I ${mrgbam} -O ${mrgbam.simpleName}.laln.bam -R ${laln_reference}
	"""

}

process mergedMarkDup {

	// Mark duplicates for merged libraries after merging using Picard MarkDuplicates
	
	publishDir "$params.outdir/03_FinalBAMs", mode: 'copy'
	
	input:
	path laln_mrg_bam from laln_merged_bam_ch
	val java_options from params.java_options
	
	output:
	path "${laln_mrg_bam.simpleName}.mrkdup.bam" into mrg_mrkdup_bam_ch
	path "${laln_mrg_bam.simpleName.split('_vs_')[0]}_vs_mt.mrkdup.bam" optional true into final_mt_ch
	path "${laln_mrg_bam.simpleName.split('_vs_')[0]}_vs_genome.mrkdup.bam" optional true into final_bam_ch
	
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools markdup -@ ${samtools_extra_threads} ${laln_mrg_bam} ${laln_mrg_bam.simpleName}.mrkdup.bam
	"""

}

process mergedFlagStats {

	// Calculate alignment statistics using SAMtools flagstat
	
	publishDir "$params.outdir/04_FinalFlagStats", mode: 'copy'
	
	input:
	path mrkdupbam from mrg_mrkdup_bam_ch
	
	output:
	path "${mrkdupbam.simpleName}.stats.txt"
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools flagstat -@ ${samtools_extra_threads} ${mrkdupbam} > ${mrkdupbam.simpleName}.stats.txt
	"""

}

// Add different channels together
final_mt_ch2 = final_mt_ch.mix(final_mt_skip_ch)
final_bam_ch2 = final_bam_ch.mix(final_bam_skip_ch)


process callMtVariants {

	// Call mtDNA variants using GATK HaplotypeCaller
	
	publishDir "$params.outdir/05_IndividualgVCFs/mt", mode: 'copy'
	
	input:
	path final_bam from final_mt_ch2.unique()
	path mtDNA from params.mtDNA
	path mtDNA_fai from fai_mtDNA_laln_ch
	
	output:
	path "${final_bam.simpleName}.vcf.gz" into gVCF_mt_ch
	path "${final_bam.simpleName}.vcf.gz.tbi" into gVCF_mt_index_ch
	
	"""
	samtools index $final_bam
	$gatk HaplotypeCaller -R ${mtDNA} -ploidy 1 -I $final_bam -O ${final_bam.simpleName}.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
	"""

}

process callGenomeVariants {

	// Call nuclear genome variants using GATK HaplotypeCaller
	
	publishDir "$params.outdir/05_IndividualgVCFs/genome", mode: 'copy'
	
	input:
	path final_bam from final_bam_ch2.unique()
	path genome from params.refseq
	path genome_fai from fai_refseq_laln_ch
	
	output:
	path "${final_bam.simpleName}.vcf.gz" into gVCF_genome_ch
	path "${final_bam.simpleName}.vcf.gz.tbi" into gVCF_genome_index_ch
	
	"""
	samtools index $final_bam
	$gatk HaplotypeCaller -R ${genome} -I $final_bam -O ${final_bam.simpleName}.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
	"""

}