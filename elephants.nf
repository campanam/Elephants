#!/usr/bin/env nextflow

// Forked from Kauai pipeline on 20 Jan 2023
// Optimizations taken from Drep_Phylogenomics pipeline dated 20 Jan 2023

gatk = 'gatk --java-options "' + params.java_options + '" ' // Simplify gatk command line

Channel
	.fromPath(params.samples)
	.splitCsv(header:true)
	.map { row -> tuple(row.Sample, row.Library, file(params.reads + row.Read1), file(params.reads + row.Read2), '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA',row.Adapter1, row.Adapter2) }
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
	val mtDNA_ID from params.mtDNA_ID
	
	output:
	path "${mtDNA.baseName}.*" into ref_mtDNA_ch // Channel to find alignment files again
	path "${mtDNA.baseName}*.{fai,dict}" into fai_mtDNA_laln_ch // Channel for fai file for LeftAlignIndels
	
	"""
	${params.bin}bwa index ${mtDNA}
	${params.bin}samtools faidx ${mtDNA}
	${params.bin}samtools dict ${mtDNA} > ${mtDNA.baseName}.dict
	"""

}

process trimAdapters {

	// Trim adapters using AdapterRemoval 2.3.1

	input:
	tuple val(sample), val(library), path(reads1), path(reads2), val(rg), val(adapter1), val(adapter2) from readpairs_ch
	
	output:
	tuple val(library), path("${library}.R1.fastq.gz"), path("${library}.R2.fastq.gz"), val(sample), val(rg) into trim_readpairs_ch
	
	"""
	AdapterRemoval --file1 $reads1 --file2 $reads2 --basename $library --adapter1 $adapter1 --adapter2 $adapter2 --gzip --minlength 30
	mv ${library}.pair1.truncated.gz ${library}.R1.fastq.gz
	mv ${library}.pair2.truncated.gz ${library}.R2.fastq.gz
	"""

}

process alignSeqs {

	// Align sequences using BWA and convert unmapped reads to FASTQ for alignment to mtDNA
	
	input:
	path refseq from params.refseq
	path "*" from ref_build_ch
	tuple val(library), path(reads1), path(reads2), val(sample), val(rg) from trim_readpairs_ch
	
	output:
	tuple file("${library}_vs_genome.bam"), val(sample) into raw_bam_ch
	tuple val(library), file("${library}.1.unmapped.fastq.gz"), file("${library}.2.unmapped.fastq.gz"), val(sample), val(rg) into mtDNA_fastq_ch

	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa mem -t ${task.cpus} -R ${rg} ${refseq} ${reads1} ${reads2} | samtools view -@ ${samtools_extra_threads} -b -o ${library}.bam -
	samtools view -@ ${samtools_extra_threads} -b -F 4 ${library}.bam - | samtools fixmate -@ ${samtools_extra_threads} -r -m - - | samtools sort -@ ${samtools_extra_threads} -o ${library}_vs_genome.bam -
	samtools view -@ ${samtools_extra_threads} -b -f 4 ${library}.bam | samtools sort -@ ${samtools_extra_threads} -o ${library}.unmapped.bam -
	samtools collate -@ ${samtools_extra_threads} -u -O ${library}.unmapped.bam | \\
	samtools fastq -@ ${samtools_extra_threads} -1 ${library}.1.unmapped.fastq -2 ${library}.2.unmapped.fastq -0 /dev/null -s /dev/null
	gzip ${library}.*.unmapped.fastq
	"""
	
}

process alignMitoSeqs {

	// Align mitochondrial sequences using BWA
	
	input:
	tuple val(library), path(mtfastq1), path(mtfastq2), val(sample), val(species), val(rg) from mtDNA_fastq_ch
	file "*" from ref_mtDNA_ch
	path mtDNA from params.mtDNA
	
	output:
	tuple path("${library}_vs_mt.bam"), val(sample) into raw_mito_bam_ch
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa mem -t ${task.cpus} -R ${rg} ${mtDNA} ${mtfastq1} ${mtfastq2} | samtools view -@ ${samtools_extra_threads} -b -F 4 - | amtools fixmate -@ ${samtools_extra_threads} -r -m - - | samtools sort -@ ${samtools_extra_threads} -o ${library}_vs_mt.bam ${library}.bam -
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
	
	"""
	java ${java_options} -jar ${params.bin}picard.jar MarkDuplicates I=${lalnbam} O=${lalnbam.simpleName}.mrkdup.bam M=${lalnbam.simpleName}.mrkdup.txt
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
	tuple path(bam), val(species), val(sample) from sample_bam_ch
	
	output:
	file("*merged.bam") // There needs to be at least one output file, but either/both is ok
	tuple file("${sample}_merged_vs_genome.merged.bam"), val(species_unique) optional true into merged_genome_bam_ch
	tuple file("${sample}_merged_vs_mt.merged.bam"), val(species_unique) optional true into merged_mt_bam_ch
	
	script:
	// First make sure that an input file exists for each type of alignment since could have been removed earlier. If no alignments exist, the whole sample should have been removed previously.
	mtbamlist = ""
	genomebamlist = ""
	species_unique = species[0] // Remove redundant species designators from combining samples
	for (i in bam) {
		category = i.simpleName.split("_vs_")[1]
		if (category == "mt")
			mtbamlist = mtbamlist + " " + i
		else
			genomebamlist = genomebamlist + " " + i
	}
	if (mtbamlist != "" && genomebamlist == "")
		"""
		samtools merge -@ ${params.samtools_extra_threads} ${sample}_merged_vs_mt.merged.bam $mtbamlist
		"""
	else if (genomebamlist != "" && mtbamlist == "")
		"""	
		samtools merge -@ ${params.samtools_extra_threads} ${sample}_merged_vs_genome.merged.bam $genomebamlist
		"""
	else
		"""
		samtools merge -@ ${params.samtools_extra_threads} ${sample}_merged_vs_mt.merged.bam $mtbamlist
		samtools merge -@ ${params.samtools_extra_threads} ${sample}_merged_vs_genome.merged.bam $genomebamlist
		"""
		
}

merged_bam_ch = merged_genome_bam_ch.mix(merged_mt_bam_ch)

process mergedLeftAlignIndels {

	// Left align indels for merged libraries
	
	input:
	tuple path(mrgbam), val(species) from merged_bam_ch
	val java_options from params.java_options
	path mtDNA from params.mtDNA
	path mtDNA_fai from fai_mtDNA_laln_ch
	path genome from params.refseq
	path genome_fai from fai_refseq_laln_ch
	
	output:
	tuple file("${mrgbam.simpleName}.laln.bam"), val(species) into laln_merged_bam_ch
	
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
	java ${java_options} -jar ${params.bin}gatk.jar LeftAlignIndels -I ${mrgbam} -O ${mrgbam.simpleName}.laln.bam -R ${laln_reference}
	"""

}

process mergedMarkDup {

	// Mark duplicates for merged libraries after merging using Picard MarkDuplicates
	
	publishDir "$params.outdir/03_FinalBAMs", mode: 'copy'
	
	input:
	tuple path(laln_mrg_bam), val(species) from laln_merged_bam_ch
	val java_options from params.java_options
	
	output:
	file("${laln_mrg_bam.simpleName}.mrkdup.bam") into mrg_mrkdup_bam_ch
	tuple file("${laln_stem}_vs_genome.mrkdup.bam"), val(species) optional true into final_bam_ch
	tuple file("${laln_stem}_vs_mt.mrkdup.bam"), val(species) optional true into final_mtbam_ch
	
	script:
	laln_stem = laln_mrg_bam.simpleName.split('_vs_')[0]
	"""
	java ${java_options} -jar ${params.bin}picard.jar MarkDuplicates I=${laln_mrg_bam} O=${laln_mrg_bam.simpleName}.mrkdup.bam M=${laln_mrg_bam.simpleName}.mrkdup.txt
	"""

}

process mergedFlagStats {

	// Calculate alignment statistics using SAMtools flagstat
	
	publishDir "$params.outdir/04_FinalFlagStats", mode: 'copy'
	
	input:
	file(mrkdupbam) from mrg_mrkdup_bam_ch
	val samtools_extra_threads from params.samtools_extra_threads
	
	output:
	file "${mrkdupbam.simpleName}.stats.txt"
	
	"""
	${params.bin}samtools flagstat -@ ${samtools_extra_threads} ${mrkdupbam} > ${mrkdupbam.simpleName}.stats.txt
	"""

}

final_species_ch = final_bam_ch.unique().groupTuple(by: 1)
final_mtspecies_ch = final_mtbam_ch.unique().groupTuple(by: 1)



process jointcallVariants {

	// Joint call genomic variants using BCFtools mpileup/call
	
	publishDir "$params.outdir/05_RawVCFs", mode: 'copy'
	
	input:
	tuple path(final_bam), val(species) from final_species_ch
	path genome from params.refseq
	path genome_fai from fai_refseq_laln_ch
	
	output:
	file "${species}_genomic_variants.raw.vcf.gz" into nuVar_ch
	
	"""
	${params.bin}bcftools mpileup -a AD,DP -f $genome -q 20 -Q 20 *.bam | ${params.bin}bcftools call -m -v -Oz -o ${species}_genomic_variants.raw.vcf.gz
	"""

}

process jointcallmtHaplotypes {

	// Joint call mitogenomic haplotypes using BCFtools mpileup/call and vcf2aln
	// Requires minimum allele depth of 3 to be included in alignment
	
	publishDir "$params.outdir/06_MtHaplotypes", mode: 'copy'
	
	input:
	tuple path(final_bam), val(species) from final_mtspecies_ch
	path mtDNA from params.mtDNA
	path mtDNA_fai from fai_mtDNA_laln_ch
	
	output:
	file "${species}_mt_*.fa.gz"
	
	"""
	${params.bin}bcftools mpileup -a AD,DP -f $mtDNA -q 20 -Q 20 *.bam | ${params.bin}bcftools call --ploidy 1 -m -Ov | ${params.bin}vcf2aln.rb --pipe -A 3 -N -o ${species}_mt
	gzip ${species}_mt_*.fa
	"""

}
	
process filternuVar {

	// Filter nuclear variants using VCFtools
	
	publishDir "$params.outdir/07_FiltVCFs", mode: 'copy'
	
	input:
	path raw_vcf from nuVar_ch
	
	output:
	file "${raw_vcf.simpleName}.filt.recode.vcf.gz" into nuVar_filt_ch
	file "${raw_vcf.simpleName}.filt.*"
	
	"""
	${params.bin}vcftools --gzvcf $raw_vcf --out ${raw_vcf.simpleName}.filt --minDP 5 --max-missing 1 --min-alleles 2 --max-alleles 2 --maf 0.01 --remove-indels --recode
	${params.bin}vcftools --vcf ${raw_vcf.simpleName}.filt.recode.vcf --het --out ${raw_vcf.simpleName}.filt
	${params.bin}vcftools --vcf ${raw_vcf.simpleName}.filt.recode.vcf --depth --out ${raw_vcf.simpleName}.filt
	gzip ${raw_vcf.simpleName}.filt.recode.vcf
	"""

}

process mapfilternuVar {
	
	// Filter nuclear variants in regions of low mappability using BEDtools
	
	publishDir "$params.outdir/08_GenMapVCFs", mode: 'copy'
	
	input:
	path filt_vcf from nuVar_filt_ch
	path gm_bed from genmap_ch
	
	output:
	file "${filt_vcf.simpleName}.gm.recode.vcf.gz" into nuVar_gm_ch
	file "${filt_vcf.simpleName}.gm.*"
	
	"""
	${params.bin}bedtools intersect -a $filt_vcf -b $gm_bed -v -header > ${filt_vcf.simpleName}.gm.recode.vcf
	${params.bin}vcftools --vcf ${filt_vcf.simpleName}.gm.recode.vcf --het --out ${filt_vcf.simpleName}.gm
	${params.bin}vcftools --vcf ${filt_vcf.simpleName}.gm.recode.vcf --depth --out ${filt_vcf.simpleName}.gm
	gzip ${filt_vcf.simpleName}.gm.recode.vcf
	"""
	
}
