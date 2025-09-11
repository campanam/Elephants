# Elephant Analysis Pipeline  
<img align="right" src="NZP-20180628-596SB_thumb.jpg">  

Michael G. Campana, 2023-2025  
Smithsonian's National Zoo & Conservation Biology Institute  

This Nextflow [1] pipeline automates the alignment of Illumina sequencing reads against a reference genome and mitogenome using BWA-MEM [2], SAMtools [3-4], and the Genome Analysis Toolkit [5]. It optionally performs PSMC analysis [6], read trimming using AdapterRemoval v2 [7] and circular mitogenome alignment using CircularMapper [8].  

## Citation  
Please cite:  
Prado, N.A., Armstrong, E.E., Brown, J.L., Goldenberg, S.Z., Leimgruber, P., Pearson, V.R., Maldonado, J.E., Campana, M.G. 2023. Genomic resources for Asian (*Elephas maximus*) and African savannah elephant (*Loxodonta africana*) conservation and health research. *Journal of Heredity*. 114(5): 529–538. DOI: [10.1093/jhered/esad034](https://doi.org/10.1093/jhered/esad034).  

## License  
This software is licensed under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Installation  
After installing [Nextflow](https://www.nextflow.io/) and [Conda](https://docs.conda.io/en/latest/), download the pipeline using:  
`nextflow pull campanam/Elephants`  

## Configuring the Pipeline  
The `nextflow.config` file included with this repository contains a standard profile for running the pipeline locally. See the Nextflow documentation for assistance in generating a configuration profile for your computing system. The parameters you will need to provide to execute the pipeline are listed in the `params` block. These are:  
`outdir`: Path to the output directory  
`refseq`: Path to the genome reference sequence  
`csi`: Use csi index (e.g. for chr > 512 Mb) (true or false)  
`circular_mtDNA`: Perform secondary independent circular alignment of mitochondrial DNA  
`mtDNA`: Path to the mitochondrial reference sequence  
`mtDNA_ID`: Name of mtDNA sequence in mitochondrial reference sequence file  
`samples`: Path to the CSV detailing sequencing libraries (See below)  
`reads`: Path to the folder containing the FASTQ read pairs  
`read_trimming`: Trim reads using AdapterRemoval v2 (true or false)  
`trimparams`: Parameters for read-trimming. This pipeline expects paired-read input, so do not collapse reads.  
`min_uniq_mapped`: Minimum number of unique mapped reads to retain an alignment file  
`java_options`: String of options for Java executables.  
`markDuplicates`: Choice of "picard", "samtools" or "sambamba" for duplicate marking  
`gatk`: Perform genotyping using GATK (true or false)  
`email`: Email to send completion status. Set to "NULL" for no email.  
`psmc`: Run PSMC analysis (true or false)  
`psmc_opts`: Parameter line for PSMC program  
`psmc_mpileup_opts`: BCFtools mpileup filter options for PSMC  
`psmc_vcfutils_opts`: vcfutils.pl filter options for PSMC  
`psmc_psmcfa_opts`: fq2psmcfa options for PSMC  
`psmc_bootstrap`: Number of PSMC bootstrap replicates  
`psmc_plot_opts`: Options for psmc_plot.pl  

## Sample CSV File  
The pipeline expects a headered CSV file listing samples and libraries with the following columns:  
`Sample,Library,Read1,Read2,Adapter1,Adapter2`  

`Sample` is the name of the sample, while `Library` is in the individual library identification (in case a sample was sequenced more than once). `Read1` and `Read2` give the forward and reverse read file names. `Adapter1` and `Adapter2` give the expected forward and reverse adapter sequences.  

*NB: The adapter columns can be omitted if no read-trimming will be performed.*  

## Executing the Pipeline  
Execute the pipeline using the following command:  
`nextflow run campanam/Elephants -r main -c <config_file.config> -profile standard`  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, [1303.3997v2](https://arxiv.org/abs/1303.3997).  
3. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, __25__, 2078-2079. DOI: [10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).  
4. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelves years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
5. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
6. Li, H., Durbin, R. (2011) Inference of human population history from individual whole-genome sequences. *Nature*, __475__, 493-496. DOI [10.1038/nature10231](https://www.nature.com/articles/nature10231).  
7. Schubert, M., Lindgreen, S., Orlando, L. (2016) AdapterRemoval v2: rapid adapter trimming, identification, and read merging. *BMC Research Notes*, __9__, 88. DOI: [10.1186/s13104-016-1900-2](https://doi.org/10.1186/s13104-016-1900-2).  
8. Peltzer, A., Jäger, G., Herbig, A., Seitz, A., Kniep, C., Krause, J., Nieselt, K. (2016) EAGER: efficient ancient genome reconstruction. *Genome Biology*, __17__, 60. DOI: [10.1186/s13059-016-0918-z](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0918-z).  

Image Credit: Skip Brown. 2018. Smithsonian's National Zoo & Conservation Biology Institute. Smithsonian Institution. https://www.si.edu/object/asian-elephant:nzp_NZP-20180628-596SB.  
