# Elephant Analysis Pipeline  
Michael G. Campana, 2023  
Smithsonian's National Zoo and Consevation Biology Institute  

This Nextflow [1] pipeline automates the alignment of Illumina sequencing reads against a reference genome and mitogenome using BWA-MEM [2], SAMtools [3-4], and the Genome Analysis Toolkit [5].  

## Citation  
Please cite:  
Prado, N.A., Armstrong, E.E., Brown, J.L., Goldenberg, S.Z., Leimgruber, P., Pearson, V.R., Maldonado, J.E., Campana, M.G. 2023. Genomic resources for Asian (*Elephas maximus*) and African savannah elephant (*Loxodonta africana*) conservation and health research. *Journal of Heredity*. DOI: [10.1093/jhered/esad034](https://doi.org/10.1093/jhered/esad034).  

## License  
This software is licensed under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Installation  
After installing [Nextflow](https://www.nextflow.io/) and [Conda](https://docs.conda.io/en/latest/), download the pipeline using:  
`nextflow pull campanam/Elephants`  

## Configuring the Pipeline  
The `nextflow.config` file included with this repository contains a standard profile for running the pipeline locally. See the Nextflow documentation for assistance in generating a configuration profile for your computing system. The parameters you will need to provide to execute the pipeline are listed in the `params` block. These are:  
`outdir`: Path to the output directory  
`refseq`: Path to the genome reference sequence  
`mtDNA`: Path to the mitochondrial reference sequence  
`samples`: Path to the CSV detailing sequencing libraries (See below)  
`reads`: Path to the folder containing the FASTQ read pairs  
`min_uniq_mapped`: Minimum number of unique mapped reads to retain an alignment file  
`java_options`: String of options for executing Java 1.8  
`email`: Email to send completion status. Set to "NULL" for no email.  

## Sample CSV File  
The pipeline expects a headered CSV file listing samples and libraries with the following columns:  
`Sample,Library,Read1,Read2`  

`Sample` is the name of the sample, while `Library` is in the individual library identification (in case a sample was sequenced more than once). `Read1` and `Read2` give the forward and reverse read file names.  

## Executing the Pipeline  
Execute the pipeline using the following command:  
`nextflow run campanam/Elephants -r main -c <config_file.config>`  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, [1303.3997v2](https://arxiv.org/abs/1303.3997).  
3. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25, 2078-2079. DOI: [10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).  
4. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelves years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
5. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
