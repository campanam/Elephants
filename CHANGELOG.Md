# Elephant Analysis Pipeline CHANGELOG  

Michael G. Campana, 2023-2025  
Smithsonian's National Zoo and Conservation Biology Institute  

### Version 0.3.0  
Conversion to DSL2  
Option to perform circular mitochondrial alignment  
Options to trim reads  
Defaults to using Mamba to install dependencies  
Standardization with other pipelines  
Bumped SAMtools/BCFtools versions  
Bumped GATK version  
Added Picard/sambamba mark duplicates options  
Bumped BWA version  
Bumpled zlib version  
Added csi index option  
Added MapQ filtering and statistics  
Statistics outputs depth and coverage  

### Version 0.2.1  
Final bams that do not need merging are published under the final names for simplicity  

### Version 0.2.0  
Added PSMC support  
GATK genotype calling can be turned off  

## Version 0.1.0a  
Updated configuration refidx paths for static use on SI/HPC  
Added SAMtools to conda environments for GATK genotyping  

## Version 0.1.0  
First released version  
