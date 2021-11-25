
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![Documentation Status](https://readthedocs.org/projects/ribofilio/badge/)](http://dammit.readthedocs.io/en/latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Snakemake Workflow for Structure Variations Calling and Gene Fusion 
=======================================================================================================================


This is a Snakemake workflow for structure variations calling using delly, tiddit, and lumpy calling for structure variations SV (More structure variations tools will be included). 

We also use genefuse for calling gene fusions. 

The pipeline uses trimgalore and cutadapt to trim adapters. Align the reads using Bowtie2. 
Then we use from GATK4 pipeline add read groups,and mark duplicates.  

We use ANNOTSV for annotating the structure variants and SURVIVOR for filtering and merging as needed. 

#### Run the pipeline 

To run the pipeline, you need to change appropriately the config file, and use:

    snakemake -jn where n is the number of cores. 

For example for using 10 cores, run:

    snakemake -j10


#### Use Conda 

For less frooodiness, to pull automatically the same versions of dependencies use: 

    snakemake -jn --use-conda 

This will pull the same versions of tools we used. Conda has to be installed in your system. 

For example, for 10 cores: 
 
    snakemake -j10 --use-conda 

#### Run rule by rule 

You can run a specific rule using the rule output, for example we can call the delly_vcf rule on a specific sample: 
 
    snakemake -j1 sample.delly.vcf --use-conda  

where sample.r_1.fq.gz and sample.s_2.r_2.fq.gz should exist in the directory. 


#### Dry run 

for a dry run use: 

    snakemake -j1 -n 

and you can see the command printed on a dry run using: 

    snakemake -j1 -n -p 

#### Keep going option 


You can try the following to keep going if any issues happen, like no variants is found by one tool: 
    
    snakemake -j1 --keep-going 


## TODO 

1. Add more parameters in the config file to replace many default parameters 

2. Add more tools to call structure variations that suite WES and WGS 

3. Add annotation rules for structure variations 

4. Add  more gene fusion tools 

5. Add more comments on the snakefile

## References

1. Brouard, Jean-Simon, Flavio Schenkel, Andrew Marete, and Nathalie Bissonnette. "The GATK joint genotyping workflow is appropriate for calling variants in RNA-seq experiments." Journal of animal science and biotechnology 10, no. 1 (2019): 1-6.

2. Van der Auwera, Geraldine A., Mauricio O. Carneiro, Christopher Hartl, Ryan Poplin, Guillermo Del Angel, Ami Levy‐Moonshine, Tadeusz Jordan et al. "From FastQ data to high‐confidence variant calls: the genome analysis toolkit best practices pipeline." Current protocols in bioinformatics 43, no. 1 (2013): 11-10.

3. Poplin, R., Ruano-Rubio, V., DePristo, M. A., Fennell, T. J., Carneiro, M. O., Van der Auwera, G. A., ... & Banks, E. (2018). Scaling accurate genetic variant discovery to tens of thousands of samples. BioRxiv, 201178.

4. Rausch, T., Zichner, T., Schlattl, A., Stütz, A. M., Benes, V., & Korbel, J. O. (2012). DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics, 28(18), i333-i339.

5. Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27(21), 2987-2993.

6. Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6.

7. Chen, S., Liu, M., Huang, T., Liao, W., Xu, M., & Gu, J. (2018). GeneFuse: detection and visualization of target gene fusions from DNA sequencing data. International journal of biological sciences, 14(8), 843–848. https://doi.org/10.7150/ijbs.24626

8. Geoffroy, V., Herenger, Y., Kress, A., Stoetzel, C., Piton, A., Dollfus, H., & Muller, J. (2018). AnnotSV: an integrated tool for structural variations annotation. Bioinformatics, 34(20), 3572-3574.

9. Jeffares, D. C., Jolly, C., Hoti, M., Speed, D., Shaw, L., Rallis, C., ... & Sedlazeck, F. J. (2017). Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast. Nature communications, 8(1), 1-11.

10. Layer, R. M., Chiang, C., Quinlan, A. R., & Hall, I. M. (2014). LUMPY: a probabilistic framework for structural variant discovery. Genome biology, 15(6), 1-19.

