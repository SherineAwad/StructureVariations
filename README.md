Snakemake Workflow for Structure Variations Calling using Delly and Tiddit 
=======================================================================================================================


This is a Snakemake workflow for structure variations calling using delly and tiddit calling for structure variations SV (More structure variations tools will be included). 

The pipeline uses trimgalore and cutadapt to trim adapters. Align the reads using Bowtie2. 
Then GATK4 pipeline follows: add read groups, mark duplicates. We are testing to see whether base recalibration in SV is better. 

The pipleine also uses the dedupped reads to call structure variations using delly. 

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

4. Add gene fusion call 

5. Add more comments on the snakefile

## References

1. Brouard, Jean-Simon, Flavio Schenkel, Andrew Marete, and Nathalie Bissonnette. "The GATK joint genotyping workflow is appropriate for calling variants in RNA-seq experiments." Journal of animal science and biotechnology 10, no. 1 (2019): 1-6.

2. Van der Auwera, Geraldine A., Mauricio O. Carneiro, Christopher Hartl, Ryan Poplin, Guillermo Del Angel, Ami Levy‐Moonshine, Tadeusz Jordan et al. "From FastQ data to high‐confidence variant calls: the genome analysis toolkit best practices pipeline." Current protocols in bioinformatics 43, no. 1 (2013): 11-10.

3. Poplin, R., Ruano-Rubio, V., DePristo, M. A., Fennell, T. J., Carneiro, M. O., Van der Auwera, G. A., ... & Banks, E. (2018). Scaling accurate genetic variant discovery to tens of thousands of samples. BioRxiv, 201178.

4. Rausch, T., Zichner, T., Schlattl, A., Stütz, A. M., Benes, V., & Korbel, J. O. (2012). DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics, 28(18), i333-i339.

5. Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27(21), 2987-2993.

6. Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6.

