
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://zenodo.org/badge/352605728.svg)](https://zenodo.org/badge/latestdoi/352605728)

Snakemake Workflow for Structure Variations Calling and Gene Fusion 
=======================================================================================================================


This is a Snakemake workflow for structure variations calling using delly, tiddit, and sniffles. (More structure variations tools will be included). 

We also use genefuse for calling gene fusions. 

The pipeline uses trimgalore and cutadapt to trim adapters. Align the reads using bwa mem. Based on the SV tool used, we either sort and index or add readgroups and mark duplicates.  
The pipeline annotates, filters, and merges samples. As well as visualize the structure variants into nice plots. 


#### Edit the configfile 

You will need to edit your config file as described below: 

| Config Variable      | Description                      |
| ---------------------| ---------------------------------|
| SAMPLES              | name of file containing your samples names, default: samples.tsv |
| TOOLS                | Leave as it is if you want to run delly, tiddit, sniffles, or remove unwanted tool   |
| BUILD                | For AnnotSV, select genome BUILD required |
| GENOME               | Path to your genome file |
| COHORT               | Name of your Cohort | 
| PAIRED               | True if your samples are paired, false otherwise | 
| RG                   | The READ Group  |
| druggable            | druggable file requried by genefuse | 
| EXCL                 | Excl file required by delly |
| MIN_BP               | Minumum Base Pairs of variant to plot. Default is 20, going lower than 20 will cause plotting issues | 


The pipeline takes samples with a suffix 'r_1.fq.gz' and 'r_2.fq.gz' if the samples are paired. Or it takes samples with suffix 'fq.gz' if the samples are single-end reads. 
Regardless your samples are paired or single-ended, SAMPLES should be listed in samples.tsv without the suffix. 

#### Output 

The pipeline shows some plots to help visualize the SV as follows:

#### Deletion Example 
 
![DEL_chr2_32916570_201284738.png](samplot-out/DEL_chr2_32916570_201284738.png)

##### Duplication Example 

![DUP_chr11_99819753_99820641.png](samplot-out/DUP_chr11_99819753_99820641.png)

##### Inversion Example 


![INV_chr2_79026861_79158701.png](samplot-out/INV_chr2_79026861_79158701.png)


For each sample you will get a a vcf, annotated tsv file in addition to a sub-folder samplot-out in your working directory which contains all the plots of the SVs in your samples and yoursample.html page listing all the SVs per sample. 

The results of gene fusion will be in an html page yoursample_report.html in your working directory.  

#### Run the pipeline

To run the pipeline use:

    snakemake -jn where n is the number of cores. 

For example for using 10 cores, run:

    snakemake -j10


#### Use Conda 

For less frooodiness, to pull automatically the same versions of dependencies use: 

    snakemake -jn --use-conda 

This will pull the same versions of tools we used. Conda has to be installed in your system. 

For example, for 10 cores: 
 
    snakemake -j10 --use-conda 

There is a conda env yaml file for all tools. However, ANNOTSV has to be installed. We couldn't find conda installation for it. 

#### Run rule by rule 

You can run a specific rule using the rule output, for example we can call the delly_vcf rule on a specific sample: 
 
    snakemake -j1 sample.delly.vcf --use-conda  


#### Dry run 

for a dry run use: 

    snakemake -j1 -n 

and you can see the command printed on a dry run using: 

    snakemake -j1 -n -p 

#### Keep going option 


You can try the following to keep going if any issues happen, like no variants is found by one tool: 
    
    snakemake -j1 --keep-going 

#### Cite Us 

If you use the pipeline, please cite us as follows: 

      Sherine Awad. (2022). SherineAwad/StructureVariations: (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.6390009


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

10. Sedlazeck, F. J., Rescheneder, P., Smolka, M., Fang, H., Nattestad, M., Von Haeseler, A., & Schatz, M. C. (2018). Accurate detection of complex structural variations using single-molecule sequencing. Nature methods, 15(6), 461-468. 
