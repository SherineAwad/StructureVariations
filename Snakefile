configfile: "config.yaml"

ruleorder: trim > tosam > AddRG > dedup > recalibrate > get_excl > delly_bcf > delly_vcf > tiddit_vcf

rule all: 
    input:
       expand("{sample}.tiddit.vcf", sample = config['SAMPLES'] ),       
       expand("{sample}.delly.vcf", sample = config['SAMPLES'] )


rule trim: 
    input: 
       r1 = "{sample}.r_1.fq.gz",
       r2 = "{sample}.r_2.fq.gz"
    output: 
      val1 = "galore/{sample}.r_1_val_1.fq.gz",
      val2 = "galore/{sample}.r_2_val_2.fq.gz"
    shell: 
        """
         trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
        """ 

rule tosam:
    input:
        genome = config['GENOME'],
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    output:
        '{sample}.sam'
    shell:
        "bowtie2 -x {input.genome} -1 {input.r1} -2 {input.r2} -S {output}"

rule AddRG: 
    input: 
       '{sample}.sam'
    output: 
       '{sample}.RG.sam' 
    params: 
        RG = config['RG']
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=@{params} RGSM={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={params}_{wildcards.sample} VALIDATION_STRINGENCY=SILENT" 


rule dedup: 
     input: 
         '{sample}.RG.sam'
     output:
       '{sample}.dedupped.bam',
       '{sample}.output.metrics'
     shell:
        "picard MarkDuplicates I={input} O={output[0]} CREATE_INDEX=true M={output[1]}"


rule recalibrate: 
    input:  
         sample = '{sample}.dedupped.bam', 
         genome = config['GENOME']
    output:
         '{sample}.recal_data.table',
         '{sample}.recalibrated.bam'
    params: 
        mem = "-Xmx800g",
        knownsites1 = config['DBSNP'],
        knownsites2 = config['INDELS'],
        knownsites3 = config['GOLD_STANDARD']
   
    shell:
       """
       gatk --java-options {params.mem} BaseRecalibrator -I {input.sample} -R {input.genome} --known-sites {params.knownsites1} --known-sites {params.knownsites2} --known-sites {params.knownsites3} -O {output[0]}
       gatk --java-options {params.mem} ApplyBQSR  -I {input.sample} -R {input.genome} --bqsr-recal-file {output[0]} -O {output[1]} 
       """ 


rule get_excl: 
      params: 
         EXCL = config['EXCL_HG38']
      output: 
         config['EXCL_HG38']
      shell: 
         """
           wget https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/{params}
         """
      
rule delly_bcf:
     input:
        "{SAMPLE}.dedupped.bam",
        genome = config['GENOME'],
        EXCL = config['EXCL_HG38']
     output:
        "{SAMPLE}.delly.bcf"
     shell:
       """
       delly call -x {input.EXCL}  -o {output} -g {input.genome} {input[0]}
       """


rule delly_vcf:
     input:
        "{SAMPLE}.delly.bcf"
     output:
        "{SAMPLE}.delly.vcf"
     shell:
         "bcftools view {input} > {output}"


rule tiddit_vcf: 
    input: 
        "{SAMPLE}.dedupped.bam"
    params: 
        "{SAMPLE}.tiddit"
    output: 
         "{SAMPLE}.tiddit.vcf"
    shell: 
      "tiddit --sv --bam {input} -o {params}" 
