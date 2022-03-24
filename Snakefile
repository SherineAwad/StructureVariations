configfile:"config.yaml" 

for f in config['TOOL']: 
    file = config['COHORT']+"."+f +".list"
    fp = open(file, "a+") 
    for i in config['SAMPLES']: 
        output = i+"."+ f +".vcf" 
        print(output, file =fp)
    fp.close()

rule all: 
    input:
       expand("{sample}.sam", sample = config['SAMPLES']),
       expand("{sample}.{sv}.vcf", sample = config['SAMPLES'], sv = config['TOOL'] ),
       expand("{COHORT}.{SV}.vcf", COHORT=config['COHORT'], SV = config['TOOL']),
       expand("{COHORT}.{SV}.list", COHORT=config['COHORT'], SV = config['TOOL']),
       expand("{sample}.{sv}.annotated.vcf.tsv", sample = config['SAMPLES'] , sv = config['TOOL']),
       expand("{sample}.{sv}_output/{sample}.{sv}.html", sample = config['SAMPLES'] , sv = config['TOOL'])
 

if config['PAIRED']:
    rule trim:
       input:
           r1 = "{sample}.r_1.fq.gz",
           r2 = "{sample}.r_2.fq.gz"
       output:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """
    rule tosam:
        input:
            genome = config['GENOME'],
            r1 = "galore/{sample}.r_1_val_1.fq.gz",
            r2 = "galore/{sample}.r_2_val_2.fq.gz"
        output:
            '{sample}.sam'
        conda: 'env/env-align.yaml'
        shell:
            "bwa mem {input.genome} {input.r1} {input.r2} > {output}"
else:
     rule trim:
       input:
           "{sample}.trimmed.fq.gz",
       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input}
           """
     rule tosam:
        input:
           "galore/{sample}_trimmed.fq.gz"
        output:
            '{sample}.sam'
        conda: 'env/env-align.yaml'
        shell:
           "bwa mem {input.genome} {input} > {output}" 

rule sam_bam: 
    input: 
        "{sample}.sam"
    output: 
        "{sample}.bam"
    shell: 
         """ 
         samtools view -S -b {input} > {output}
         samtools index {input}
         """
rule sort_index: 
     input: 
       "{sample}.bam" 
     output: 
       "{sample}.sorted.bam"
     shell: 
         """
         samtools sort {input} -o {output}
         samtools index {output}
         """  
rule AddRG: 
    input: 
       '{sample}.sam'
    output: 
       '{sample}.RG.sam' 
    params: 
        RG = config['RG']
    conda: 'env/env-picard.yaml'
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=@{params} RGSM={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={params}_{wildcards.sample} VALIDATION_STRINGENCY=SILENT" 


rule dedup: 
     input: 
         '{sample}.RG.sam'
     output:
       '{sample}.dedupped.bam',
       '{sample}.output.metrics'
     conda: 'env/env-picard.yaml'
     shell:
        "picard MarkDuplicates I={input} O={output[0]} CREATE_INDEX=true M={output[1]}"


rule get_excl: 
      params: 
         EXCL = config['EXCL_HG38']
      output: 
         config['EXCL_HG38']
      shell: 
         "wget https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/{params}"
 
rule delly_bcf:
     input:
        "{SAMPLE}.dedupped.bam",
        genome = config['GENOME'],
        EXCL = config['EXCL_HG38']
     output:
        "{SAMPLE}.delly.bcf"
     conda: 'env/env-delly.yaml'
     shell:
       "delly call -x {input.EXCL}  -o {output} -g {input.genome} {input[0]}"

rule delly_vcf:
     input:
        "{SAMPLE}.delly.bcf"
     output:
        "{SAMPLE}.delly.vcf"
     conda: 'env/env-delly.yaml'
     shell:
         "bcftools view {input} > {output}"


rule tiddit_vcf: 
    input: 
        "{SAMPLE}.dedupped.bam"
    params: 
        "{SAMPLE}.tiddit"
    output: 
         "{SAMPLE}.tiddit.vcf",
    conda: 'env/env-tiddit.yaml'
    shell: 
      """ 
         tiddit --sv --bam {input} -o {params}
      """

rule sniffles_vcf:
    input:
       "{SAMPLE}.sorted.bam"
    output:
       "{SAMPLE}.sniffles.vcf"
    shell:
        """
        sniffles -m {input} -v  {output}
        """

rule gene_fuse: 
    input:
       val1 = "galore/{sample}.r_1_val_1.fq.gz",
       val2 = "galore/{sample}.r_2_val_2.fq.gz",
       genome = config['GENOME']
    params:   
       druggable = config['druggable'] 
    output: 
        html = "{sample}_report.html",
        result = "{sample}_result"
    conda: 'env/env-fuse.yaml' 
    shell:
        """
        wget https://raw.githubusercontent.com/OpenGene/GeneFuse/master/genes/{params.druggable}
        genefuse -r {input.genome}  -f {params.druggable} -1 {input.val1} -2 {input.val2} -h {output.html} > {output.result}
        """

rule annotate: 
    input: 
          "{sample}.{sv}.vcf" 
    output: 
          "{sample}.{sv}.annotated.vcf.tsv"
    params:
        build = config['BUILD']
    shell:
        "$ANNOTSV/bin/AnnotSV -SVinputFile ./{input} -outputFile ./{output} -genomeBuild {params.build}"         

rule SURVIVOR_LIST:
   input:  
      "{COHORT}.{SV}.list",
   output: 
      "{COHORT}.{SV}.vcf"
   shell: 
      """
        SURVIVOR merge {input[0]}  1000 10 1 1 0 30 {output[0]}
      """ 


rule plot: 
    input: 
       "{sample}.{sv}.vcf", 
       "{sample}.dedupped.bam" 
    params: 
       "{sample}.{sv}_output"
    output: 
        "{sample}.{sv}_output/{sample}.{sv}.html"
    shell: 
       """
       samplot vcf --vcf {input[0]} -O png -b {input[1]} 
       mv samplot-out {params}  
       mv {params}/index.html {output}
       """ 
