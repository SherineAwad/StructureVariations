configfile:"config.yaml" 

with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)


for f in config['TOOL']: 
    file = config['COHORT']+"."+f +".list"
    fp = open(file, "w+") 
    for i in SAMPLES: 
        output = i+"."+ f +".vcf" 
        print(output, file =fp)
    fp.close()

rule all: 
    input:
       expand("{sample}.sam", sample = SAMPLES),
       expand("{sample}.sorted.bam", sample = SAMPLES),
       expand("{sample}_report.html",  sample =SAMPLES), 
       expand("{sample}.{sv}.vcf", sample = SAMPLES, sv = config['TOOL'] ),
       expand("{sample}.{sv}.annotated.vcf.tsv", sample = SAMPLES , sv = config['TOOL']),
       expand("samplot-out/{sample}.{sv}.html", sample = SAMPLES , sv = config['TOOL']),
       expand("{COHORT}.{SV}.vcf", COHORT = config['COHORT'], SV= config['TOOL'])

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
           "galore/{sample}_trimmed.fq.gz",
           genome = config['GENOME']
        output:
            '{sample}.sam'
        conda: 'env/env/env-align.yaml'
        shell:
           "bwa mem {input.genome} {input[0]} > {output}" 

rule sam_bam: 
    input: 
        "{sample}.sam"
    output: 
        "{sample}.bam"
    conda: 'env/env-tools.yaml'
    shell: 
         """ 
         samtools view -S -b {input} > {output}
         """

rule sort_index: 
     input: 
       "{sample}.bam" 
     output: 
       "{sample}.sorted.bam"
     conda: 'env/env-tools.yaml'
     shell: 
         """
         samtools sort {input} -o {output}
         samtools index {output}
         """  
rule AddRG:
    input:
       '{sample}.sam',
       '{sample}_r1.fastq.gz'
    output:
       '{sample}.RG.sam'
    log: "logs/{sample}.addRG.log"
    benchmark: "logs/{sample}.AddRG.benchmark"
    conda: 'env/env-picard.yaml'
    params:
        PL = config['PL']
    shell:
        """
         python src/RG.py -s {input[0]} -f {input[1]} -o {output} -p = {params.PL}
        """

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
         EXCL = config['EXCL']
      output: 
         config['EXCL']
      shell: 
         "wget https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/{params}"
 
rule delly_bcf:
     input:
        "{SAMPLE}.recalibrated.bam",
        genome = config['GENOME'],
        EXCL = config['EXCL']
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
     conda: 'env/env-tools.yaml'
     shell:
         "bcftools view {input} > {output}"


rule tiddit_vcf: 
    input: 
        "{SAMPLE}.recalibrated.bam"
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
       "{SAMPLE}.recalibrated.bam"
    output:
       "{SAMPLE}.sniffles.vcf"
    conda: "env/env-sniffles.yaml"
    shell:
        """
        sniffles -i {input} --vcf {output}
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
      expand("{sample}.{sv}.vcf" ,sample = SAMPLES, sv =config['TOOL'])
   params:  
      "{COHORT}.{sv}.list" 
   output: 
      "{COHORT}.{sv}.vcf"
   conda: 'env/env-survivor.yaml'
   shell: 
      """
        SURVIVOR merge {params[0]}  1000 10 1 1 0 30 {output[0]}
      """ 

rule plot: 
    input: 
       "{sample}.{sv}.vcf",
       "{sample}.recalibrated.bam" 
    params: 
        config['MIN_BP']
    output:
        "samplot-out/{sample}.{sv}.html"
    conda: 'env/env-plot.yaml'
    shell: 
       """
       samplot vcf --vcf {input[0]} --min_bp {params} -O png -b {input[1]} && mv samplot-out/index.html {output}   
       """ 

