#Snakefile for the Integrated Prioritization Workflow pipeline
#
#Author: Massimiliano Cocca
#
#
configfile: "config.yaml"
#First we need to phase our data
#preferred input files format are vcf, but we will handle also plink formatted files
rule phase:
    input:
        # g_map=config["genetic_map_chr"],
        lambda wildcards: config["chr_to_phase"][wildcards.chr]
        # chr_to_phase="{input_folder}/{chr}.vcf.gz"
    params:
        input_f=config["input_folder"],
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr{chr}_combined_b37.txt"
    output:
        chr_phased=config["output_folder"]"/"config["pop"]"/chr{chr}.haps.gz"
        samples=config["output_folder"]"/"config["pop"]"/chr{chr}.samples"
    threads: 8
    shell:
        "shapeit -V {input_f}/{input} -M {g_map} -O {chr_phased} {samples} -T {threads}"
