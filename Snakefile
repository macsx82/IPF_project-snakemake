#Snakefile for the Integrated Prioritization Workflow pipeline
#
#Author: Massimiliano Cocca
#
#
configfile: "config.yaml"

print (config['chr_to_phase'].keys())


def generate_shapeit_out_files(key):
    chr_phased= "%s/%s/chr%s.haps.gz" % (config["output_folder"],config["pop"], key)
    samples= "%s/%s/chr%s.samples" % (config["output_folder"],config["pop"],key)

    return chr_phased,samples

def generate_end_of_pipeline_files(key):
    return "%s/%s/chr%s.pipe.done" % (config["output_folder"],config["pop"],key)


rule all:
    input:
        # lambda wildcards: config["chr_to_phase"][wildcards.chr],
        generate_end_of_pipeline_files("{chr}")
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
        generate_shapeit_out_files("{chr}")
        # chr_phased=config["output_folder"]"/"config["pop"]"/chr{chr}.haps.gz",
        # samples=config["output_folder"]"/"config["pop"]"/chr{chr}.samples"
    threads: 8
    shell:
        "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"

rule pipe_finish:
    input:
        # lambda wildcards: config["chr_to_phase"][wildcards.chr],
        generate_shapeit_out_files("{chr}")
    output:
        generate_end_of_pipeline_files("{chr}")
    shell:
        "touch {output}"
