#Snakefile for the Integrated Prioritization Workflow pipeline
#
#Author: Massimiliano Cocca
#
#
configfile: "config.yaml"

print (config['chr_to_phase'].keys())


def generate_shapeit_out_files(key):
    chr_phased= "%s/%s/%s/chr%s.haps.gz" % (config["output_folder"],config["pop"],key,key)
    samples= "%s/%s/%s/chr%s.samples" % (config["output_folder"],config["pop"],key,key)

    return chr_phased,samples

def generate_end_of_pipeline_files(key):
    return "%s/%s/chr%s.pipe.done" % (config["output_folder"],config["pop"],key)


rule all:
    input:
        # config["output_folder"]+"/"+config["pop"]+"/{chr}.pipe.done"
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"
        # generate_end_of_pipeline_files("{chr}")
        # lambda wildcards: config["chr_to_phase"][wildcards.chr]
#First we need to phase our data
#preferred input files format are vcf, but we will handle also plink formatted files
rule phase:
    input:
        # g_map=config["genetic_map_chr"],
        # lambda wildcards: config["chr_to_phase"][wildcards.chr]
        # chr_to_phase="{input_folder}/{chr}.vcf.gz"
        # "{params.input_f}/{chr}.vcf.gz"
        # config["input_folder"] + "/{chr}.vcf.gz"
        config["input_folder"] + "/" + config["chr"]+ ".vcf.gz"
    params:
        input_f=config["input_folder"],
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        # base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]
    output:
        # generate_shapeit_out_files("{input.chr}")
        generate_shapeit_out_files(config["chr"])
        # chr_phased=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+".haps.gz",
        # samples= params[base_out] + "/chr"+config["chr"]+".samples"
    # threads: 2
    shell:
        # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        # "{config[shapeit_path]} -V {input} -M {params.g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        "{config[shapeit_path]} -V {input} -M {params.g_map} -O {output[0]} {output[1]} -T {threads}"
        # "touch {output.chr_phased} {output.samples}"

rule relate_poplabels:
    input:
        config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+".samples"
    params:
        input_f=config["input_folder"],
        base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]
    output:
        # generate_shapeit_out_files("{input.chr}")
        # generate_shapeit_out_files("{chr}")
        config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+".poplabels"
    shell:
        # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        "(echo \"sample population group sex\";tail -n+3 {input} | awk '{{OFS=" "}}{{print $1,{config[pop_group]},{config[pop]},$6}}') > {params.base_out}/{output}"

rule relate:
    input:
        generate_shapeit_out_files(config["chr"])
        # chr_phased=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".haps.gz",
        # samples=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".samples"
    params:
        input_f=config["input_folder"],
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"],
        # out_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate"
        out_prefix="chr"+config["chr"]+"_relate"
    output:
        # generate_shapeit_out_files("{input.chr}")
        # generate_shapeit_out_files("{chr}")
        # base_out + "/chr"+config["chr"]+"_relate"
        # multiext(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate", ".anc", ".mut")
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate{ext}", ext=[".anc", ".mut"])
    shell:
        # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        # "{config[relate_path]}/bin/Relate --mode All --m 1.25e-8 -N 30000 --haps {input.chr_phased} --sample {input.samples} --map {params.g_map} --seed {config[relate_seed]} -o {params.out_prefix}"
        "cd {params.base_out};{config[relate_path]}/bin/Relate --mode All --m 1.25e-8 -N 30000 --haps {input[0]} --sample {input[1]} --map {params.g_map} --seed {config[relate_seed]} -o {params.out_prefix}"

rule relate_pop_s_est:
    input:
        # relate_files=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+"_relate",
        relate_files=expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate{ext}", ext=[".anc", ".mut"]),
        poplabel_file=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+".poplabels"
    params:
        input_f=config["input_folder"],
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        # base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]
        in_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate",
        out_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate_popsize"
    output:
        # generate_shapeit_out_files("{input.chr}")
        # generate_shapeit_out_files("{chr}")
        # multiext(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize", ".pdf",".anc.gz",".mut.gz",".dist",".coal",".bin","_avg.rate")
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize{ext}", ext=[".pdf",".anc.gz",".mut.gz",".dist",".coal",".bin","_avg.rate"])
        # base_out + "/chr"+config["chr"]+"_relate_popsize"
    shell:
        # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        "{config[relate_path]}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i {params.in_prefix} --poplabels{input.poplabel_file} --m 1.25e-8 \ "
        " --seed {config[relate_seed]} -o {params.out_prefix}"

# rule relate_mut_rate_est:
#     input:
#         chr_phased=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".haps.gz",
#         samples=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".samples"
#     params:
#         input_f=config["input_folder"],
#         g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt"
#     output:
#         # generate_shapeit_out_files("{input.chr}")
#         # generate_shapeit_out_files("{chr}")
#         config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+"_relate",
#     shell:
#         # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
#         "{config[relate_path]} --mode All --m 1.25e-8 -N 30000 --haps {input.chr_phased} --sample {input.samples} \ "
#         "--map {params.g_map} --seed 1 -o {output}"

rule pipe_finish:
    # params:
    #     base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]
    input:
        # lambda wildcards: config["chr_to_phase"][wildcards.chr],
        # generate_shapeit_out_files("{chr}")
        # chr_phased=config["output_folder"] + "/" + config["pop"] + "/chr{chr}.haps.gz",
        # samples=config["output_folder"] + "/" + config["pop"] + "/chr{chr}.samples"
        # chr_phased=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".haps.gz",
        # samples=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".samples"
        # config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+"_relate"
        # base_out + "/chr"+config["chr"]+"_relate_popsize"
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize{ext}", ext=[".pdf",".anc.gz",".mut.gz",".dist",".coal",".bin","_avg.rate"])
    output:
        # generate_end_of_pipeline_files("{chr}")
        # config["output_folder"]+"/"+config["pop"]+"/{chr}.pipe.done"
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"
    shell:
        "touch {output}"
