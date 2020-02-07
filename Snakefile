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

def generate_relate_clean_in_files(key):
    chr_clean_phased= "%s/%s/%s/chr%s.relate_clean.haps.gz" % (config["output_folder"],config["pop"],key,key)
    clean_samples= "%s/%s/%s/chr%s.relate_clean.samples" % (config["output_folder"],config["pop"],key,key)

    return chr_clean_phased,clean_samples

def generate_pop_size_threshold_est(sample_file):
    #we need to use the haplotype file, to count how many haps we have.
    #we can use the sample file, which is lighter
    return (sum(1 for line in open('%s' %(sample_file),'r')) - 1 )*2


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

#section to prepare relate input files
rule relate_prepare_input:
    input:
        generate_shapeit_out_files(config["chr"]),
        poplabel_file=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+".poplabels"
        # chr_phased=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".haps.gz",
        # samples=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+".samples"
    params:
        out_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+".relate_clean",
        ancestor_chr=config["ancestor_ref_path"]+"/human_ancestor_"+config["chr"]+".fa"
    output:
        generate_relate_clean_in_files(config["chr"])
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+".relate_clean.{ext}", ext=["haps.gz", "samples"])
    shell:
        """
        set +e
        {config[relate_path]}/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps {input[0]} --sample {input[0]} \
         --ancestor {params.ancestor_chr} \
         --poplabels {input.poplabel_file} -o {params.out_prefix}
        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            echo "No error found..exiting correctly"
            exit 0
        else
            echo "WARNING....The software raised some errors or warning, be careful and check the results."
            exit 0
        fi
        """


rule relate_poplabels:
    input:
        # config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+".samples"
        generate_shapeit_out_files(config["chr"])
    params:
        input_f=config["input_folder"],
        # base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]
    output:
        # generate_shapeit_out_files("{input.chr}")
        # generate_shapeit_out_files("{chr}")
        config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+".poplabels"
    priority: 1
    shell:
        # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        "(echo \"sample population group sex\";tail -n+3 {input[1]} | awk '{{OFS=\" \"}}{{print $1,\"{config[pop]}\",\"{config[pop_group]}\",0}}') > {output}"

rule relate:
    input:
        # generate_shapeit_out_files(config["chr"])
        generate_relate_clean_in_files(config["chr"])
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
        """
        set +e
        cd {params.base_out};
        {config[relate_path]}/bin/Relate --mode All -m 1.25e-8 -N 30000 --haps {input[0]} --sample {input[1]} --map {params.g_map} --seed {config[relate_seed]} -o {params.out_prefix}
        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            echo "No error found..exiting correctly"
            exit 0
        else
            echo "WARNING....The software raised some errors or warning, be careful and check the results."
            exit 0
        fi
        """
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
        out_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate_popsize",
        # relate_threshold=generate_pop_size_threshold_est("{input.poplabel_file}")
        relate_threshold=0
    output:
        # generate_shapeit_out_files("{input.chr}")
        # generate_shapeit_out_files("{chr}")
        # multiext(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize", ".pdf",".anc.gz",".mut.gz",".dist",".coal",".bin","_avg.rate")
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize{ext}", ext=[".pdf",".anc.gz",".mut.gz",".dist",".coal",".bin","_avg.rate"])
        # base_out + "/chr"+config["chr"]+"_relate_popsize"
    shell:
        # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        """
        set +e
        {config[relate_path]}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i {params.in_prefix} --poplabels {input.poplabel_file} -m 1.25e-8 \
         --threshold {params.relate_threshold} --seed {config[relate_seed]} -o {params.out_prefix}
        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            echo "No error found..exiting correctly"
            exit 0
        else
            echo "WARNING....The software raised some errors or warning, be careful and check the results."
            exit 0
        fi
        """

rule relate_mut_rate_est:
    input:
        # relate_files=config["output_folder"] + "/" + config["pop"] + "/chr"+config["chr"]+"_relate",
        relate_files=expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate{ext}", ext=[".anc", ".mut"]),
    params:
        input_f=config["input_folder"],
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        # base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]
        in_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate",
        out_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate_mut_rate"
    output:
        # generate_shapeit_out_files("{input.chr}")
        # generate_shapeit_out_files("{chr}")
        # multiext(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize", ".pdf",".anc.gz",".mut.gz",".dist",".coal",".bin","_avg.rate")
        config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_mut_rate_avg.rate"
        # base_out + "/chr"+config["chr"]+"_relate_popsize"
    shell:
        # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        """
        set +e
        {config[relate_path]}/bin/RelateMutationRate --mode Avg -i {params.in_prefix} -o {params.out_prefix}
        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            echo "No error found..exiting correctly"
            exit 0
        else
            echo "WARNING....The software raised some errors or warning, be careful and check the results."
            exit 0
        fi        
        """

rule relate_detect_selection:
    input:
        relate_files=expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize{ext}", ext=[".anc.gz", ".mut.gz"])
    params:
        in_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate_popsize",
        out_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+"_relate_pos_sel"
    output:
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_pos_sel{ext}", ext=[".freq",".lin",".sele"])
    shell:
        """
        set +e
        {config[relate_path]}/scripts/DetectSelection/DetectSelection.sh \
                 -i {params.in_prefix} \
                 -o {params.out_prefix} \
                 -m 1.25e-8 \
                 --years_per_gen 28
        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            echo "No error found..exiting correctly"
            exit 0
        else
            echo "WARNING....The software raised some errors or warning, be careful and check the results."
            exit 0
        fi        
        """

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
        config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_mut_rate_avg.rate",
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_popsize{ext}", ext=[".pdf",".anc.gz",".mut.gz",".dist",".coal",".bin","_avg.rate"]),
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_pos_sel{ext}", ext=[".freq",".lin",".sele"])
    output:
        # generate_end_of_pipeline_files("{chr}")
        # config["output_folder"]+"/"+config["pop"]+"/{chr}.pipe.done"
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"
    shell:
        "touch {output}"
