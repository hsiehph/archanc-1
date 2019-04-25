import os

# shell.prefix("source config.sh; set -eo pipefail ; ")

configfile: "config.json"

if not os.path.exists("log"):
	os.makedirs("log")

modelNAME = config["modelNAME"]
MNM_dist = config["MNM_dist"]
MNM_frac = config["MNM_frac"]
replicates = config["replicates"]
GeneticMap = config["GeneticMap"]
popIDs = config["popIDs"]
outgrp = config["outgrp"]
rmPopID = config["rmPopID"]
sampleSize = config["sampleSize"]
segmentLength = config["segmentLength"]

MODELS = ["GutenkunstThreePop", "TennessenTwoPop", "RagsdaleArchaic"]

#-----------------------------------------------------------------------------
# Set final targets
#-----------------------------------------------------------------------------
rule all:
	wildcard_constraints:
		replicate = "\d+",
		pop = "[afr|eur]"
	input:
		expand("output/ArchIE/{model_name}/{model_name}_{mnm_status}_{pop}_predicted.txt", # , #, #% (modelNAME, modelNAME),
			model_name = MODELS,
			mnm_status=["mnm%s-%s" % (MNM_dist, MNM_frac), "womnm"],
			pop=["afr", "eur"]),
		# expand("output/ArchIE/{model_name}/{model_name}_womnm_{pop}_predicted.txt" , #, #% (modelNAME, modelNAME), pop=["afr", "eur"]),
		# expand("output/ArchIE/{model_name}/{model_name}_mnm%s-%s_{pop}_test_data.txt" % (modelNAME, modelNAME, MNM_dist, MNM_frac), pop=["afr", "eur"]),
		expand("output/ArchIE/{model_name}/{model_name}_{mnm_status}_{pop}_test_data.txt" , #, #% (modelNAME, modelNAME),
			model_name = MODELS,
			mnm_status=["mnm%s-%s" % (MNM_dist, MNM_frac), "womnm"],
			pop=["afr", "eur"]),
		sprime_plot = "output/sprime/%s/sinaplot_topSegment_wMNM_vs_womnm.out.score.pdf" , #% modelNAME
    	# "output/ArchIE/{model_name}/{model_name}_rep{replicate}_mnm%s-%s.txt" % (modelNAME, modelNAME, MNM_dist, MNM_frac),
    	# expand("output/ArchIE/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_afr.txt" % (modelNAME, modelNAME, MNM_dist, MNM_frac), replicate=range(1, replicates+1)),
    	# "output/ArchIE/{model_name}/{model_name}_mnm%s-%s_eur_predicted.txt" % (modelNAME, modelNAME, MNM_dist, MNM_frac),
    	# "output/ArchIE/{model_name}/{model_name}_womnm_eur_predicted.txt" , #, #% (modelNAME, modelNAME),


#-----------------------------------------------------------------------------
# Simulate and output in EIGENSTRAT format
#-----------------------------------------------------------------------------
rule sim_data_archie:
	# input:
	# 	outG = "output/msprime/{model_name}/{model_name}.%s.indID" % (modelNAME, modelNAME, outgrp)
	output:
		snp = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s.snp" % (MNM_dist, MNM_frac),
		geno_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_afr.geno" % (MNM_dist, MNM_frac),
		geno_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_eur.geno" % (MNM_dist, MNM_frac),
		ind_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_afr.ind" % (MNM_dist, MNM_frac),
		ind_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_eur.ind" % (MNM_dist, MNM_frac),
		snp_womnm = "output/msprime/{model_name}/{model_name}_rep{replicate}.snp" , #, #% (modelNAME, modelNAME),
		geno_afr_womnm = "output/msprime/{model_name}/{model_name}_rep{replicate}_afr.geno" , #, #% (modelNAME, modelNAME),
		geno_eur_womnm = "output/msprime/{model_name}/{model_name}_rep{replicate}_eur.geno" , #, #% (modelNAME, modelNAME),
		ind_afr_womnm = "output/msprime/{model_name}/{model_name}_rep{replicate}_afr.ind" , #, #% (modelNAME, modelNAME),
		ind_eur_womnm = "output/msprime/{model_name}/{model_name}_rep{replicate}_eur.ind" , #% (modelNAME, modelNAME)
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		shell(" python eval_models.py --replicates 1 --replicate_ID {wildcards.replicate} --length {segmentLength} --mnm_dist {MNM_dist} --mnm_frac {MNM_frac} --demographic_model {wildcards.model_name} --method archie")


#-----------------------------------------------------------------------------
# Run ArchIE on EUR data
# to-do: merge MNM/womnm rules by population
#-----------------------------------------------------------------------------
rule archie_stats_MNM_eur:
	input:
		snp = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s.snp" % (MNM_dist, MNM_frac),
		geno_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_eur.geno" % (MNM_dist, MNM_frac),
		geno_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_afr.geno" % (MNM_dist, MNM_frac),
		ind_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_eur.ind" % (MNM_dist, MNM_frac),
	output:
		archie_out_eur = "output/ArchIE/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_eur.txt" % (MNM_dist, MNM_frac)
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	wildcard_constraints:
		replicate="\d+"
	run:
		shell("python src/ArchIE/data/calc_stats_window_data.py -s {input.snp} -i {input.ind_eur} -a {input.geno_eur} -r {input.geno_afr} -c 1 -b 0 -e 50000 -w 50000 -z 50000 > {output.archie_out_eur} ")

rule archie_stats_womnm_eur:
	input:
		snp = "output/msprime/{model_name}/{model_name}_rep{replicate}.snp" , #, #% (modelNAME, modelNAME),
		geno_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_eur.geno" , #, #% (modelNAME, modelNAME),
		geno_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_afr.geno" , #, #% (modelNAME, modelNAME),
		ind_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_eur.ind" , #% (modelNAME, modelNAME)
	output:
		archie_out_eur = "output/ArchIE/{model_name}/{model_name}_rep{replicate}_womnm_eur.txt" , #% (modelNAME, modelNAME)
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	wildcard_constraints:
		replicate="\d+"
	run:
		shell("python src/ArchIE/data/calc_stats_window_data.py -s {input.snp} -i {input.ind_eur} -a {input.geno_eur} -r {input.geno_afr} -c 1 -b 0 -e 50000 -w 50000 -z 50000 > {output.archie_out_eur} ")


#-----------------------------------------------------------------------------
# Run ArchIE on AFR data
#-----------------------------------------------------------------------------
rule archie_stats_MNM_afr:
	input:
		snp = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s.snp" % (MNM_dist, MNM_frac),
		geno_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_afr.geno" % (MNM_dist, MNM_frac),
		geno_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_eur.geno" % (MNM_dist, MNM_frac),
		ind_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_afr.ind" % (MNM_dist, MNM_frac),
	output:
		archie_out_afr = "output/ArchIE/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_afr.txt" % (MNM_dist, MNM_frac),
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	wildcard_constraints:
		replicate="\d+"
	run:
		shell("python src/ArchIE/data/calc_stats_window_data.py -s {input.snp} -i {input.ind_afr} -a {input.geno_afr} -r {input.geno_eur} -c 1 -b 0 -e 50000 -w 50000 -z 50000 > {output.archie_out_afr}")

rule archie_stats_womnm_afr:
	input:
		snp = "output/msprime/{model_name}/{model_name}_rep{replicate}.snp" , #, #% (modelNAME, modelNAME),
		geno_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_afr.geno" , #, #% (modelNAME, modelNAME),
		geno_eur = "output/msprime/{model_name}/{model_name}_rep{replicate}_eur.geno" , #, #% (modelNAME, modelNAME),
		ind_afr = "output/msprime/{model_name}/{model_name}_rep{replicate}_afr.ind" , #% (modelNAME, modelNAME)
	output:
		archie_out_afr = "output/ArchIE/{model_name}/{model_name}_rep{replicate}_womnm_afr.txt" , #% (modelNAME, modelNAME)
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	wildcard_constraints:
		replicate="\d+"
	run:
		shell("python src/ArchIE/data/calc_stats_window_data.py -s {input.snp} -i {input.ind_afr} -a {input.geno_afr} -r {input.geno_eur} -c 1 -b 0 -e 50000 -w 50000 -z 50000 > {output.archie_out_afr}")


#-----------------------------------------------------------------------------
# merge Archie output
#-----------------------------------------------------------------------------
rule archie_merge:
	input:
		# "output/ArchIE/{model_name}/{model_name}_rep{replicate}_mnm%s-%s_{pop}.txt" % (modelNAME, modelNAME, MNM_dist, MNM_frac)
		expand("output/ArchIE/{model_name}/{model_name}_rep{replicate}_{mnm_status}_{pop}.txt" , #, #% (modelNAME, modelNAME),
			model_name = MODELS,
			replicate=range(1, replicates+1),
			mnm_status=["mnm%s-%s" % (MNM_dist, MNM_frac), "womnm"],
			pop=["afr", "eur"])
		# "output/ArchIE/{model_name}/{model_name}_rep{replicate}_{mnm_status}_{pop}.txt" , #% (modelNAME, modelNAME)
	output:
		"output/ArchIE/{model_name}/{model_name}_{mnm_status}_{pop}_test_data.txt" , #% (modelNAME, modelNAME)
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	wildcard_constraints:
		replicate = "\d+"
	run:
		shell("cat {input} > {output}")


#-----------------------------------------------------------------------------
#  create and train ArchIE model
#-----------------------------------------------------------------------------
rule archie_create_training:
	output:
		"src/ArchIE/simulations/training_data.txt"
		# "output/ArchIE/archie_trained_model.Rdata"
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		shell("bash src/ArchIE/simulations/create_training.sh")


rule archie_run_training:
	input:
		"src/ArchIE/simulations/training_data.txt"
	output:
		"output/ArchIE/archie_trained_model.Rdata"
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		shell("Rscript ArchIE_train.R {input} {output}")


#-----------------------------------------------------------------------------
#  Use ArchIE to predict archaic ancestry under each population with and without MNMs
#-----------------------------------------------------------------------------
rule archie_predict:
	input:
		training_data = "output/ArchIE/archie_trained_model.Rdata",
		testing_data = "output/ArchIE/{model_name}/{model_name}_{mnm_status}_{pop}_test_data.txt", #% (modelNAME, modelNAME)
	output:
		predicted_data = "output/ArchIE/{model_name}/{model_name}_{mnm_status}_{pop}_predicted.txt", #% (modelNAME, modelNAME)
	params:
		sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		shell("Rscript ArchIE_predict.R {input.training_data} {input.testing_data} {output.predicted_data}")


#-----------------------------------------------------------------------------
# Get sample IDs per group
#-----------------------------------------------------------------------------
rule output_IDfiles:
	output:
		expand("output/msprime/{model_name}/{model_name}.{ID}.indID",
			model_name = MODELS,
			ID=popIDs)
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		for outfile in output:
			ID = outfile.split(".")[1]
			with open(outfile, "w") as fout:
				for i in range(sampleSize):
					fout.write("%s_ind%s" % (ID, i) + "\n")


#-----------------------------------------------------------------------------
# Simulate and output in VCF format
#-----------------------------------------------------------------------------
rule sim_data_sprime:
	input:
		outG = "output/msprime/{model_name}/{model_name}.%s.indID" % outgrp
	output:
		vcf = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s.vcf" % (MNM_dist, MNM_frac),
		vcf_womnm = "output/msprime/{model_name}/{model_name}_rep{replicate}.vcf", #% (modelNAME, modelNAME)
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		if modelNAME == "GutenkunstThreePop":
			shell(" python eval_models.py --replicates 1 --replicate_ID {wildcards.replicate} --length {segmentLength} --mnm_dist {MNM_dist} --mnm_frac {MNM_frac} --demographic_model {wildcards.model_name} --method sprime")

		# else:
		# 	shell(" python eval_models.py {wildcards.replicate} {segmentLength} {MNM_dist} {MNM_frac} | \
		# 	    awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' | \
		# 	    bgzip -c > {output.vcf}")


#-----------------------------------------------------------------------------
# clean up non-MNM VCF
#-----------------------------------------------------------------------------
rule process_vcf_womnm:
    input:
        vcf = "output/msprime/{model_name}/{model_name}_rep{replicate}.vcf", #, #% (modelNAME, modelNAME),
        rm_ids = "output/msprime/{model_name}/{model_name}.%s.indID"  % rmPopID
    output:
        vcf = "output/msprime/{model_name}/{model_name}_rep{replicate}_womnm.vcf.gz" , #% (modelNAME, modelNAME)
    params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
    run:
        shell("""awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' {input.vcf} | \
			    vcftools --vcf - --remove {input.rm_ids} --recode --recode-INFO-all --stdout | \
			    bgzip -c > {output.vcf}""")


#-----------------------------------------------------------------------------
# clean up MNM VCF
#-----------------------------------------------------------------------------
rule process_vcf_wMNM:
    input:
        vcf = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s.vcf" % (MNM_dist, MNM_frac),
        rm_ids = "output/msprime/{model_name}/{model_name}.%s.indID"  % rmPopID
    output:
        vcf = "output/msprime/{model_name}/{model_name}_rep{replicate}_mnm%s-%s.vcf.gz" % (MNM_dist, MNM_frac)
    params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
    run:
        shell("""awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' {input.vcf} | \
			    vcftools --vcf - --remove {input.rm_ids} --recode --recode-INFO-all --stdout | \
			    bgzip -c > {output.vcf}""")


#-----------------------------------------------------------------------------
# run sprime on non-MNM data
#-----------------------------------------------------------------------------
rule sprime_run:
	input:
	    vcf = "output/msprime/{model_name}/{model_name}_rep{replicate}_{mnm_status}.vcf.gz" , #, #% (modelNAME, modelNAME),
		outG = "output/msprime/{model_name}/{model_name}.%s.indID" % (outgrp)
	output:
	    "output/sprime/{model_name}/{model_name}_rep{replicate}_{mnm_status}.sprime.out.score" , #% (modelNAME, modelNAME)
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		prefix = ".".join(output[0].split(".")[:-1])
		shell(" java -jar src/sprime.jar gt={input.vcf} outgroup={input.outG} map={GeneticMap} out={prefix} minscore=1 ")


rule sprime_pullTopSegment:
	input:
		"output/sprime/{model_name}/{model_name}_rep{replicate}_{mnm_status}.sprime.out.score" , #% (modelNAME, modelNAME)
	output:
		"output/sprime/{model_name}/{model_name}_rep{replicate}_{mnm_status}.sprime.out.score.top" , #% (modelNAME, modelNAME)
	params:
	    sge_opts="-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	shell:
		" set +o pipefail ; sort -n -r -k8 {input} | head -1 > {output} "


rule sprime_mergeTopSegment:
	input:
		# expand("output/sprime/{model_name}/{model_name}_rep{replicate}_{mnm_status}.sprime.out.score.top" , #, #% (modelNAME, modelNAME), replicate=range(1, replicates+1))
		expand("output/sprime/{model_name}/{model_name}_rep{replicate}_{mnm_status}.sprime.out.score.top" , #, #% (modelNAME, modelNAME),
			model_name = MODELS,
			replicate=range(1, replicates+1),
			mnm_status=["mnm%s-%s" % (MNM_dist, MNM_frac), "womnm"])
		# "output/sprime/{model_name}/{model_name}_rep{replicate}_{mnm_status}.sprime.out.score.top" , #% (modelNAME, modelNAME)
	output:
		"output/sprime/%s/topSegment_{mnm_status}.out.score" , #% modelNAME
	params:
	    sge_opts="-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
#	shell:
#		" cat {input} > {output} ; rm {input} "
	run:
		shell(""" cat {input} | grep -v "CHROM" > {output} """)
#		for f in input:
#			shell(" cat %s >> {output} " % (f) )
# 		shell(" echo output/sprime/{model_name}/{model_name}_womnm_rep{{1..%s}}.sprime.out.score.top | \
# 		    xargs cat | \
# 		    grep -v "CHROM" > {output} ; echo output/sprime/{model_name}/{model_name}_womnm_rep{{1..%s}}.sprime.out.score.top | xargs rm " % (modelNAME, modelNAME, replicates, modelNAME, modelNAME, replicates) )
# 		shell(" echo output/sprime/{model_name}/{model_name}_womnm_rep{{1..%s}}.sprime.out.score.top | \
# 		    xargs cat | \
# 		    grep -v "CHROM" > {output} ; echo output/sprime/{model_name}/{model_name}_womnm_rep{{1..%s}}.sprime.out.score.top | xargs rm " % (modelNAME, modelNAME, replicates, modelNAME, modelNAME, replicates)
		# shell(""" echo output/sprime/{model_name}/{model_name}_womnm_rep{{1..%s}}.sprime.out.score.top | \
		    # xargs cat | grep -v "CHROM" > {output}""" % (modelNAME, modelNAME, replicates))

rule sprime_plot:
	input:
	    top_womnm = "output/sprime/%s/topSegment_womnm.out.score" , #, #% modelNAME,
		top_mnm = "output/sprime/%s/topSegment_wMNM.out.score" , #% modelNAME
	output:
		"output/sprime/%s/sinaplot_topSegment_wMNM_vs_womnm.out.score.pdf" , #% modelNAME
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	shell:
		" Rscript plot_sprime_wMNM_vs_woMNM.r {input.top_mnm} {input.top_womnm} {output} "























