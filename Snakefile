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


rule dummy:
	input:
		"output/sprime/%s/sinaplot_topSegment_wMNM_vs_woMNM.out.score.pdf" % modelNAME


rule plot:
	input:
	    top_womnm = "output/sprime/%s/topSegment_woMNM.out.score" % modelNAME,
		top_mnm = "output/sprime/%s/topSegment_wMNM.out.score" % modelNAME
	output:
		"output/sprime/%s/sinaplot_topSegment_wMNM_vs_woMNM.out.score.pdf" % modelNAME
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	shell:
		""" Rscript plot_sprime_wMNM_vs_woMNM.r {input.top_mnm} {input.top_womnm} {output} """


rule mergeTopSegment_woMNM:
	input:
		expand("output/sprime/%s/%s_womnm_rep{replicate}.sprime.out.score.top" % (modelNAME, modelNAME), replicate=range(1, replicates+1))
	output:
		"output/sprime/%s/topSegment_woMNM.out.score" % modelNAME
	params:
	    sge_opts="-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
#	shell:
#		""" cat {input} > {output} ; rm {input} """
	run:
#		for f in input:
#			shell(""" cat %s >> {output} """ % (f) )
# 		shell(""" echo output/sprime/%s/%s_womnm_rep{{1..%s}}.sprime.out.score.top | \
# 		    xargs cat | \
# 		    grep -v "CHROM" > {output} ; echo output/sprime/%s/%s_womnm_rep{{1..%s}}.sprime.out.score.top | xargs rm """ % (modelNAME, modelNAME, replicates, modelNAME, modelNAME, replicates) )
# 		shell(""" echo output/sprime/%s/%s_womnm_rep{{1..%s}}.sprime.out.score.top | \
# 		    xargs cat | \
# 		    grep -v "CHROM" > {output} ; echo output/sprime/%s/%s_womnm_rep{{1..%s}}.sprime.out.score.top | xargs rm """ % (modelNAME, modelNAME, replicates, modelNAME, modelNAME, replicates)
		shell(""" echo output/sprime/%s/%s_womnm_rep{{1..%s}}.sprime.out.score.top | \
		    xargs cat | grep -v "CHROM" > {output}""" % (modelNAME, modelNAME, replicates))


rule pullTopSegment_woMNM:
	input:
		"output/sprime/%s/%s_womnm_rep{replicate}.sprime.out.score" % (modelNAME, modelNAME)
	output:
		"output/sprime/%s/%s_womnm_rep{replicate}.sprime.out.score.top" % (modelNAME, modelNAME)
	params:
	    sge_opts="-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	shell:
		""" set +o pipefail ; sort -n -r -k8 {input} | head -1 > {output} """


rule mergeTopSegment_wMNM:
	input:
		expand("output/sprime/%s/%s_rep{replicate}_mnm%s-%s.sprime.out.score.top" % (modelNAME, modelNAME, MNM_dist, MNM_frac), replicate=range(1, replicates+1))
	output:
		"output/sprime/%s/topSegment_wMNM.out.score" % modelNAME
	params:
	    sge_opts="-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
# 		shell(""" echo output/sprime/%s/%s_rep{{1..%s}}_mnm%s-%s.sprime.out.score.top | \
# 		    xargs cat | \
# 		    grep -v "CHROM" > {output} ; echo output/sprime/%s/%s_rep{{1..%s}}_mnm%s-%s.sprime.out.score.top | \
# 		    xargs rm """ % (modelNAME, modelNAME, replicates, MNM_dist, MNM_frac, modelNAME, modelNAME, replicates, MNM_dist, MNM_frac) )
		shell(""" echo output/sprime/%s/%s_rep{{1..%s}}_mnm%s-%s.sprime.out.score.top | \
	    xargs cat | \
	    grep -v "CHROM" > {output}""" % (modelNAME, modelNAME, replicates, MNM_dist, MNM_frac) )


rule pullTopSegment_wMNM:
	input:
		"output/sprime/%s/%s_rep{replicate}_mnm%s-%s.sprime.out.score" % (modelNAME, modelNAME, MNM_dist, MNM_frac)
	output:
		"output/sprime/%s/%s_rep{replicate}_mnm%s-%s.sprime.out.score.top" % (modelNAME, modelNAME, MNM_dist, MNM_frac)
	params:
	    sge_opts="-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	shell:
		""" set +o pipefail ; sort -n -r -k8 {input} | head -1 > {output} """

#-----------------------------------------------------------------------------
# run sprime on non-MNM data
#-----------------------------------------------------------------------------
rule sprime_woMNM:
	input:
	    vcf = "output/msprime/%s/%s_womnm_rep{replicate}.vcf.gz" % (modelNAME, modelNAME),
		outG = "output/msprime/%s/%s.%s.indID" % (modelNAME, modelNAME, outgrp)
	output:
	    "output/sprime/%s/%s_womnm_rep{replicate}.sprime.out.score" % (modelNAME, modelNAME)
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		prefix = ".".join(output[0].split(".")[:-1])
		shell(""" java -jar src/sprime/sprime.jar gt={input.vcf} outgroup={input.outG} map={GeneticMap} out={prefix} minscore=1 """)


#-----------------------------------------------------------------------------
# run sprime on MNM data
#-----------------------------------------------------------------------------
rule sprime_wMNM:
	input:
	    vcf = "output/msprime/%s/%s_rep{replicate}_mnm%s-%s.vcf.gz" % (modelNAME, modelNAME, MNM_dist, MNM_frac),
		outG = "output/msprime/%s/%s.%s.indID" % (modelNAME, modelNAME, outgrp)
	output:
	    "output/sprime/%s/%s_rep{replicate}_mnm%s-%s.sprime.out.score" % (modelNAME, modelNAME, MNM_dist, MNM_frac)
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		prefix = ".".join(output[0].split(".")[:-1])
		shell(""" java -jar src/sprime/sprime.jar gt={input.vcf} outgroup={input.outG} map={GeneticMap} out={prefix} minscore=1 """)


#-----------------------------------------------------------------------------
# clean up non-MNM VCF
#-----------------------------------------------------------------------------
rule process_vcf_woMNM:
    input:
        vcf = "output/msprime/%s/%s_rep{replicate}.vcf" % (modelNAME, modelNAME),
        rm_ids = "output/msprime/%s/%s.%s.indID"  % (modelNAME, modelNAME, rmPopID)
    output:
        vcf = "output/msprime/%s/%s_womnm_rep{replicate}.vcf.gz" % (modelNAME, modelNAME)
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
        vcf = "output/msprime/%s/%s_rep{replicate}_mnm%s-%s.vcf" % (modelNAME, modelNAME, MNM_dist, MNM_frac),
        rm_ids = "output/msprime/%s/%s.%s.indID"  % (modelNAME, modelNAME, rmPopID)
    output:
        vcf = "output/msprime/%s/%s_rep{replicate}_mnm%s-%s.vcf.gz" % (modelNAME, modelNAME, MNM_dist, MNM_frac)
    params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
    run:
        shell("""awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' {input.vcf} | \
			    vcftools --vcf - --remove {input.rm_ids} --recode --recode-INFO-all --stdout | \
			    bgzip -c > {output.vcf}""")


#-----------------------------------------------------------------------------
# Simulate data using msprime
#
#-----------------------------------------------------------------------------
rule sim_data:
	input:
		outG = "output/msprime/%s/%s.%s.indID" % (modelNAME, modelNAME, outgrp)
	output:
		vcf = "output/msprime/%s/%s_rep{replicate}_mnm%s-%s.vcf" % (modelNAME, modelNAME, MNM_dist, MNM_frac),
		vcf_woMNM = "output/msprime/%s/%s_rep{replicate}.vcf" % (modelNAME, modelNAME)
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		if modelNAME == "GutenkunstThreePop":
			shell(""" python eval_models.py --replicates 1 --replicate_ID {wildcards.replicate} --length {segmentLength} --mnm_dist {MNM_dist} --mnm_frac {MNM_frac} --demographic_model {modelNAME} --method sprime""")
# 		else:
# 			shell(""" python eval_models.py {wildcards.replicate} {segmentLength} {MNM_dist} {MNM_frac} | \
# 			    awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' | \
# 			    bgzip -c > {output.vcf}""")

#

#-----------------------------------------------------------------------------
# Get sample IDs per group
#-----------------------------------------------------------------------------
rule output_IDfiles:
	output:
		expand("output/msprime/%s/%s.{ID}.indID" % (modelNAME, modelNAME), ID=popIDs)
	params:
	    sge_opts = "-l h_rt=120:00:00 -l mfree=4G -l gpfsstate=0"
	run:
		for ID in popIDs:
			with open("output/msprime/%s/%s.%s.indID" % (modelNAME, modelNAME, ID), "w") as fout:
				for i in range(sampleSize):
					fout.write("%s_ind%s" % (ID, i) + "\n")

