#####################################################
## samples.tsv
# sample	assembly	descriptive_name	cell
# CMC_ATAC	hg38	CMC_ATAC	CMC
# CMC_H3K27ac	hg38	CMC_H3K27ac	CMC
# KRT_p300	hg38	KRT_p300	KRT


#####################################################
## Snakefile

import os
import pandas as pd
configfile: "config.yaml"

# do onstart/onexit things
sample_schemas = ['sample', 'assembly', 'strandedness', 'condition']
config_schemas = ['general', 'download', 'alignment_general', 'peakcalling', 'trackhub']
include: "../../rules/configuration.smk"

# load all the relevant rules
include: f"{config['rule_dir']}/alignment.smk"
include: f"{config['rule_dir']}/call_peak.smk"
include: f"{config['rule_dir']}/get_genome.smk"
include: f"{config['rule_dir']}/get_fastq.smk"
include: f"{config['rule_dir']}/merge_replicates.smk"
include: f"{config['rule_dir']}/peak_count.smk"
include: f"{config['rule_dir']}/qc.smk"
include: f"{config['rule_dir']}/trackhub.smk"
include: f"{config['rule_dir']}/trimming.smk"
include: f"{config['rule_dir']}/pananse.smk"


# set the quality_control functions
quality_control = [get_trimming_qc, get_alignment_qc, get_peak_calling_qc]


rule call_enahncer_all:
    """
    call enahncers for each sample (or condition if applies)
    """
    input:
        expand((["{result_dir}/trackhub"] if config['create_trackhub'] else []) +
               (["{qc_dir}/multiqc_{assemblies}.html"] if config["create_qc_report"] else []) +
               ["{result_dir}/enhancer/{assemblies}-{cell}_enahncer.bed"],
                **{**config,
                   **{'assemblies': set(samples['assembly']), 'peak_caller': config['peak_caller'].keys(), 'cell': set(samples.cell)}})


#####################################################
## Code in pananse.smk


import os

## How can I merge sort_bdg_p300 and sort_bdg_h2k27ac
rule sort_bdg_p300:
    """
    Sort bdg file from macs2.
    """
    input:
        bdgfile = expand("{result_dir}/macs2/{{assembly}}-{{cell}}_p300_treat_pileup.bdg", **config)
    output:
        sortbdgfile = expand("{result_dir}/tmp/{{assembly}}-{{cell}}_sort.bdg", **config)
    shell:
        "sort -k1,1 -k2,2n {input.bdgfile} > {output.sortbdgfile}"

rule sort_bdg_h2k27ac:
    """
    Sort bdg file from macs2.
    """
    input:
        bdgfile = expand("{result_dir}/macs2/{{assembly}}-{{cell}}_H3K27ac_treat_pileup.bdg", **config)
    output:
        sortbdgfile = expand("{result_dir}/tmp/{{assembly}}-{{cell}}_sort.bdg", **config)
    shell:
        "sort -k1,1 -k2,2n {input.bdgfile} > {output.sortbdgfile}"


rule bdg2wig:
    """
    Switch bdg file to wig file.
    """
    input:
        sortbdgfile = expand("{result_dir}/tmp/{{assembly}}-{{cell}}_sort.bdg", **config),
        genome_size = expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        wigfile = expand("{result_dir}/tmp/{{assembly}}-{{cell}}.wig", **config)
    conda:
        "../envs/pananse.yaml"
    shell:
        "bedGraphToBigWig {input.sortbdgfile} {input.genome_size} {output.wigfile}"


rule rename_peak_atac:
    input:
        narrowpeak_atac = expand("{result_dir}/macs2/{{assembly}}-{{cell}}_ATAC_peaks.narrowPeak", **config)
    output:
        narrowpeak = expand("{result_dir}/tmp/{{assembly}}-{{cell}}_peaks.narrowPeak", **config)
    shell:
        "mv {input.narrowpeak_atac} {output.narrowpeak}"


rule rename_peak_p300:
    input:
        narrowpeak_p300 = expand("{result_dir}/macs2/{{assembly}}-{{cell}}_p300_peaks.narrowPeak", **config)
    output:
        narrowpeak = expand("{result_dir}/tmp/{{assembly}}-{{cell}}_peaks.narrowPeak", **config)
    shell:
        "mv {input.narrowpeak_p300} {output.narrowpeak}"


def read_chrsize(sizefile):
    sdic = {}
    with open (sizefile) as sizes:
        for chrom in sizes:
            sdic[chrom.split()[0]] = int(chrom.split()[1])
    return sdic


rule make_enhancer:
    input:
        narrowpeak = expand("{result_dir}/tmp/{{assembly}}-{{cell}}_peaks.narrowPeak", **config),
        wigfile = expand("{result_dir}/tmp/{{assembly}}-{{cell}}.wig", **config),
        genome_size = expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        enhancerbed = expand("{result_dir}/enhancer/{{assembly}}-{{cell}}_enahncer.bed", **config)
    # conda:
    #     "../envs/pananse.yaml"
    run:
        chrsizedic = read_chrsize(input.genome_size[0])
        with open(str(input.narrowpeak)) as bdgf, open(str(output.enhancerbed),"w") as enh_bed:
            for line in bdgf:
                a = line.split()
                peak = int(a[9]) + int(a[1])
                start = peak - 100
                end = peak + 100
                if start>0 and chrsizedic[a[0]]>end:
                    commd = "/home/qxu/bin/ucsc/bigWigSummary -type=max {wig} {chrname} {start} {end} 1".format(wig=str(input.wigfile), chrname=a[0], start=start, end=end)
                    commd_result = os.popen(commd)
                    r = commd_result.read()
                    if r != "":
                        enh_bed.write("{chrname}\t{start}\t{end}\t{score}".format(chrname=a[0], start=start, end=end, score=str(r)))

