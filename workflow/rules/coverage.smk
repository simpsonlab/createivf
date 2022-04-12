import os
import sys

configfile: "config.yaml"

# helper functions
def get_sample():
    """
    Return the sample name from the config.yaml
    """
    return config['sample']

def get_merged_bam(wildcards):
    """
    Return the full path to the final merged BAM file for a sample
    """
    return f'{get_sample()}/merged.sorted.bam'

def get_depth_files(wildcards):
    """
    Return a list of depth of coverage files
    """
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    depth_files = list()
    for chrom in chromosomes:
        depth_files.append(f'{get_sample()}/coverage/{get_sample()}.{chrom}.tsv')
    return depth_files

def get_nanoplot_stats(wildcards):
    """
    Return the NanoPlot stats file
    """
    return f'{get_sample()}/coverage/nanoplot/{get_sample()}NanoStats.txt'

def get_coverage_by_window_files(wildcards):
    """
    Return a list containing the coverage.txt files
    """
    coverage_files = list()
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    for chrom in chromosomes:
        coverage_files.append(f'{get_sample()}/coverage/{get_sample()}.{chrom}.coverage.txt')
    return coverage_files

def get_concatenated_coverage_file(wildcards):
    """
    Return the concatenated coverage file
    """
    return f'{get_sample()}/coverage/{get_sample()}.depth.txt'

#
# rules for coverage analysis
#

rule run_nanoplot:
    input:
        bam="{sample}/merged.sorted.bam"
    output:
        protected("{sample}/coverage/nanoplot/{sample}NanoPlot-report.html"),
        protected("{sample}/coverage/nanoplot/{sample}NanoStats.txt")
    params:
        program='NanoPlot',
        outdir="{sample}/coverage/nanoplot"
    shell:
        "{params.program} -t 8 -p {wildcards.sample} --title {wildcards.sample} --bam {input.bam} -o {params.outdir}"

rule run_samtools_depth:
    input:
        bam="{sample}/merged.sorted.bam"
    output:
        "{sample}/coverage/{sample}.{chromosomes}.tsv"
    params:
        program='samtools depth'
    shell:
        "{params.program} -r {wildcards.chromosomes} -a {input.bam} > {output}"

rule run_coverage_window:
    input:
        depth='{sample}/coverage/{sample}.{chromosomes}.tsv',
        stats='{sample}/coverage/nanoplot/{sample}NanoStats.txt'
    output:
        protected('{sample}/coverage/{sample}.{chromosomes}.coverage.txt')
    params:
        program=srcdir('../scripts/coverage_by_window.py')
    shell:
        'python {params.program} --input {input.depth} --output {output} --nanostats {input.stats}'

rule concatenate_coverage_files:
    input:
        coverage_files=get_coverage_by_window_files
    output:
        protected("{sample}/coverage/{sample}.depth.txt")
    params:
        program='cat'
    shell:
        '{params.program} {input.coverage_files} > {output}'

