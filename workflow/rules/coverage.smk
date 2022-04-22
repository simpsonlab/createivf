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
        coverage_files.append(f'{get_sample()}/coverage/{get_sample()}.{chrom}.coverage.txt.gz')
    return coverage_files

def get_concatenated_coverage_file(wildcards):
    """
    Return the concatenated coverage file
    """
    return f'{get_sample()}/coverage/{get_sample()}.depth.txt.gz'

def get_bed_option(wildcards):
    """
    Return a BED file which lists regions to include.
    """
    if os.path.isfile(config['bed']):
        return f' -b {config["bed"]}'
    else:
        return ""

def get_gaps_file(wildcards):
    """
    Return a BED file containing gaps from the config
    """
    if os.path.isfile(config['gaps']):
        return config['gaps']
    else:
        print(f'Missing gaps file')
        sys.exit(1)

def get_size_file(wildcards):
    """
    Return the chromosome sizes file from the config
    """
    if os.path.isfile(config['size']):
        return config['size']
    else:
        print(f'Missing size file')
        sys.exit(1)

def get_coverage(sample_name):
    """
    Get the genome coverage from the NanoStats.txt file
    """
    file_name = str(sample_name) + "/coverage/nanoplot/" + str(sample_name) + "NanoStats.txt"
    if os.path.exists(file_name):
        f = open(file_name)
        pattern = "Total bases"
        for line in f:
            if re.search(pattern, line):
                total_bases = line.split(":")[1].strip().replace(',', '')
                return float(total_bases)/3100000000

def get_copy_number_bed(wildcards):
    """
    Return the copy number variant BED file
    """
    return f'{get_sample()}/{get_sample()}.CNVs.bed'

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
        "{sample}/coverage/{sample}.{chromosomes}.tsv.gz"
    params:
        program='samtools depth',
        bed=get_bed_option
    shell:
        "{params.program}  -r {wildcards.chromosomes} {params.bed} -a {input.bam} | gzip -c - > {output}"

rule run_coverage_window:
    input:
        depth='{sample}/coverage/{sample}.{chromosomes}.tsv.gz',
        stats='{sample}/coverage/nanoplot/{sample}NanoStats.txt'
    output:
        protected('{sample}/coverage/{sample}.{chromosomes}.coverage.txt.gz')
    params:
        program=srcdir('../scripts/coverage_by_window.py')
    shell:
        'python {params.program} --input {input.depth} --output {output} --nanostats {input.stats}'

rule concatenate_coverage_files:
    input:
        coverage_files=get_coverage_by_window_files
    output:
        protected("{sample}/coverage/{sample}.depth.txt.gz")
    params:
        program='cat'
    shell:
        '{params.program} {input.coverage_files} > {output}'

rule get_genome_copy_number:
    input:
        stats='{sample}/coverage/nanoplot/{sample}NanoStats.txt',
        depth='{sample}/coverage/{sample}.depth.txt.gz'
    output:
        "{sample}/{sample}.CNVs.bed"
    params:
        program=srcdir('../scripts/filter_CNs.py'),
        gaps=get_gaps_file,
        size=get_size_file,
        cov=lambda wildcards: get_coverage(wildcards.sample)
    shell:
        'python {params.program} -depth {input.depth} -gaps {params.gaps} -cnv_out {output} -size {params.size} -cov {params.cov}'
        
