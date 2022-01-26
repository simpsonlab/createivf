import sys
import os
import glob

configfile: "config.yaml"

#
# helper functions
#
def get_dir_index_pattern():
    """
    Return the index directory pattern
    """
    return config.get('index_pattern', '{run_root}/{run_name}/fast5/{index}')

def get_index_dir(wildcards):
    """
    Return the FAST5 directory for a given index
    """
    index_pattern = get_dir_index_pattern()
    index_dir = index_pattern.format(run_root=config['run_root'], run_name=config['run_name'], index=wildcards.index)
    return index_dir

def get_save_path_pattern():
    """
    Return basecall save path pattern
    """
    return config.get('save_pattern', '{analysis_root}/{sample}/basecalling/{run_name}/{index}-basecalled')

def get_save_path(wildcards):
    save_path_pattern = get_save_path_pattern()
    save_path_dir = save_path_pattern.format(analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name'], index=wildcards.index)

def get_basecall_output(wildcards):
    """
    Return a list of indexed basecall output directories
    """
    basecall_dirs = list()
    index_dirs_path = '/'.join([config['run_root'], config['run_name'], 'fast5', '*'])
    for index_path in glob.glob(index_dirs_path):
        index_dir = os.path.basename(index_path)
        if os.path.isdir(index_path) and index_dir.isnumeric():
            basecall_dirs.append('/'.join([config['analysis_root'], config['sample'], 'basecalling', config['run_name'], '-'.join([index_dir, 'basecalled'])])) 
    return basecall_dirs

def get_index_fastq_pattern(wildcards):
    """
    Return the basecalled FASTQ files generated from the indices in the
    fast5 directory
    """
    return "{analysis_root}/{sample}/basecalling/{run_name}/{index}-basecalled/pass/{id}.fastq"

def get_index_fastq(wildcards):
    """
    Return a list of basecalled FASTQ files from each of the indexed directories
    """
    index_fastq_files = list()
    index_fastq_paths = glob.glob('/'.join([config['analysis_root'], config['sample'], 'basecalling', config['run_name'], '*', 'pass']))
    for path in index_fastq_paths:
        if os.path.isdir(path):
            for file in os.listdir(path):
                if file.endswith('.fastq'):
                    index_fastq_files.append('/'.join([path, file]))
    return index_fastq_files

def get_sample_fastq(wildcards):
    """
    Return the full path and filename for the merged FASTQ file for a given sample
    from a run
    """
    sample_fastq_file = '.'.join([config['sample'], 'fastq'])
    sample_fastq_path = '/'.join([config['analysis_root'], config['sample'], config['run_name'], sample_fastq_file])
    return sample_fastq_path


def get_reference(wildcards):
    """
    Return the path to the reference genome (this should be in the config.yaml)
    """
    return config['reference']

def get_sample_run_bam_files(wildcards):
    """
    Return a list of sample BAM files for a given run
    """
    bam_files = glob.glob('/'.join([config['analysis_root'], config['sample'], 'mapped-pipeline', '*.minimap.sorted.bam']))
    return bam_files

def get_merged_sample_bam_file(wildcards):
    """
    Return the merged BAM file across all runs for a given sample.
    """
    return '/'.join([config['analysis_root'], config['sample'], 'merged.sorted.bam'])

#
# rules
#
rule all:
    input: get_basecall_output

rule all_map:
    input:
        expand("{analysis_root}/{sample}/mapped-pipeline/{run_name}.minimap.sorted.bam.bai", analysis_root=config['analysis_root'], run_name=config['run_name'], sample=config['sample'])

rule all_merge:
    input:
        expand("{analysis_root}/{sample}/merged.sorted.bam.bai", analysis_root=config['analysis_root'], sample=config['sample'])

rule all_stats:
    input: 
        expand("{analysis_root}/{sample}/basecalling/{run_name}_read_stats.tsv", analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name'])

#
# dependency rules
#        
rule run_guppy_on_dir:
    input:
        get_index_dir
    output:
        directory("{analysis_root}/{sample}/basecalling/{run_name}/{index}-basecalled")
        #get_basecall_output
    params:
        program="guppy_basecaller",
        config="dna_r9.4.1_450bps_hac.cfg",
        save_path=get_save_path,
        num_callers="8",
        gpu_runner_per_device="4",
        chunks_per_runner="512",
        device="cuda:0 cuda:1"
    shell:
        "{params.program} -c {params.config} -i {input} -s {output} --num_callers {params.num_callers} --chunks_per_runner {params.chunks_per_runner} -x {params.device} --disable_pings"

rule merge_index_fastq_files:
    input:
        get_index_fastq
    output:
        expand('{analysis_root}/{sample}/basecalling/{run_name}.fastq', analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name'])
    params:
        program="cat"
    shell:
        "{params.program} {input} > {output}"

rule map_sample_run_fastq_file:
    input:
        fastq_files=expand('{analysis_root}/{sample}/basecalling/{run_name}.fastq', analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name']),
        ref=get_reference
    output:
        expand("{analysis_root}/{sample}/mapped-pipeline/{run_name}.minimap.sorted.bam", analysis_root=config['analysis_root'], run_name=config['run_name'], sample=config['sample'])
    params:
        aligner="minimap2",
        postprocess="samtools sort"
    shell:
        "{params.aligner} -t 8 -ax map-ont {input.ref} {input.fastq_files} | {params.postprocess} -T {output}.tmp -o {output} -"

rule index_mapped_bam:
    input:
        expand("{analysis_root}/{sample}/mapped-pipeline/{run_name}.minimap.sorted.bam", analysis_root=config['analysis_root'], run_name=config['run_name'], sample=config['sample'])
    output:
        expand("{analysis_root}/{sample}/mapped-pipeline/{run_name}.minimap.sorted.bam.bai", analysis_root=config['analysis_root'], run_name=config['run_name'], sample=config['sample'])
    params:
        program="samtools index"
    shell:
        "{params.program} {input}"

rule merge_run_sample_bam_files:
    input:
        get_sample_run_bam_files        
    output:
        expand("{analysis_root}/{sample}/merged.sorted.bam", analysis_root=config['analysis_root'], sample=config['sample'])
    params:
        program='samtools merge'
    shell:
        "{params.program} {output} {input}"

rule index_merged_bam_file:
    input:
        get_merged_sample_bam_file
    output:
        expand("{analysis_root}/{sample}/merged.sorted.bam.bai", analysis_root=config['analysis_root'], sample=config['sample'])
    params:
        program="samtools index"
    shell:
        "{params.program} {input}"

rule get_fastq_stats:
    input:
        get_index_fastq
    output:
        expand("{analysis_root}/{sample}/basecalling/{run_name}_read_stats.tsv", analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name'])
    params:
        program=srcdir("../scripts/abyss-fac.pl")
    shell:
        "cat {input} | grep -Ev '^$' | awk 'NR % 4 == 1 || NR % 4 == 2' | tr '@' '>' | {params.program} > {output}"

