import sys
import os
import glob

#
# rules
#
rule merge_index_fastq_files:
    input:
        get_run_fastq_files
    output:
        expand('{analysis_root}/{sample}/basecalling/{run_name}.fastq.gz', analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name'])
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
        postprocess="samtools sort",
        threads=workflow.cores,
        preset="map-ont"
    shell:
        "{params.aligner} -t {params.threads} -ax {params.preset} {input.ref} {input.fastq_files} | {params.postprocess} -T {output}.tmp -o {output} -"

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
