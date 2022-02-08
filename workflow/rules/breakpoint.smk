import sys
import os
import glob

#
# rules
#
rule get_breakpoints:
    input:
        bam=expand("{analysis_root}/{sample}/merged.sorted.bam", analysis_root=config['analysis_root'], sample=config['sample']), 
        bai=expand("{analysis_root}/{sample}/merged.sorted.bam.bai", analysis_root=config['analysis_root'], sample=config['sample'])
    output:
        "{sample}/bp{bp_id}.reads"
    params:
        sample=config['sample'],
        program=srcdir('../scripts/createivf_breakpoints.py'),
        cytobands=get_cytobands,
        metadata=get_metadata
    shell:
        "python {params.program} --sample {params.sample} --metadata {params.metadata} --cytobands {params.cytobands} --bam {input.bam} --outprefix '{params.sample}/bp'"

