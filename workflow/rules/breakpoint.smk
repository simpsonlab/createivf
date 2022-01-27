import sys
import os
import glob

#
# rules
#
rule get_breakpoints:
    input:
        expand("{analysis_root}/{sample}/merged.sorted.bam.bai", analysis_root=config['analysis_root'], sample=config['sample'])
    output:
        "{sample}/bp{bp_id}.reads"
    params:
        sample=config['sample'],
        program=srcdir('../scripts/get_breakpoint_for_sample.sh')
    shell:
        "bash {params.program} {params.sample}"

