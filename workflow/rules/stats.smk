import sys
import os
import glob

#
# rules
#

rule get_fastq_stats:
    input:
        expand('{analysis_root}/{sample}/basecalling/{run_name}.fastq.gz', analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name'])
    output:
        expand("{analysis_root}/{sample}/basecalling/{run_name}_read_stats.tsv", analysis_root=config['analysis_root'], sample=config['sample'], run_name=config['run_name'])
    params:
        program=srcdir("../scripts/abyss-fac.pl")
    shell:
        "cat {input} | grep -Ev '^$' | awk 'NR % 4 == 1 || NR % 4 == 2' | tr '@' '>' | {params.program} > {output}"

