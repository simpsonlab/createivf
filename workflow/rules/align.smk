
configfile: "config.yaml"

import os

#
# helper functions
#
def get_fastq_file():
    """
    Return the full path to the FASTQ file to process
    """
    pass

def get_reference_fasta():
    """
    Return the reference genome from the config.yaml file
    """
    return config.get('reference', '')


def get_minimap():
    """
    Return the aligner to use
    """
    return config['aligner']


def get_align_threads():
    """
    Return the number of threads to use (default: 8)
    """
    return config.get('aligner_threads', '8')

def get_minimap_preset():
    """
    Return the presets used by the aligner.
    """
    return config.get('minimap_preset', 'map-ont') 


#
# rules
#
rule minimap2_align:
# $MINIMAP2_ROOT/bin/minimap2 -t 8 -ax map-ont /.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa $IN | samtools sort -T $OUT.tmp -o $OUT -
    input:
        fastq=get_fastq_file,
        reference=get_reference_fasta
    output:
        sam=get_mapped_sam
    params:
        program=get_minimap,
        threads=get_align_threads,
        preset=get_minimap_preset 
    shell:
        "{params.program} -t {params.threads} -a -x {params.preset} {input.reference} {input.fastq}"
