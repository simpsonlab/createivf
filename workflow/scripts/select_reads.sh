#! /bin/bash -x
INPUT=$1
REGION=$2

#DIR=/.mounts/labs/ont/createivf/pipeline/bin
DIR=/.mounts/labs/simpsonlab/sw/createivf/bin
samtools view -h $INPUT $REGION | python $DIR/alignment_length.py | sort | uniq
