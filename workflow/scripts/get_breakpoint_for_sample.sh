#! /bin/bash

SAMPLE=$1
#cd $SAMPLE/mapped-pipeline
cd $SAMPLE
DIR=/.mounts/labs/simpsonlab/sw/createivf/bin
BP1=`python $DIR/get_breakpoint_regions.py $SAMPLE 1`
BP2=`python $DIR/get_breakpoint_regions.py $SAMPLE 2`

$DIR/select_reads.sh merged.sorted.bam $BP1 > bp1.reads
$DIR/select_reads.sh merged.sorted.bam $BP2 > bp2.reads
python3 $DIR/intersect_reads.py bp1.reads bp2.reads
