import sys
import os
import glob

#
# rules
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

