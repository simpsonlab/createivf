import sys
import os
import glob

#
# rules
#
rule run_guppy_on_dir:
    input:
        get_fast5_dir
    output:
        pass_dir=directory("{analysis_root}/{sample}/basecalling/{run_name}/pass")
    params:
        program=get_basecaller,
        config=get_guppy_config,
        save_path=get_save_path,
        num_callers=get_guppy_num_callers,
        gpu_runners_per_device=get_gpu_runners_per_device,
        chunks_per_runner=get_chunks_per_runner,
        device=get_gpu_device
    shell:
        "{params.program} -c {params.config} -i {input} -s {params.save_path} --num_callers {params.num_callers} --gpu_runners_per_device {params.gpu_runners_per_device} --chunks_per_runner {params.chunks_per_runner} -x {params.device} --disable_pings --compress_fastq"
