# createivf

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The `createivf` uses `snakemake` to automate the generation of
breakpoint read support from Nanopore sequenced reads.

## Installation
Download the package:
```
git clone https://github.com/simpsonlab/createivf.git
```

The pipeline uses `conda` to manage dependencies:
```
conda env create -f createivf/workflows/envs/environment.yaml
```

Once the `conda` dependencies have completed installing, activate
the conda environment:
```
conda activate createivf
```

Basecalling is performed with `guppy`.  A Nvidia GPU is required
on the server to speed up the process of generating FASTQ files
using the CUDA set of libraries.


## Usage
To run the workflow, individual run steps should be executed followed
by running breakpoints which occur on the merged run data for a
single sample.

### Configuration
To configure the pipeline, a `config.yaml` file needs to be created with
various parameters defined:
```
run_name: "run_name"
run_root: "/path/to/ont/runs"
index_pattern: "{run_root}/{run_name}/fast5/{index}"
analysis_root: "/path/to/analysis/root/directory"
sample: "samplename"
reference: "/path/to/reference/genome.fa"
basecaller: "guppy_basecaller"
guppy_config: "dna_r9.4.1_450bps_hac.cfg"
num_callers: "8"
gpu_runner_per_device: "4"
chunks_per_runner: "512"
device: "'cuda:0 cuda:1'"
metadata: "/path/to/sample_cytogenetics.tsv"
cytobands: "/path/to/hg38.cytoBand.composite.txt"
```

The `guppy_basecaller` is not available as part of the `conda` package
and must be installed separately.  The following parameters are
set for the GPU version of `guppy`:
* `basecaller`
* `guppy_config`
* `num_callers`
* `gpu_runner_per_device`
* `chunks_per_runner`
* `device`

There are two required files for running breakpoint analysis:
* `metadata` which is a tab separated file containing `sample`, band for `region1` and band for `region2`
* `cytobands` which contains the genomic regions and their corresponding cytogenetic bands

### Run the workflow
The basecaller uses `guppy` and a GPU to convert FAST5 files to FASTQ files.
```
snakemake -s /path/to/createivf/workflow/Snakefile --cores 1 all_basecall
```

Once basecalling has completed, the remainder of the pipeline can
be executed on a standard server.  Individual runs can be executed
using the following to generate BAM files for each run for a given
sample:
```
snakemake -s /path/to/createivf/workflow/Snakefile --cores 8 all_map
```
Once alignment has completed, the `all_breakpoint` rule will execute
the merging of individual runs per sample into a single sample
BAM file and run breakpoint detection on the merged file:
```
snakemake -s /path/to/createivf/workflow/Snakefile --cores 8 all_breakpoint
```
This will merge each run BAM file for a given sample into a single sample
BAM file and breakpoint read information:
```
sample/merged.sorted.bam
sample/bp1.reads
sample/bp2.reads
```


## Credits and Acknowledgements
The script (`abyss-fac.pl`) to generate read stats was obtained from:
```
https://github.com/bcgsc/abyss/blob/master/bin/abyss-fac.pl
```
The cytogenetic band file was obtained from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz).  A column was appended to the file that prefixed the band with `chr`.

## License
MIT
