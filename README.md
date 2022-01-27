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
```

### Run the basecaller
The basecaller uses `guppy` and a GPU to convert FAST5 files to FASTQ files.
```
snakemake -s /path/to/createivf/workflow/Snakefile --cores 1 all_basecall
```

Once basecalling has completed, the remainder of the pipeline can
be executed on a standard server:
```
snakemake -s /path/to/createivf/workflow/Snakefile --cores 8 all_breakpoint
```

## Credits and Acknowledgements
The script (`abyss-fac.pl`) to generate read stats was obtained from:
```
https://github.com/bcgsc/abyss/blob/master/bin/abyss-fac.pl
```


## License
MIT
