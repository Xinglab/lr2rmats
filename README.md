# lr2gtf: Long read to GTF
## <a name="lr2gtf"></a>What is lr2gtf ?
lr2gtf is a [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based light-weight pipeline which is designed to utilize both third-generation long-read and second-generation short-read RNA-seq data to generate an enhanced gene annotation file.


## Table of Contents

- [What is lr2gtf ?](#lr2gtf)
- [Installation](#install)
  - [Operating system](#os)
  - [Cloning and building lr2gtf pipeline](#build)
  - [Dependencies](#depen)
- [Getting started with provided toy example in `test_data`](#start)
- [Input and output](#input_output)
  - [Input files](#input)
  - [Output file](#output)
- [Running lr2gtf on local machine](#local)
- [Running lr2gtf on computer cluster](#cluster)
- [More about `snakemake` and configuration file](#snakemake)
- [FAQ](#FAQ)
- [Contact](#contact)


## <a name="install"></a>Installation
### <a name="os"></a>Operating system
lr2gtf currently can only be complied and run on Linux/Unix systems.

### <a name="build"></a>Clone and build lr2gtf pipeline
```
git clone https://github.com/Xinglab/lr2gtf.git --recursive
cd lr2gtf
make dependencies
make
```
`make dependencies` command will build all dependencies that are needed by lr2gtf. `make` will build main program of lr2gtf pipeline.

After building is done, path of `lr2gtf/bin` needs to be add to the environment variable PATH.

### <a name="depen"></a>Dependencies 
lr2gtf is dependent on following open-source software: [minimap2](https://github.com/lh3/minimap2), [STAR](https://github.com/alexdobin/STAR), [samtools](https://github.com/samtools/samtools) and [Snakemake](https://snakemake.readthedocs.io/en/stable/).

They will be automatically downloaded and built via `make dependencies` command in above section.

If some of them are already installed in your machine, you can choose to only build one of them separately, like `make minimap2` or `make snakemake`.


## <a name="start"></a>Getting started with provided toy example in `test_data`
```
snakemake -p --configfile ./config.yaml updated.gtf
``` 
Enhanced gene annotation file `updated.gtf` will be generated in current working directory, along with some intermediate files.

## <a name="input_output"></a>Input and output
All the input and output files are specified in the configuration file `config.yaml`: 
```
# input
genome:
    fasta: test_data/genome/genome.fa
    gtf: test_data/gtf/original.gtf
    rRNA: test_data/gtf/rRNA.gtf

sample:
    long_read:
        samp1: test_data/read/samp1_long.fa
    short_read:
        samp1:
            first: test_data/read/samp1_short_1.fa
            second: test_data/read/samp1_short_2.fa
# output
output:
    gtf: updated.gtf 
```

### <a name="input"></a>Input files
`genome.fa` is the genome file in FASTA format.

`original.gtf` is the original existing gene annotation file, `rRNA.gtf` is the GTF file for all the ribosomal RNA.

`samp1_long.fa` is long-read data of sample #1.
 
`samp1_short_1.fa` and `samp1_short_2.fa` are paired-end short-read data of sample #1 

### <a name="output"></a>Output file
`updated.gtf` is the generated enhanced gene annotation file, which contains both known annotation(`original.gtf`) and reliable novel transcript information extracted from long- and short-read data.


## <a name="local"></a>Running lr2gtf on local machine
`snakemake -p --snakefile ./Snakefile --configfile ./config.yaml --jobs 8`

8 threads will be used in parallel on the local machine. For more details about how to run lr2gtf on local machine, refer to [More about `snakemake` and configuration file](#snakemake).

 
## <a name="cluster"></a>Running lr2gtf on computer cluster
```
snakemake -p --snakefile ./Snakefile --configfile ./config.yaml  --cluster-config ./config.yaml \
--cluster "qsub -cwd -V -l h_data={cluster.h_data},h_rt={cluster.h_rt} -pe shared {cluster.threads} \"
```
Computing jobs will be automatically submitted to the computer cluster. Allocation information are also specified in `config.yaml`. For more details about how to run lr2gtf on computer cluster, refer to [More about `snakemake` and configuration file](#snakemake).   

## <a name="snakemake"></a>More about `snakemake` and configuration file
1. All the input and output files that are in relative path format will use current working directory or directory specified by `snakemake` argument `--directory`(`-d`) as origin.
2. Maximum cores to be used on local machine or computer cluster could be specified with `snakemake` argument `--jobs` (`--cores`, `-j`).
3. Computing resources for each job could be set in configuration file. 
    * `h_data`: memory size
    * `h_rt`: running time
    * `threads`: multi-threads

## <a name="FAQ"></a>FAQ

1. **Q**: How many samples are needed for lr2gtf?

   **A**: No limit on sample size. As long as one sample has matched long & short-read data, it can be provided to lr2gtf.

2. **Q**: How to specify the directory of output and intermediate files?
   **A**: Use `snakemake` argument `--directory`(`-d`) to specify working directory. Note that when working directory is set, all the relative path in configuration file will use it as origin. Or, you could just use absolute path instead.
   
3. **Q**: How to specify single-end short-read data in configuration file?

   **A**: Write single-end data file path after `first:` and leave `second:` empty like this:
   ```
   first: test_data/read/short_single.fa
   second: []
   ```
   
4. **Q**: How to specify that long and short-read data match with each other?

   **A**: Use uniform name for data from same sample, like this:
   ```
   long_read:
      samp1: test_data/read/samp1_long.fa
      samp2: test_data/read/samp2_long.fa
   short_read:
       samp1:
           first: test_data/read/samp1_short_1.fa
           second: test_data/read/samp1_short_2.fa
       samp2:
           first: test_data/read/samp2_short_single.fa
           second: []
   ```
   Here `samp1` and `samp2` are the names of two samples, they should be consistent in `long_read` and `short_read`. 


## <a name="contact"></a>Contact
Yan Gao yangao@ucla.edu

Yi Xing yxing@ucla.edu

[github issues](https://github.com/Xinglab/lr2gtf/issues)

  
