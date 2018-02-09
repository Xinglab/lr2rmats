# lr2rmats: Long read to rMATS
## <a name="lr2rmats"></a>What is lr2rmats ?
lr2rmats is a [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based light-weight pipeline which is designed to utilize both third-generation long-read and second-generation short-read RNA-seq data to generate an enhanced gene annotation file. The newly generated annotation file could be provided to [rMATS](http://rnaseq-mats.sourceforge.net/) for differential alternative splicing analysis.


## Table of Contents

- [What is lr2rmats ?](#lr2rmats)
- [Installation](#install)
  - [Operating system](#os)
  - [Cloning and building lr2rmats pipeline](#build)
  - [Dependencies](#depen)
- [Getting started with toy example in `test_data`](#start)
- [Input and output](#input_output)
  - [Input files](#input)
  - [Output files](#output)
  - [Intermediate and log files](#intermediate)
- [Running lr2rmats on a local machine](#local)
- [Running lr2rmats on a computer cluster](#cluster)
- [More about `snakemake` and configuration file](#snakemake)
- [FAQ](#FAQ)
- [Contact](#contact)


## <a name="install"></a>Installation
### <a name="os"></a>Operating system
lr2rmats currently can only be built and run on Linux/Unix systems.

### <a name="build"></a>Cloning and building lr2rmats pipeline
```
git clone https://github.com/Xinglab/lr2rmats.git --recursive
cd lr2rmats
make dependencies
make lr2rmats
```
`make dependencies` command will build all dependencies that are needed by lr2rmats. `make lr2rmats` will build the main program of the lr2rmats pipeline. 

Alternatively, you could simply type `make` to build everything you need.

After the building is done, the path of `lr2rmats/bin` needs to be added to the environment variable PATH.

### <a name="depen"></a>Dependencies 
lr2rmats is dependent on following open-source software: [minimap2](https://github.com/lh3/minimap2), [STAR](https://github.com/alexdobin/STAR), [samtools](https://github.com/samtools/samtools) and [Snakemake](https://snakemake.readthedocs.io/en/stable/).

They will be automatically downloaded and built via `make dependencies` command if they are not installed in your computer.

You can choose to build them separately, like `make minimap2`, `make snakemake`.


## <a name="start"></a>Getting started with toy example in `test_data`
```
snakemake -p --snakefile ./Snakefile --configfile ./config.yaml
``` 
Enhanced gene annotation file `output/updated.gtf` will be generated in current working directory, along with some intermediate and log files.

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
    updated_gtf: output/updated.gtf
    
```

### <a name="input"></a>Input files
- `genome.fa`: genome file in FASTA format (***required***).

- `original.gtf`: original existing gene annotation file (***required***).

- `rRNA.gtf`: GTF file for all the ribosomal RNA. Long-read alignment that overlaps with rRNA will be skipped (***optional***).

- `samp1_long.fa`: long-read data of sample #1 (***required***).
 
- `samp1_short_1.fa` and `samp1_short_2.fa`: paired-end short-read data of sample #1 (***required***).

For single-end read multiple samples data, please refer to [FAQ](#FAQ). 

### <a name="output"></a>Output files
- `updated.gtf`: enhanced gene annotation file. It contains both known annotation(`original.gtf`) and reliable novel transcript information extracted from long and short-read data.

In addition to the updated annotation file, more detailed information are also provided:

- `known.gtf`: all known transcripts from long-read that already exist in the input annotation file.

- `novel.gtf`: all novel transcripts from long-read that are different from any of the annotation transcript. Here, only novel transcripts that have at least one known splice site are kept. 

- `unrecog.gtf`: all unrecognized transcripts from long-read that have no known splice site.

- `detail.txt`: detailed information for each long-read derived transcript.

- `summary.txt`: summary statistics of the whole pipeline. It includes the summary of annotation and updated/known/novel/unrecognized transcript information. 

""

### <a name="intermediate"></a>Intermediate and log files
Intermediate files will be generated in four folders: `alignment`, `gtf`, `logs`, and `benchmark`.
* `alignment` contains all the alignment result files of long-read and short-read data.
* `gtf` contains all the intermediate GTF files.
* `logs` contains the log file of each job.
* `benchmark` contains each job's running time, memory usage and other measurements.

## <a name="local"></a>Running lr2rmats on a local machine
`snakemake -p --snakefile ./Snakefile --configfile ./config.yaml --cores 8`

8 threads will be used in parallel on the local machine. For more details about how to run lr2rmats on a local machine, refer to [More about `snakemake` and configuration file](#snakemake).

 
## <a name="cluster"></a>Running lr2rmats on a computer cluster
```
snakemake -p --snakefile ./Snakefile --configfile ./config.yaml  --cluster-config ./config.yaml \
--cluster "qsub -cwd -V -l h_data={cluster.h_data},h_rt={cluster.h_rt} -pe shared {cluster.threads}"
```
Computing jobs will be automatically submitted to the computer cluster. Allocation information are specified in the configuration file `config.yaml`. For more details about how to run lr2rmats on a computer cluster, refer to [More about `snakemake` and configuration file](#snakemake).   

## <a name="snakemake"></a>More about `snakemake` and configuration file
1. All the input and output files that are in relative path format will use current working directory or directory specified by `snakemake` argument `--directory`(`-d`) as the origin.
2. Maximum number of cores to be used on a local machine and number of cluster nodes to be used on a computer cluster could be specified with `snakemake` argument `--jobs`(`--cores`, `-j`).
3. Computing resources for each job could be set in the configuration file. 
    * `h_data`: memory size
    * `h_rt`: running time
    * `threads`: number of threads

## <a name="FAQ"></a>FAQ

1. **Q**: How many samples are needed for lr2rmats?

   **A**: There is no limit on sample amount. As long as one sample has matched long and short-read data, it can be provided to lr2rmats.

2. **Q**: Do I need to provide read data in FASTA format?

   **A**: lr2rmats works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats. 
   
3. **Q**: How to specify the directory of output and intermediate files?

   **A**: You can use `snakemake` argument `--directory`(`-d`) to specify working directory. Note that when working directory is set, all the relative paths in the configuration file will use it as the origin. Or, you could just use the absolute path instead.
   
4. **Q**: How to use single-end short-read data as input?

   **A**: You can write single-end data file path after `first:`, and leave `second:` empty like this:
   ```
   first: test_data/read/short_single.fa
   second: []
   ```
   
5. **Q**: How to specify matched pairs of long and short-read data?

   **A**: You should use uniform name for data from same sample, like this:
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

[github issues](https://github.com/Xinglab/lr2rmats/issues)

  
