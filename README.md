# **Scrooge:** A fast and memory-frugal genomic sequence aligner for CPUs, GPUs and ASICs

_Scrooge_ is a fast pairwise genomic sequence aligner. It efficiently aligns short and long genomic sequence pairs on multiple computing platforms. It is based on the GenASM algorithm ([Senol Cali+, 2020](https://arxiv.org/abs/2009.07692)), and adds multiple algorithmic improvements that significantly improve the throughput and resource efficiency for CPUs, GPUs and ASICs. For long reads, the CPU version of Scrooge achieves a 20.1x, 1.7x, and 2.1x speedup over KSW2, Edlib, and a CPU implementation of GenASM, respectively. The GPU version of Scrooge achieves a 4.0x 80.4x, 6.8x, 12.6x and 5.9x speedup over the CPU version of Scrooge, KSW2, Edlib, Darwin-GPU, and a GPU implementation of GenASM, respectively. We estimate an ASIC implementation of Scrooge to use 3.6x less chip area and 2.1x less power than a GenASM ASIC while maintaining the same throughput.

This repository contains Scrooge's CPU and GPU implementations, and several evaluation scripts. We describe Scrooge in our paper [on arXiv](https://doi.org/10.48550/arXiv.2208.09985) and [in Bioinformatics](https://doi.org/10.1093/bioinformatics/btad151).

## **Citing Scrooge**

If you use Scrooge in your work, please cite:

> Joël Lindegger, Damla Senol Cali, Mohammed Alser, Juan Gómez-Luna, Nika Mansouri Ghiasi, and Onur Mutlu.
> ["Scrooge: A Fast and Memory-Frugal Genomic Sequence Aligner for CPUs, GPUs and ASICs."](https://doi.org/10.1093/bioinformatics/btad151)
> Bioinformatics (2023).

> Joël Lindegger, Damla Senol Cali, Mohammed Alser, Juan Gómez-Luna, and Onur Mutlu. 
> ["Algorithmic Improvement and GPU Acceleration of the GenASM Algorithm."](https://arxiv.org/abs/2203.15561) 
> HiCOMB (2022).

Below are the citations in bibtex format.

```bibtex
@article{lindegger2023scrooge,
  title={{Scrooge: A Fast and Memory-Frugal Genomic Sequence Aligner for CPUs, GPUs, and ASICs}},
  author={Lindegger, Jo{\"e}l and Senol Cali, Damla and Alser, Mohammed and G{\'o}mez-Luna, Juan and Mansouri Ghiasi, Nika and Mutlu, Onur},
  journal={Bioinformatics},
  year={2023}
}
@article{lindegger2022algorithmic,
  title={{Algorithmic Improvement and GPU Acceleration of the GenASM Algorithm}},
  author={Lindegger, Jo{\"e}l and Senol Cali, Damla and Alser, Mohammed and G{\'o}mez-Luna, Juan and Mutlu, Onur},
  journal={HiCOMB},
  year={2022}
}
```

## **Repository Structure**

```
.
└── 1. baseline_algorithms
└── 2. scripts
└── 3. src
└── 4. profile #after running scripts/profile.py or scripts/download_profile.py
└── 5. cacti #after running scripts/asic_numbers.py --cacti
└── 6. datasets #after running scripts/download_datasets.py
```

1. the `baseline_algorithms` directory contains the source code of various CPU and GPU baseline tools to compare against
2. the `scripts` directory contains various python and bash tools to automate experiments, plotting, and data download/generation
3. The `src` directory contains the CPU and GPU code for Scrooge's alignment, testing and example components
4. The `profile` directory contains the `.csv` results of parameter sweeps done by `scripts/profile.py`. Alternatively, our profiling results can be downloaded from Zenodo with `scripts/download_profile.py`
5. The `cacti` directory contains a clone of the CACTI git repository after running `scripts/asic_numbers.py --cacti`
6. The datasets directory contains several example datasets after running `scripts/download_datasets.py`

## **Requirements**

- **Python 3:** At least version 3.6, tested with 3.9 and 3.10
    - packages: matplotlib, pandas
- **g++:** At least support for C++17 and OpenMP, tested with g++ 9.4
- **CUDA:** To run GPU tools, tested with CUDA 10 and CUDA 11
- **Docker:** To run the PBSIM2 docker we provide, if desired

## **Using Scrooge as a Library**

Scrooge can be used for pairwise sequence alignment in a variety of use-cases by calling it as a library. `library_example.cu` gives an example for each supported library interface, it can be compiled and run with `make library_example_linux`. All components of Scrooge run on Linux. The CUDA and C++ components are supported on Windows and have their own make rules (e.g., `make library_example_windows`).

```bash
git clone https://github.com/CMU-SAFARI/Scrooge && cd Scrooge
make library_example_linux
./library_example
```

There are two types of interfaces supported, *unstructured pairwise alignment* and *pairwise alignment for read mapping*, each on both CPUs and GPUs. Both calculate the semiglobal edit distance and corresponding alignment between pairs of strings, where the entire query/read sequence must be consumed, but only a prefix of the target/reference sequence.

### **Unstructured Pairwise Alignment**

This is the more general of the two supported interfaces. It simply accepts a list of pairs of strings, and calculates the edit distance and CIGAR string between each pair.

### **Pairwise Alignment for Read Mapping**

This is the more specialized, but potentially more efficient of the two supported interfaces. It accepts a reference genome of one or multiple chromosomes, a list of reads, and a list of candidate locations for each read. It then performs pairwise alignment between each read and the reference at each of the reads' candidate locations.

The key advantage of this interface is that it does not create any redundant copies of the sequence pairs in memory: Each read is stored once, and a single copy of the reference genome is stored. This can improve the memory footprint, required memory bandwidth, and cache hitrate. In contrast, when read mapping with the more general unstructured interface, each candidate location has a separate copy of the read and reference segment in memory.

## **Performance and Accuracy Evaluation**

We provide the `tests` command line utility to evaluate Scrooge's throughput and accuracy in a read mapping use-case, it can be compiled as follows:

```bash
git clone https://github.com/CMU-SAFARI/Scrooge && cd Scrooge
make tests_linux
./tests --unit_tests
```

The `tests` utility accepts any short or long read dataset as a tuple of &lt;FASTA, FASTQ, MAF/PAF&gt; files, i.e., a reference genome, a read set, and a candidate location list. This simulates the case that seeding and chaining were already executed, and the remaining pairs must be evaluated using pairwise sequence alignment. The only requirement on the dataset is that the fasta and fastq files must contain only upper- and lowercase "ACGT" characters, but no "N" (nondetermined) bases. For example, it can be invoked on one of our prepared datasets:

```bash
python3 scripts/download_datasets.py
./tests --reference=./datasets/pbsim_groundtruth/reference.fasta --reads=./datasets/pbsim_groundtruth/reads.fastq  --seeds=./datasets/pbsim_groundtruth/candidates.maf
```
This runs a GPU performance test by default, and produces an output in the following form:
```
1 visible GPU(s):
idx=0 name="NVIDIA GeForce RTX 3060" SMs=28 smem=100kiB

align_all() took 96202ms (data transfers, conversion, gpu kernel and post-processing)
GPU kernel tool 4608ms
GPU kernel ran at 25004 aligns/second
```
Alternatively, the `tests` utility can also run a CPU performance test, producing a similar output:
```bash
./tests --reference=./datasets/pbsim_groundtruth/reference.fasta --reads=./datasets/pbsim_groundtruth/reads.fastq  --seeds=./datasets/pbsim_groundtruth/candidates.maf --cpu_performance_test
```

We provide the `cpu_baseline` utility to evaluate the performance of CPU baselines, and compare Scrooge's accuracy to the accuracy of baseline tools.

```bash
make cpu_baseline_linux
./cpu_baseline --help
```

## **Advanced Compilation**

When running Scrooge on a GPU, the correct CUDA GPU architecture should be supplied to nvcc. The makefile we provide checks the architecture set in the `NVCC_ARCH` environment variable, e.g.:
```
export NVCC_ARCH=sm_86
make tests_linux
```
If `NVCC_ARCH` is not set, the makefile will default to sm_86. If `nvcc` is not in the `PATH`, the nvcc binary can be specified using the `NVCC` environment variable.

## **Automated Profiling**

Alternatively, the `profile` script automatically compiles and runs `tests` across multiple passes and algorithm configurations for a dataset. It assumes each dataset is in a separate subdirectory of `datasets` and contains a reference, read and candidate location file. The `download_datasets` script produces this file structure. For example, it can be run as follows:
```bash
python3 scripts/download_datasets.py
python3 scripts/profile.py cpu pbsim_groundtruth
```

The `profile` script then produces CSV files with different columns, depending on the profiling target. For example, `python3 scripts/profile.py cpu [...]` produces the following:
```csv
W,O,SENE,DENT,early termination,threads,aligns/second
64,2,False,False,False,48,11607.4
...
```

For all options, see 
```bash
python3 scripts/profile.py --help
```
## **Obtaining Datasets**

All our prepared datasets can also be downloaded from [Zenodo](https://zenodo.org/record/7013734/files/scrooge_datasets.tar.gz). 

The `download_datasets` script automatically downloads and untars our prepared datasets:
```bash
python3 scripts/download_datasets.py
```
The gzipped datasets are 13.9GiB in total, and 63GiB when unzipped.

In [DATASETS.md](DATASETS.md) we describe in detail how we prepared our datasets, and how other real and simulated datasets can be prepared.

 For reproducibility, we also provide the outputs we obtained from `profile.py` on [Zenodo](https://zenodo.org/record/6736836/files/scrooge_profile_results.tar.gz). The `download_profile` script automatically downloads and extracts the files:
 ```
 python3 scripts/download_profile.py
 ```
 The gzipped `profile` directory is 3.1GiB in total, and 17GiB when unzipped.

## **Plotting**

The `plot` script produces several figures that analyze Scrooge's throughput and accuracy. It assumes all files in our profile results are available, i.e. should be used in conjunction with `download_profile`.

```bash
python3 scripts/download_profile.py
python3 scripts/plot.py
```

## **CIGAR inspector**

We provide a specialized plotting script called `cigar_inspector` to visualize the CIGAR strings different algorithms produce for each sequence pair. This is particularly helpful to understand under what circumstances GenASM produces different results than e.g. KSW2 and Edlib.

This should we run on a csv file generated by the `profile` script that includes the CIGAR strings:
```bash
python3 scripts/profile.py accuracy_cpu --cigar [dataset]
python3 scripts/cigar_inspector.py profile/[...]_all_accuracy_cigar.csv 10
```
The above command starts the CIGAR inspector with the sequence pair for which Scrooge/GenASM generated the 10th-worst alignment score. You can zoom into regions of interest using the matplotlib GUI functions. At higher zoomlevels, the axes ticks are labeled with the sequence contents.

## **Getting Help**

If you have any suggestion for improvement, please contact jmlindegger at gmail dot com. If you encounter bugs or have further questions or requests, you can raise an issue at the [issue page](https://github.com/CMU-SAFARI/Scrooge/issues).

## **Limitations**

Like GenASM, Scrooge may calculate an edit distance and alignment that is very close to the actual edit distance, but sometimes overestimates the actual edit distance, and finds a sub-optimal alignment. We evaluate GenASM/Scrooge's accuracy behavior in detail in our paper.

## **Reproducing Results**

All results from our paper can be reproduced by repeating our exact experiments with the following commandlines. Note that (1) some of these experiments run 10s of hours, and (2) the GPU command lines require an NVIDIA GPU to be available.

```bash
#download datasets
python3 scripts/download_datasets.py

#calculate ASIC area and power
python3 ./scripts/asic_numbers.py --cacti

#accuracy
python3 ./scripts/profile.py accuracy_cpu illumina_chained --override_W 32
python3 ./scripts/profile.py accuracy_cpu pbsim_chained
python3 ./scripts/profile.py accuracy_cpu pbsim_groundtruth
python3 ./scripts/profile.py accuracy_cpu pbsim_groundtruth --cigar

#CPU performance
python3 ./scripts/profile.py cpu illumina_chained --override_W 32
python3 ./scripts/profile.py cpu pbsim_chained
python3 ./scripts/profile.py cpu_baselines illumina_chained
python3 ./scripts/profile.py cpu_baselines pbsim_chained

#GPU performance
python3 ./scripts/profile.py gpu illumina_chained --arch sm_86 --override_W 32
python3 ./scripts/profile.py gpu pbsim_chained --arch sm_86
python3 ./scripts/profile.py gpu_baselines illumina_chained --tmp_dir tmp --arch sm_86
python3 ./scripts/profile.py gpu_baselines pbsim_chained --tmp_dir tmp --arch sm_86

#plot
python3 ./scripts/plot.py
```
