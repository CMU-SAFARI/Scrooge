
# **Preparing Datasets**

## **Reference Genome**

Reference genome files (.fasta) can be obtained from various sources, e.g., the [human reference genome GRCh38](https://www.ncbi.nlm.nih.gov/genome/guide/human/).
To ensure only ACGT bases remain (see above), we provide the `GenConverter` script:

```bash
python3 GenConverter.py --restrict ACGT --prune_titles --if FASTA GCF_000001405.39_GRCh38.p13_genomic.fna filtered_genome.fasta
```

## Reads and Candidate Locations

Reads and candidate locations can be obtained in two ways:
1. Simulate reads from the reference genome, using a tool such as [PBSIM2](https://github.com/yukiteruono/pbsim2) (outputs .fastq and .maf)
2. Obtain real reads (.fastq) from a database such as the [NCBI Sequence Read Archive](https://trace.ncbi.nlm.nih.gov/Traces/sra/), and then obtain candidate locations from a mapper such as [minimap2](https://github.com/lh3/minimap2) (outputs .paf)

### **Simulating Reads Using PBSIM2**

To facilitate dataset generation with PBSIM2, we provide a docker container that comes with PBSIM2 installed. Alternatively, PBSIM2 can be installed from its [repository](https://github.com/yukiteruono/pbsim2).

```bash
#open an interactive bash shell inside the docker
bash ./start_pbsim2_docker.sh
pbsim --help
```

PBSIM2 produces a reads file (.fastq) and a file with the ground truth locations of the reads in the reference genome (.maf). Scrooge can be evaluated on those ground truth locations, or the reads can be mapped to the reference genome with minimap2 (results in a .paf).

### **Downloading Reads from NCBI**

Genomic read files (.fastq) can be obtained from various sources, e.g., the [SRR13278681](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13278681) read set from a human cell line.

Read read sets should be sanitized as follows:
```bash
python3 GenConverter.py --restrict ACGT --prune_titles SRR13278681_1.fastq SRR13278681_1_sanitized.fastq
```

### **Mapping Reads with minimap2**`

Scrooge requires candidate locations to evaluate its pairwise alignment performance. Candidate locations can be obtained from minimap2 as a [.paf file](https://lh3.github.io/minimap2/minimap2.html#10). Depending on how minimap2 is configured, these locations can be either mapped locations (default), or all candidate locations remaining after minimap2's chaining step (with the [-P flag](https://lh3.github.io/minimap2/minimap2.html#5)).


# Reproducing our Datasets

## PBSIM2 groundtruth
1. Sanitize the [GRCh38.p14](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) human reference genome as shown above
2. Run PBSIM2:
    ```
    pbsim --prefix pacbio-chr1-simulated-m10k-k5 --depth 10.0 --length-min 10000 --length-max 10000 --length-mean 10000.0 --length-sd 0 --accuracy-min 0.95 --accuracy-max 0.95 --accuracy-mean 0.95 --difference-ratio 6:50:54 --hmm_model /pbsim2/data/P6C4.model filtered_genome.fasta
    ```
3. Use the generated reads of the first chromosome and corresponding ground truth locations in the .maf file

## PBSIM2 chained
1. Prepare `PBSIM2 groundtruth` as above
2. Truncate the reads file to the first 5000 reads, e.g., using:
    ```bash
    head -5000 [...] > reads.fastq
    ```
3. Run minimap2 (with 48 threads):
    ```bash
    minimap2 -x map-pb reference.fasta reads.fastq candidates.paf -P -t 48
    ```

## Illumina mapped

1. Sanitize the [GRCh38.p14](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) human reference genome as shown above
2. Download and sanitize the [SRR13278681_1 dataset](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13278681)
3. Run minimap2 (with 48 threads):
    ```bash
    minimap2 -x sr reference.fasta reads.fastq -o candidates.paf -t 48
    ```

## Illumina chained
1. Prepare the Illumina files as above
2. Truncate the reads file to the first 100,000 reads, e.g., using:
    ```bash
    head -100000 [...] > reads.fastq
    ```
3. Run minimap2 (with 48 threads). Note the addition of `--secondary==yes` to override [minimap2's default](https://lh3.github.io/minimap2/minimap2.html#8) `--secondary==no` for the short read mapping preset:
    ```bash
    minimap2 -x sr reference.fasta reads.fastq -o candidates.paf -t 48  --secondary=yes -P
    ```

## HiFi mapped
1. Sanitize the [GRCh38.p14](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) human reference genome as shown above
2. Download and sanitize [SRR12519035.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRR12519035)
3. Truncate the reads file to the first 400,000 reads, e.g., using:
    ```bash
    head -400000 [...] > reads.fastq
    ```
4. Run minimap2 (with 48 threads):
    ```bash
    minimap2 -x map-hifi reference.fasta reads.fastq -o candidates.paf -t 48
    ```

## HiFi chained
1. Prepare the HiFi files as above
2. Truncate the reads file to the first 1,000 reads, e.g., using:
    ```bash
    head -1000 [...] > reads.fastq
    ```
3. Run minimap2 (with 48 threads):
    ```bash
    minimap2 -x map-hifi reference.fasta reads.fastq -o candidates.paf -t 48 -P
    ```

## ONT mapped
1. Sanitize the [GRCh38.p14](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) human reference genome as shown above
2. Download and sanitize [SRR12564436.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRR12564436)
3. Run minimap2 (with 48 threads):
    ```bash
    minimap2 -x map-ont reference.fasta reads.fastq -o candidates.paf -t 48

## ONT chained
1. Prepare the ONT files as above
2. Truncate the reads file to the first 100 reads, e.g., using:
    ```bash
    head -100 [...] > reads.fastq
    ```
3. Run minimap2 (with 48 threads):
    ```bash
    minimap2 -x map-ont reference.fasta reads.fastq -o candidates.paf -t 48 -P
