# DeepSAP

**DeepSAP** is a transformer-based workflow designed to enhance splice junction detection in RNA-seq data. By default, DeepSAP utilizes a highly sensitive GPU-accelerated GSNAP TGGA aligner for FASTQ inputs. Alternatively, it can also score pre-aligned BAM files directly — either from GSNAP itself or from any other aligner whose SAM records carry the `XA` (alternative alignments) tag. <br>

<!-- INTERNAL_README_START
We evaluated the performance of DeepSAP in our *Genome Biology* article: [***DeepSAP: improved RNA-seq alignment by integrating transcriptome guidance with transformer-based splice junction scoring***](https://doi.org/10.1186/s13059-026-04100-3) (Berakdar, Wu, Zhu, Samadi, Vats, 2026). In our benchmark, DeepSAP demonstrated strong performance, achieving consistently outstanding results across all evaluated metrics using Baruzzo et al. datasets. For additional resources, including data, detailed analyses, and other supplementary materials accompanying the **DeepSAP** article, please refer to the [DeepSAP GitHub repository](https://github.com/clara-parabricks-workflows/DeepSAP/tree/main).
INTERNAL_README_END -->

<!-- PUBLIC_GITHUB_README_START -->
We evaluated the performance of DeepSAP in our *Genome Biology* article: [***DeepSAP: improved RNA-seq alignment by integrating transcriptome guidance with transformer-based splice junction scoring***](https://doi.org/10.1186/s13059-026-04100-3) (Berakdar, Wu, Zhu, Samadi, Vats, 2026). In our benchmark, DeepSAP demonstrated strong performance, achieving consistently outstanding results across all evaluated metrics using Baruzzo et al. datasets.
<br>

<img src="manuscript_data_code/figures_plots_and_data/Figure_3/Figure_3_a_Baruzzo_spirder_recalls_vs_precision.png" width="1000">

For additional resources, including data, detailed analyses, and supplementary materials accompanying the **DeepSAP** article, please refer to [`manuscript_data_code/README.md`](manuscript_data_code/README.md) in this repository.<!-- PUBLIC_GITHUB_README_END -->

For questions, bug reports, or other DeepSAP support requests, please use the [Parabricks developer forum](https://forums.developer.nvidia.com/c/healthcare/parabricks/).
<br>

## Table of Contents
- [Requirements](#requirements)
- [Usage](#usage)
  - [Least-privilege container execution](#least-privilege-container-execution)
- [Pipeline Modes](#pipeline-modes)
- [Command-line Arguments](#command-line-arguments)
- [Version History](#version-history)
- [License/Terms of Use](#licenseterms-of-use)

## Requirements
### System Software:
* Docker with GPU support<br>
### System Hardware:
Sizing below is for a human genome–scale reference (GRCh38). The two pipeline stages run sequentially, so peak GPU memory is the **maximum** of the alignment-stage and TSJS-stage footprints — not their sum.

**CPU & RAM:**
* **CPU:** 24 cores recommended (drives GSNAP's pipeline-parallel stages — reader / solver / writer threads — and DeepSAP's TSJS scoring stage).
* **System RAM:** 64 GB minimum.

**GPU memory:**
* **Minimum recommended:** 40 GB (validated on NVIDIA A100 PCIe 40 GB, H100 PCIe, and RTX A6000 48 GB).
* The alignment stage sets the floor; the TSJS stage's GPU memory scales with `--batch` and `--fp16`.

*Alignment stage (GPU-accelerated GSNAP):*
* GSNAP transcriptome-guided genome index resident on device: **~24 GB**.
* `--localdb-scratch` (Stage-2 localdb GPU scratch buffer): default **`12G`**, tunable.
* **Default total: ~36 GB.** Setting `--localdb-scratch=1G` brings the alignment-stage footprint down to **~25 GB** (fits a 24 GB card with little headroom — closer to 32 GB is comfortable).

*TSJS (transformer splice-junction scoring) stage:*

GPU memory here is dominated by two parameters:
* **`--batch`**: number of candidate splice junctions scored per transformer forward pass. Larger batches significantly improve throughput but require more GPU memory.
* **`--fp16`**: half-precision floating-point inference is enabled by default and roughly halves GPU memory versus fp32. Disable with `--no-fp16`, which approximately doubles the per-batch memory shown below.

  | **`--batch`** | **Approximate GPU memory (with `--fp16`)** |
  | :------------ | :----------------------------------- |
  | 64            | ~1.2 GB                              |
  | 128           | ~1.6 GB                              |
  | 256           | ~2.2 GB                              |
  | 2048 (default)| ~10.4 GB                             |
  | 8192          | ~39.5 GB                             |

### Input Data:
* RNA-seq reads in FASTQ format.
* Reference file in FASTA format.
* Annotation file in GTF format.
* Optionally, a path to a GSNAP index.

## Usage
This guide demonstrates how to quickly test DeepSAP's functionality using the **`malaria_short_pe`** dataset. Follow these steps to set up your environment and run DeepSAP:

### Least-privilege container execution

For production runs, mount only the input, index, and output directories that DeepSAP needs. Keep reference data, FASTQ/BAM inputs, and reusable indexes read-only; use a separate writable output directory for results.

### Step 1: Prepare Environment and Download Test Data
This step downloads the latest DeepSAP Docker container and all required reference files and test sequencing data.

```bash
# Pull the DeepSAP Parabricks Docker image
docker pull nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest

# Create writable result directories as the invoking host user
mkdir -p test/outputdir test/results/prebuilt test/results/rescored

# Download reference genome and annotation files
wget -P test/malaria_short_pe/ https://raw.githubusercontent.com/clara-parabricks-workflows/DeepSAP/main/test/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf
wget -P test/malaria_short_pe/ https://raw.githubusercontent.com/clara-parabricks-workflows/DeepSAP/main/test/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa

# Download downsampled FASTQ sequence reads (10K) from DeepSAP GitHub
wget -P test/malaria_short_pe/ https://raw.githubusercontent.com/clara-parabricks-workflows/DeepSAP/main/test/malaria_short_pe/SRR14793977_10K_1.fastq.gz
wget -P test/malaria_short_pe/ https://raw.githubusercontent.com/clara-parabricks-workflows/DeepSAP/main/test/malaria_short_pe/SRR14793977_10K_2.fastq.gz
```

### Step 2: Build a GSNAP Index (`--mode index`)
This command builds a standalone, reusable GSNAP TGGA index from the FASTA + GTF and writes it under `<out>/<prefix>/`. Useful when you plan to score many samples against the same reference — build the index once, then reuse it in Step 4.

```bash
# Build a reusable GSNAP index from the malaria reference
docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir:/outputdir:rw"                                  \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode index                                                                    \
    --out /outputdir/                                                               \
    --prefix malaria_idx                                                            \
    --gtf /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf           \
    --fasta /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa
```

### Step 3: Run DeepSAP End-to-End (auto-build index, `--mode GSNAP+TSJS`)
This command executes the full DeepSAP pipeline (GSNAP alignment + transformer splice-junction scoring) on the downloaded test dataset using the default `--mode GSNAP+TSJS`. Since `--gsnap_idx` is not specified, DeepSAP auto-builds a GSNAP index inline at `<out>/gsnap_idx/` before alignment. Pick this path for one-shot runs where you don't need to reuse the index later.

```bash
# Run DeepSAP end-to-end (GSNAP index will be auto-generated)
docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir:/outputdir:rw"                                  \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_10K                                                           \
    --mate_1 /workdir/malaria_short_pe/SRR14793977_10K_1.fastq.gz                   \
    --mate_2 /workdir/malaria_short_pe/SRR14793977_10K_2.fastq.gz                   \
    --gtf /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf           \
    --fasta /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa
```

### Step 4: Run DeepSAP with a Pre-existing GSNAP Index (`--mode GSNAP+TSJS` + `--gsnap_idx`)
If you have already generated a GSNAP index (e.g., from Step 2, a previous DeepSAP run, or shared infrastructure), point DeepSAP at it via `--gsnap_idx`. This skips index construction, runs the same GSNAP alignment pipeline as an auto-indexed run, writes `<prefix>_gsnap.bam`, and then applies TSJS scoring to that BAM.

```bash
# Run DeepSAP using the index built in Step 2
docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir/malaria_idx:/gsnap_idx:ro"                      \
    --volume "$(pwd)/test/results/prebuilt:/outputdir:rw"                           \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_10K                                                           \
    --mate_1 /workdir/malaria_short_pe/SRR14793977_10K_1.fastq.gz                   \
    --mate_2 /workdir/malaria_short_pe/SRR14793977_10K_2.fastq.gz                   \
    --gtf /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf           \
    --fasta /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa\
    --gsnap_idx /gsnap_idx/
```

### Step 5: Score an Existing BAM with TSJS Only (`--mode GSNAP+TSJS` + `--sam`)
If you already have a GSNAP-aligned BAM (e.g., from a prior GSNAP alignment run, or from any other aligner whose SAM records carry the `XA` (alternative alignments) tag), pass it via `--sam` and DeepSAP skips alignment entirely — running transformer splice-junction scoring directly on the BAM. The output is a new BAM with TSJS-derived MAPQ adjustments and junction-scoring metadata.

```bash
# Score a pre-aligned BAM (no GSNAP step)
docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir/test_run_10K_gsnap.bam:/input.bam:ro"           \
    --volume "$(pwd)/test/results/rescored:/outputdir:rw"                           \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_10K_rescored                                                  \
    --sam /input.bam                                                                \
    --gtf /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf           \
    --fasta /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa
```

> Note: `--sam` and `--mate_1`/`--mate_2` are mutually exclusive — DeepSAP either aligns *or* scores an existing alignment, never both in the same run.

<!-- ```bash
# Run DeepSAP with the test dataset and pre-generated GSNAP index
docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir/malaria_idx:/gsnap_idx:ro"                      \
    --volume "$(pwd)/test/results/prebuilt:/outputdir:rw"                           \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_10K                                                           \
    --mate_1 /workdir/malaria_short_pe/SRR14793977_10K_1.fastq.gz                   \
    --mate_2 /workdir/malaria_short_pe/SRR14793977_10K_2.fastq.gz                   \
    --gtf /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf           \
    --fasta /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa\
    --gsnap_idx /gsnap_idx/

docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir:/outputdir:rw"                                  \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_human_ERG                                                     \
    --gtf /workdir/human_ERG_short_pe/chromosome_21.gtf                             \
    --fasta /workdir/human_ERG_short_pe/chromosome_21.fa                            \
    --sam /workdir/human_ERG_short_pe/ERG_novel_sj.sam

docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir:/outputdir:rw"                                  \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_danio                                                         \
    --gtf /workdir/zebra_fish_short_pe/GCF_000002035.6_GRCz11_genomic.gtf           \
    --fasta /workdir/zebra_fish_short_pe/GCF_000002035.6_GRCz11_genomic.fna         \
    --mate_1 /workdir/zebra_fish_short_pe/SRR14793977_10K_1.fastq.gz                \
    --mate_2 /workdir/zebra_fish_short_pe/SRR14793977_10K_2.fastq.gz


docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir:/outputdir:rw"                                  \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_human_genecode                                                \
    --gtf /workdir/human_genecode_short_pe/gencode.v48.annotation.gtf               \
    --fasta /workdir/human_genecode_short_pe/GRCh38.p14.genome.fa                   \
    --mate_1 /workdir/human_genecode_short_pe/SRR14793977_10K_1.fastq.gz            \
    --mate_2 /workdir/human_genecode_short_pe/SRR14793977_10K_2.fastq.gz
``` -->


## Pipeline Modes

DeepSAP's `--mode` flag selects which pipeline mode to run. The default `GSNAP+TSJS` reproduces the v0.0.x end-to-end behavior; `index` lets you pre-build a GSNAP index in isolation (useful for sharing a pre-built index across many samples).

| `--mode` | Required inputs | Optional inputs | Outputs |
|---|---|---|---|
| `index` | `--fasta`, `--gtf` | — | GSNAP index at `<out>/<prefix>/` |
| `GSNAP+TSJS` (default) | `--fasta`, `--gtf`, **and** either `--mate_1`+`--mate_2` (optionally with `--gsnap_idx`) **or** `--sam` | model / batching flags, `--score_method` | scored BAM at `<out>/<prefix>.bam` (+ intermediate datasets) |

### Mode 1: Build a GSNAP index only

```bash
docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir:/outputdir:rw"                                  \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode index                                                                    \
    --out /outputdir/                                                               \
    --prefix malaria_idx                                                            \
    --gtf /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf           \
    --fasta /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa
```

### Mode 2: Score an existing BAM with TSJS only

```bash
docker run --gpus 1 --ulimit memlock=-1 --ulimit stack=67108864 --rm                \
    --user "$(id -u):$(id -g)"                                                      \
    --cap-drop=ALL --security-opt=no-new-privileges                                 \
    --volume "$(pwd)/test:/workdir:ro"                                              \
    --volume "$(pwd)/test/outputdir/test_run_10K_gsnap.bam:/input.bam:ro"           \
    --volume "$(pwd)/test/results/rescored:/outputdir:rw"                           \
    nvcr.io/nvidia/clara/clara-parabricks-deepsap:latest                            \
    --mode GSNAP+TSJS                                                               \
    --out /outputdir/                                                               \
    --prefix test_run_10K_rescored                                                  \
    --sam /input.bam                                                                \
    --gtf /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf           \
    --fasta /workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa
# -> /outputdir/test_run_10K_rescored.bam (TSJS-scored)
```

### DeepSAP Expected Output
```
[2025-07-18 12:51:27]   [INFO]  Running DeepSAP v0.1.0
[2025-07-18 12:51:32]   [LOG]   Running GSNAP
[2025-07-18 12:51:32]   [LOG]   Building GSNAP TGGA index
[2025-07-18 12:52:44]   [LOG]   Running GSNAP TGGA 
[2025-07-18 12:52:46]   [LOG]   Parsing FASTA file '/workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa'
[2025-07-18 12:52:46]   [LOG]   Parsing GTF file '/workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.60.gtf'
[2025-07-18 12:52:47]   [LOG]   Transcript information: 
Number of transcripts:             5767
Shortest transcript:               67   EPT00050203058
Longest transcript:                30863        CAG25094
Transcripts length mean:           2456.79
Transcripts length median:         1618
Transcripts length mode:           71
Shortest intron:                   1    PF3D7_1478200: 14__-__3219919__3220323 -> 14__-__3220325__3220534
Longest intron:                    2425 CZU00099: 14__+__1639681__1639728 -> 14__+__1642154__1642455
Introns length mean:               163.03
Introns length median:             141.0
Introns length mode:               1
Number of multi exons transcripts: 3064 53.13%
Number of mono exon transcripts:   2703 46.87%

Type of transcripts:
              BioType  Count  Percentage
0      protein_coding   5358       92.91
1          pseudogene    153        2.65
3               ncRNA    102        1.77
4                tRNA     79        1.37
5                rRNA     44        0.76
7                sRNA     17        0.29
6               snRNA     10        0.17
2  nontranslating_CDS      4        0.07
[2025-07-18 12:52:47]   [LOG]   Collecting splice junctions from GTF
[2025-07-18 12:52:47]   [LOG]   Collecting splice junctions in mode=NotStrict and window=150
[2025-07-18 12:52:47]   [LOG]   Collecting splice junctions from transcript types: All
Number of duplicated junctions:        328
Number of short junctions (intron):    0
Number of short junctions (donor):     0
Number of short junctions (acceptor):  0
Number of junctions contains N:        0
Number of accepted junctions:          8764
The First 10 Splicing Signals Types: 
Signal  Forward  Reverse  Percentage
  GTAG     4096     4431       97.30
  AAAA       18       17        0.40
  TATA       12        8        0.23
  GCAG        9        9        0.21
  TTTT        6        9        0.17
  ATAT        4        7        0.13
  GAGA        5        6        0.13
  AGAG        3        6        0.10
  TATT        3        6        0.10
  TAAT        4        5        0.10
[2025-07-18 12:52:47]   [LOG]   Collecting splice junctions from SAM/BAM file '/outputdir/test_run_10K_gsnap.bam'
[2025-07-18 12:52:47]   [INFO]  Sense junctions 518
[2025-07-18 12:52:47]   [INFO]  Antisense junctions 551
[2025-07-18 12:52:47]   [INFO]  Total number of reads 20479
[2025-07-18 12:52:47]   [INFO]  Total number of spliced reads 2233 10.903852727183946%
[2025-07-18 12:52:47]   [LOG]   Finished parsing a SAM file, len(found_junctions_table)= 1069
[2025-07-18 12:52:47]   [LOG]   Generating splice-junction prediction dataset batch: 1
[2025-07-18 12:52:47]   [LOG]   Writting dev.csv file for predicting into '/outputdir/test_run_10K_prediction_batch_1/'
[2025-07-18 12:52:47]   [LOG]   dev.csv file contains:   0: 1069, 1: 1069
[2025-07-18 12:52:47]   [LOG]   Predicting found splice junctions using DNABERT MS150
100%|██████████| 67/67 [00:01<00:00, 58.23it/s]
[2025-07-18 12:52:51]   [LOG]   Generating genome regions 
[2025-07-18 12:52:51]   [LOG]   Parsing FASTA file '/workdir/malaria_short_pe/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa'
[2025-07-18 12:52:53]   [LOG]   Finished writing BAM successfully into '/outputdir/test_run_10K'
[2025-07-18 12:52:53]   [LOG]   Number of SAM records: 20479 
[2025-07-18 12:52:53]   [LOG]   Number of reads IDs:   12644 
[2025-07-18 12:52:53]   [LOG]   Number of processed reads IDs: 1405  11.11% 

[2025-07-18 12:52:54]   [LOG]   Finished successfully
```
<br>

## Command-line Arguments

| Argument          | Description                                                                                          | Required       | Default |
|-------------------|------------------------------------------------------------------------------------------------------|----------------|---------|
| `--mode`          | Pipeline mode to run: `index` or `GSNAP+TSJS`. See [Pipeline Modes](#pipeline-modes).                | No             | `GSNAP+TSJS` |
| `-o, --out`       | Path to the output folder                                                                            | Yes            | — |
| `--prefix`        | Output files prefix string                                                                           | Yes            | — |
| `-g, --gtf`       | Path to the GTF annotation file compatible with the BAM file                                         | Yes            | — |
| `-f, --fasta`     | Path to the FASTA genome file compatible with the BAM file                                           | Yes            | — |
| `-s, --sam`       | Path to the SAM/BAM file or directory of files                                                       | Yes (if BAM)   | — |
| `--mate_1`        | Path to FASTQ file of mate 1 (for paired-end reads)                                                  | Yes (if FASTQ) | — |
| `--mate_2`        | Path to FASTQ file of mate 2 (for paired-end reads)                                                  | Yes (if FASTQ) | — |
| `--gsnap_idx`     | Path to GSNAP index. If omitted in `GSNAP+TSJS` mode, one is auto-built from `--fasta`+`--gtf`.      | No             | auto-build at `<out>/gsnap_idx/` |
| `--gsnap_idx_flags` | Extra flags passed to `gmap_build` and `gsnap`                                                     | No             | `-d index -c transcriptome` |
| `--gsnap_aln_flags` | Extra flags passed to `gsnap` at alignment time. See [GSNAP accelerated parameters](#gsnap-accelerated-parameters) below for GPU-acceleration knobs you can wire in here. | No             | `--gunzip -A sam --novelsplicing 1` |
| `-c, --config`    | Config `.json` file to control DeepSAP internal parameters                                           | No             | `/scripts/parameters_config.json` |
| `--batch`         | Number of candidate splice junctions scored per transformer forward pass. Larger values raise throughput but increase GPU memory use (see [Requirements](#requirements) for a memory-vs-batch reference). | No             | `2048` |
| `--score_method`  | Scoring method applied to junctions and reads: `Add` or `Multiply`                                      | No             | `Add` |
| `--no-fp16`       | Don't use fp16 half-precision floating-point                                                         | No             | fp16 enabled |
| `--set_size`      | Set size to split datasets for inference                                                             | No             | `102400` (= 1024 × 100) |
| `-t, --threads`   | Number of threads                                                                                    | No             | host `os.cpu_count()` |
| <a id="gsnap-accelerated-parameters"></a> `--localdb-batch`   | [GSNAP accelerated, passed through to `gsnap` only if set] Requests packed into each GPU kernel launch on the accelerated `--localdb=GPU` path. | No | unset (gsnap default `24000`) |
| `--localdb-scratch` | [GSNAP accelerated, passed through to `gsnap` only if set] Unified GPU device-byte budget for localdb scratch (accepts `K`/`M`/`G` suffixes, e.g. `8G`). | No | unset (gsnap default `12G`) |
| `--batch-nreads`    | [GSNAP accelerated, passed through to `gsnap` only if set] Max individual reads per frame. Paired-end input requires an even value ≥ 2. | No | unset (gsnap default `250`) |
---------------------------------------------------------------------------------------------------------------------------------------------
<br>

## Version History
### v0.1.0
* Added GPU-accelerated GSNAP v2025-04-19. The runtime image now ships a CUDA-accelerated GSNAP build with both Stage-1 (r2d) and Stage-2 (localdb) running on the GPU by default; tunable passthrough knobs are exposed via `--localdb-batch`, `--localdb-scratch`, and `--batch-nreads` (see [Command-line Arguments](#command-line-arguments)).
* Added `--mode` flag to explicitly select pipeline mode (`index`, `GSNAP+TSJS`). The default `GSNAP+TSJS` preserves the v0.0.x end-to-end behaviour, including auto-building a GSNAP index when `--gsnap_idx` is omitted.
* Bug fix: output BAM is now correctly suffixed with `.bam`.
* Bug fix: SAM records with empty CIGAR strings are now normalised to `*` before being written to the BAM stream.
* Bug fix: stricter logits-shape validation in `predict.py` (previously masked by a bare `except`).
### v0.0.3
* Fixed key error in parsing FASTA files.
* Fixed gene_id pattern error in parsing GTF files. 
### v0.0.2
* Updated GSNAP aligner to version `2025-04-19`.
### v0.0.1
* Initial release.

## License/Terms of Use
By pulling and using the Parabricks DeepSAP container, you accept the governing terms: The software and materials are governed by the NVIDIA Software License Agreement (found at https://www.nvidia.com/en-us/agreements/enterprise-software/nvidia-software-license-agreement/) and the Product-Specific Terms for NVIDIA AI Products (found at https://www.nvidia.com/en-us/agreements/enterprise-software/product-specific-terms-for-ai-products/); except for the model which is governed by the NVIDIA Models Community License Agreement(found at NVIDIA Community Model License). ADDITIONAL INFORMATION: Apache 2.0.
