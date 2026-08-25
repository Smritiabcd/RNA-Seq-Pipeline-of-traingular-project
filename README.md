# Bee Virome Pipeline

Paired-end read QC, trimming, sequential host-read removal, and viral alignment for a honeybee/Varroa virome study — three host organisms (worker bees, queen bees, and *Varroa destructor* mites), 14 samples across two hives.

## What it does

For every sample in `raw_fastq/`, the pipeline runs:

```
Raw reads
   │
   ├─▶ FastQC (raw)                         → qc/fastqc_raw/
   │
   ▼
fastp (adapter trim, quality filter, Q20, min length 50bp)
   │                                        → trimmed/
   │                                        → qc/fastp_reports/
   ▼
FastQC (trimmed)                            → qc/fastqc_trimmed/
   │
   ▼
Host removal (BWA-MEM, sequential subtraction)
   ├─ 1. Apis cerana genome     → host_removed/*_no_cerana_*
   ├─ 2. Apis mellifera genome  → host_removed/*_no_mellifera_*
   └─ 3. Varroa genome          → host_removed/*_clean_*
   │
   ▼
Alignment of clean reads → target virus (BWA-MEM)
   │                                        → aligned/newvirus/*.bam
   ▼
flagstat / coverage / depth per sample      → aligned/newvirus/*_flagstat.txt
                                               aligned/newvirus/*_coverage.txt
                                               aligned/newvirus/*_depth.txt

Once all samples are done:
MultiQC aggregates raw FastQC + trimmed FastQC + fastp   → qc/multiqc_final/multiqc_final_report.html
```

Host reads are removed **sequentially, not in parallel** — each species' unmapped reads feed into the next alignment. This matters because *Apis mellifera* and *Apis cerana* references can cross-map some reads; removing cerana first, then mellifera, then Varroa, avoids double-counting and keeps only reads unmapped to all three hosts as "clean" (candidate viral) reads.

## Requirements

| Tool | Used for |
|---|---|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | raw & trimmed read quality reports |
| [fastp](https://github.com/OpenGene/fastp) | adapter/quality trimming |
| [BWA](https://github.com/lh3/bwa) | host removal + viral alignment (`bwa mem`) |
| [samtools](https://www.htslib.org/) | BAM sorting/filtering, flagstat, coverage, depth |
| [MultiQC](https://multiqc.info/) | aggregated QC report |

Install via conda:

```bash
conda install -c bioconda fastqc fastp bwa samtools multiqc -y
```

## Input layout

Paired-end files named `<sample>_1.fq.gz` / `<sample>_2.fq.gz` in `raw_fastq/`:

```
raw_fastq/
├── Honeybee-C-1_1.fq.gz   Honeybee-C-1_2.fq.gz
├── Honeybee-C-2_1.fq.gz   Honeybee-C-2_2.fq.gz
├── Honeybee-D-1_1.fq.gz   Honeybee-D-1_2.fq.gz
├── Queenbee-E-1_1.fq.gz   Queenbee-E-1_2.fq.gz
├── Varroa-A-1_1.fq.gz     Varroa-A-1_2.fq.gz
└── ...
```

The sample name is derived automatically from the `_1.fq.gz` suffix — no sample sheet needed. Any file missing its `_2` mate is skipped with a warning.

## Reference genomes

Place these in `references/` (or edit the paths at the top of the script):

```
references/
├── apis_cerana_genome.fa
├── apis_mellifera.fa
├── varroa_genome.fa
└── newvirus.fasta
```

BWA indices are built automatically on first run if not already present (`.bwt` file check), and skipped on subsequent runs.

## Usage

```bash
chmod +x run_full_pipeline.sh
./run_full_pipeline.sh
```

Edit the config block at the top of the script to match your setup:

```bash
PROJECT_DIR="/mnt/d/bee_immune_analysis"
THREADS=8            # bwa/fastp/samtools threads
FASTQC_THREADS=4      # parallel FastQC file handles
```

## Output structure

```
qc/
├── fastqc_raw/              # per-file HTML+zip, raw reads
├── fastqc_trimmed/          # per-file HTML+zip, trimmed reads
├── fastp_reports/           # per-sample HTML+JSON trimming reports
└── multiqc_final/
    └── multiqc_final_report.html   ← combined QC report (start here)

trimmed/                     # fastp-trimmed reads
host_removed/                # intermediate + final host-subtracted reads
  ├── *_no_cerana_[12].fq.gz
  ├── *_no_mellifera_[12].fq.gz
  └── *_clean_[12].fq.gz     ← final non-host reads, input to viral alignment

aligned/newvirus/
  ├── *_vs_newvirus.bam(.bai)
  ├── *_flagstat.txt         # mapping summary
  ├── *_coverage.txt         # per-reference % breadth + mean depth
  └── *_depth.txt            # per-base depth

logs/                        # stderr from every tool call, per sample
```

## Interpreting virus hits

Reads are called as a positive detection when, per sample:

- **≥10% genome coverage breadth**, and
- **≥10 mapped reads** 

These are printed live during the run and also derivable from `aligned/newvirus/*_coverage.txt`:

```bash
awk 'NR>1 && $6>10 && $7>1' aligned/newvirus/*_coverage.txt
```

