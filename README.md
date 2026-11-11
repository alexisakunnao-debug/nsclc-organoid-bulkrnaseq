# NSCLC organoid RNA-seq quality and quantification

I built this workflow to look at public NSCLC organoid RNA-seq data from [SRA study SRP158722](https://www.ncbi.nlm.nih.gov/sra?term=SRP158722) and [GEO series GSE119004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119004).

The goal is to check the reads with FastQC and MultiQC, then use Salmon to estimate transcript abundance. I wrote Bash entry points for running the steps on an HPC cluster and Python code for checking the sample metadata, downloads and results.

**Current status:** the workflow is in place, and its software checks have been run on small synthetic examples. I have not run FastQC or Salmon on the study reads yet. Reference indexing and the full analysis still need an HPC run, so there are no mapping rates or expression results to report.

## Cohort

I found nine GEO records marked as organoids. Each is linked to one paired-end RNA-seq run through its SRA experiment and BioSample identifiers. The compressed reads total **35.55 GiB**.

I chose three samples for the pilot. Together, their reads total **7.15 GiB**.

| GEO sample | Submitted title | Run |
|---|---|---|
| GSM3356477 | PDO54 | SRR7749979 |
| GSM4063806 | LPTO85_org | SRR10075656 |
| GSM4063812 | XDO137 | SRR10075662 |

I picked these to try the workflow on a smaller part of the organoid set first. This selection does not show that samples came from the same patient or establish a biological replicate design.

## What the workflow does

`mermaid
flowchart LR
    A[GEO and ENA metadata] --> B[Validated sample sheet]
    C[GENCODE 47 transcripts and GRCh38 genome] --> D[Decoy-aware Salmon index]
    B --> E[Checksum and paired-read checks]
    E --> F[FastQC]
    F --> G[Salmon quantification]
    D --> G
    G --> H[MultiQC and cohort tables]
    F --> H
`

In order, the workflow:

1. Uses the pinned reference URLs and checksums in `config/reference.json`.
2. Builds a Salmon index from GENCODE 47 transcripts and GRCh38 primary-genome decoys.
3. Checks each compressed FASTQ against the ENA checksum and file size, then checks that paired reads match.
4. Runs FastQC on the untrimmed reads.
5. Runs Salmon with automatic library detection, sequence and GC bias correction, and 30 bootstrap replicates.
6. Checks the Salmon summary and abundance tables, then creates cohort tables and a MultiQC report.

I kept read trimming out of the automatic steps so I can first review the actual quality reports. The workflow also saves checksums for completed stages, so a later run can check whether an earlier result is safe to reuse.

## Set up on Linux or an HPC cluster

The environment file pins FastQC 0.12.1, Salmon 1.10.2, MultiQC 1.25.1 and Python 3.11.11. It still needs to be solved and checked on the Linux system where the analysis will run.

`bash
conda env create -f environment.yml
conda activate nsclc-rnaseq
python3 -m unittest discover -s tests -v
bash run.sh preflight --require-tools
`

Before starting, run the preflight step to see the required tools and available storage. Current planning estimates are about **68 GiB for the pilot** and **96 GiB for the full cohort**, plus space for the software environment and working room. The Slurm profile requests 48 GiB of memory for reference indexing; actual memory use has not been measured.

Use shared storage that all jobs can reach. Check your cluster's scratch cleanup rules and resource limits before submitting jobs.

## Submit the pilot to Slurm

Activate the environment first. Set the account and partition only if your cluster requires them.

`bash
export SLURM_ACCOUNT=your_account
export SLURM_PARTITION=your_partition

bash submit-slurm.sh config/pilot.tsv /shared/scratch/your_user/nsclc-rnaseq
`

This submits the reference step, a quantification array with up to two workers, and a report step. Each step waits for its prerequisites to finish successfully. No HPC jobs have been submitted yet.

To use all nine runs, replace `config/pilot.tsv` with `config/cohort.tsv` and adjust the time, memory and concurrency settings in `submit-slurm.sh` for your cluster.

## Run steps individually

`bash
bash run.sh reference --work-dir /shared/scratch/your_user/nsclc-rnaseq
bash run.sh quant --task 1 --samples config/pilot.tsv --work-dir /shared/scratch/your_user/nsclc-rnaseq
# Run tasks 2 and 3 before making the report
bash run.sh summarize --samples config/pilot.tsv --work-dir /shared/scratch/your_user/nsclc-rnaseq
`

`--task` is the one-based row number in the sample sheet. Use the same work directory and configuration for each step.

## Expected outputs

| Location under the work directory | Contents |
|---|---|
| `reference/<signature>/` | Reference FASTA, decoys, Salmon index and provenance |
| `quant/<sample>/<signature>/fastqc/` | FastQC HTML and ZIP reports for each mate |
| `quant/<sample>/<signature>/salmon/` | `quant.sf`, inferential replicates, library format and mapping metadata |
| `quant/<sample>/<signature>/metrics.json` | Abundance and fragment checks, with review flags |
| `reports/<signature>/multiqc_report.html` | Cohort quality summary |
| `reports/<signature>/sample_qc.tsv` | Run IDs, mapping metrics and FastQC status counts |
| `reports/<signature>/transcript_tpm.tsv` | Transcript TPM matrix |
| `reports/<signature>/transcript_estimated_counts.tsv` | Fractional estimated fragment-count matrix |
| `complete.json` within each completed stage | Input/configuration provenance and output checksums |

## Checks and interpretation

The tests cover metadata joins, incomplete or corrupt downloads, mismatched read pairs, reference IDs, numerical output checks, checkpoints and stage handoffs. The stage test uses small stand-in tools; it does not run FastQC or Salmon on study reads. See [validation status](docs/validation_status.md) for what has and has not been run.

FastQC warnings and low mapping rates should be reviewed in context. Duplication or sequence-composition warnings can have more than one explanation, including the library preparation. The thresholds in `config/qc.json` only flag samples for review; they do not remove samples automatically.

TPM is useful for looking at relative transcript abundance within a sample. Estimated counts can be fractional. Comparing samples for differential expression would need a suitable experimental design and a separate statistical analysis. A human reference alone also cannot assess mouse contamination in xenograft-derived material.

See [methods and decisions](docs/methods.md) for more detail on the cohort, reference and interpretation. The repository contains public metadata and project files only; raw reads, reference files, generated results and my private study notes are kept out.

## Refresh the metadata

`bash
python3 workflow/discover.py --out /path/to/a/review-copy
`

Review the new tables before replacing the pinned sample sheets. The discovery script checks for the expected nine organoid records, confirms the BioSample links and stops if a sample has an unexpected number of runs. Source URLs and SHA-256 checksums are saved in `data/metadata/provenance.json`.
