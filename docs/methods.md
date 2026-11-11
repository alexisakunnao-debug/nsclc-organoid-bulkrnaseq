# Methods and decisions

## Choosing the samples

I started with the samples in GSE119004 whose `sample type` characteristic says `organoid`. I matched their SRA experiment and BioSample IDs to the ENA file report for SRP158722, which provides paired FASTQ links, file sizes and MD5 checksums. The discovery step checks that the experiment, BioSample, project PRJNA487654, human organism, RNA-seq strategy and paired-end layout agree.

The saved cohort has nine samples, with one run per sample. I keep the submitted titles as annotations, but do not use numbers in the titles to guess patient identity or connect models. If I find a sample with multiple runs later, I will need to decide how to handle lanes and runs before adding it.

I chose GSM3356477, GSM4063806 and GSM4063812 for the pilot. They are part of the documented organoid set and let me try the workflow on three samples first. I have not assigned a treatment comparison or biological replicate structure.

## Reference and Salmon settings

I pinned GENCODE release 47 transcript sequences and the matching GRCh38 primary genome assembly for decoys. Release 47 is associated with GRCh38.p14. The exact URLs and provider MD5 checksums are in `config/reference.json`. I chose this fixed release so the reference stays the same when the workflow is rerun. [GENCODE release 47](https://www.gencodegenes.org/human/release_47.html)

The reference builder combines the transcript and genome sequences. Genome sequence IDs are written to `decoys.txt`; empty references and duplicate IDs stop the build. The index uses Salmon's GENCODE header handling and k=31.

Decoys let genomic sequences compete with annotated transcripts during selective alignment, which can reduce assignments to similar genomic regions that are not annotated as transcripts. Salmon uses mapping validation, automatic library detection, sequence and GC bias correction, and 30 bootstrap replicates. The bootstraps estimate uncertainty under the model; they are not biological replicates. [Salmon documentation](https://salmon.readthedocs.io/en/stable/salmon.html)

The worker checks software versions before running. The environment file pins package versions, but it does not lock every platform-specific build or transitive dependency. After solving the environment on HPC, I should save an explicit environment export with the run records.

## Read checks

Before reusing a download, the workflow checks its compressed byte count and ENA MD5 checksum. New files are written to a temporary path and moved into place only after those checks pass.

Both FASTQ files are streamed together from beginning to end. The checks cover record structure, sequence and quality lengths, nucleotide symbols, quality character range, and matching read names. Unequal file lengths, mismatched pairs or empty data stop the run. Streaming also checks gzip integrity. These checks confirm file structure; FastQC is used to review sequencing quality.

I run FastQC on the untrimmed reads first. Trimming is not automatic because I want to review adapter and quality reports before deciding whether it is needed. If I add a trimmed workflow later, I will keep the original reports and record the trimming settings. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Quantification checks

For each sample, Salmon's `num_processed` should match the input fragment count. Mapped fragments must be between zero and the processed count. The reported mapping percentage is checked against that denominator, allowing a difference of 0.1 percentage points.

The transcript table is checked for the expected columns, unique IDs, and finite non-negative values. When fragments map, the TPM total should be close to one million (within 100, allowing for rounding). This is a file consistency check; it does not show that the result is biologically correct. Fractional estimated fragment counts are retained.

All samples need the same transcript IDs in the same order and must use the same reference signature. Transcript version numbers are kept.

## Reviewing samples

The workflow flags mapping below 70% or fewer than one million processed fragments for review. These are project settings in `config/qc.json`, not clinical cutoffs. They do not automatically exclude a sample. The report also counts FastQC WARN and FAIL modules for each mate.

Low mapping can have several causes, such as contamination, a reference mismatch, adapters or low read quality. A high mapping rate does not confirm sample identity. If xenograft-derived material is involved, I may need an additional human/mouse check. Duplication, GC content and sequence bias should be considered alongside the library protocol and sample context.

MultiQC only receives completed FastQC and Salmon outputs from the selected cohort, so reports from earlier or unrelated work are less likely to be mixed in. [MultiQC documentation](https://docs.seqera.io/multiqc/)

## Keeping and resuming work

Reference and sample signatures include the configuration and pipeline code checksum. Sample signatures also include the run's file URLs, provider checksums and byte counts. Completed outputs are checked by SHA-256 before reuse. A configuration change creates a separate output signature.

Each step writes to its own temporary build directory. Once required files pass checks, the workflow writes a completion manifest and moves the directory into place. Failed builds keep their logs but are not treated as completed results.

Directory locks prevent two jobs from writing to the same step at once. If a job stops unexpectedly, check the Slurm status and lock owner record before removing a stale lock. If a completed output was changed, investigate it or use a new work directory.

All jobs need access to the same shared work directory. The Slurm profile currently requests 8 CPUs and 48 GiB for indexing, 8 CPUs and 32 GiB per quantification task, and 2 CPUs and 8 GiB for reporting. At most two quantification tasks run together. These are starting requests; I have not measured actual runtime or peak memory on the target cluster.

## Scope

This project is for read quality checks and transcript quantification. I have not defined a differential-expression comparison, clinical endpoint or patient-to-model pairing. I have not reported mapping rates, expression findings or sample exclusions because the real analysis has not run yet.
