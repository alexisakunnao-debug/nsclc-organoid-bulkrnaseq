# Methods and decisions

## Source selection

The cohort begins with sample records in GSE119004 whose `sample type` characteristic explicitly states `organoid`. Experiments and BioSample identifiers come from the submitted database relations. The ENA file report for SRP158722 provides paired FASTQ locations, file sizes and MD5 checksums. Discovery requires agreement between the experiment, BioSample, project PRJNA487654, human organism, RNA-seq strategy and paired layout.

The pinned cohort contains nine records and one run per record. Titles are retained as annotations. Title digits are never used to infer patient identity or model relationships. A future experiment with multiple sequencing runs requires a reviewed lane/run aggregation policy before entering this cohort.

Pilot samples GSM3356477, GSM4063806 and GSM4063812 provide three records from the documented organoid set. They support workflow evaluation, with no treatment contrast or inferred replicate design.

## Reference and quantification

GENCODE release 47 supplies the chromosome transcript sequences and the GRCh38 primary genome assembly used as decoys. Release 47 is associated with GRCh38.p14; the chosen genome FASTA specifically covers the primary assembly. The exact file URLs and provider MD5 values are pinned in `config/reference.json`. The release is a reproducible fixed choice, independent of which release is newest when the project is rerun. [GENCODE release 47](https://www.gencodegenes.org/human/release_47.html)

Genome FASTA records follow the transcript records in the combined reference. Their sequence identifiers form `decoys.txt`. The builder rejects empty references and duplicate identifiers, and uses Salmon's GENCODE header handling and k=31.

Genome decoys allow genomic sequences to compete with annotated transcripts during selective alignment. This reduces spurious transcript assignments from similar unannotated genomic regions. Salmon is run with mapping validation, automatic library detection, sequence bias correction, GC bias correction and 30 bootstrap replicates. Bootstraps represent inferential uncertainty conditional on the data and model; they do not supply biological replicates. [Salmon documentation](https://salmon.readthedocs.io/en/stable/salmon.html)

The reference and cohort software are versioned together. The worker checks executable versions before computation. The environment specification pins package versions but does not lock platform-specific builds or all transitive dependencies. A successful first HPC run should retain an explicit environment export alongside the run records.

## Read checks

Downloads must match both expected compressed byte counts and ENA MD5 checksums. A partial download is written to a separate temporary path and published only after verification. An existing file is verified before reuse.

Paired FASTQ files are read together from beginning to end. Record structure, sequence/quality lengths, nucleotide symbols, quality-character range and normalized read names are checked. Unequal record counts, mismatched pairs and empty data stop processing. Full streaming also exercises gzip integrity checks. This validates structural compatibility while FastQC assesses sequencing quality.

FastQC is applied to the untrimmed reads. Trimming is not automatic because it should follow inspection of actual adapter and quality evidence. A future trimmed branch should retain pre-trimming reports and document its parameters. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Quantification checks

Salmon's `num_processed` must equal the input fragment count. Mapped fragments must fall between zero and processed fragments. The reported mapping percentage must agree with that denominator within 0.1 percentage points.

The transcript table must contain the expected columns, unique transcript identifiers and finite non-negative numerical values. When fragments map, the TPM sum must be within 100 of one million, allowing output rounding. This is a numerical consistency check, independent of biological validity. Fractional estimated fragment counts are retained.

The cohort matrices require identical transcript identifiers and ordering across samples. Each sample was quantified against the same reference signature. Transcript versions are preserved.

## Sample review

Mapping below 70% and fewer than one million processed fragments trigger review flags. These thresholds are project conventions, configurable in `config/qc.json`. They are not clinical acceptance criteria. The workflow also reports per-mate counts of FastQC WARN and FAIL modules.

Low mapping can have several causes, including contamination, reference mismatch, adapters or poor-quality reads. A high mapping rate does not establish sample identity. Xenograft-derived organoids may require an additional human/mouse assessment. RNA-seq duplication, GC content and sequence bias should be interpreted with the library protocol and sample context.

MultiQC aggregates only the selected cohort's completed FastQC and Salmon directories, reducing the chance that reports from an earlier or unrelated run enter the summary. [MultiQC documentation](https://docs.seqera.io/multiqc/)

## Execution and recovery

Reference and sample signatures include the configuration and pipeline code checksum. Sample signatures additionally include the run's file URLs, provider checksums and byte counts. Completed outputs are checked by SHA-256 before reuse. A configuration change produces a separate output signature.

Each stage writes into a uniquely named build directory. The completion manifest is written after required output checks, and the directory is then renamed into its final location. Failed build directories retain logs for diagnosis. They are never used as completed results.

Directory locks prevent two workers from writing the same stage. A terminated process can leave a stale lock. Check Slurm job state and the lock's owner record before manually removing it. Modified completed outputs produce a clear error and require investigation or a new work directory.

All jobs need the same shared work directory. The Slurm profile requests 8 CPUs/48 GiB for indexing, 8 CPUs/32 GiB per quantification task, and 2 CPUs/8 GiB for reporting. The quantification array allows two concurrent workers. These are starting resource requests; runtime and peak memory require measurement on the target system.

## Scope

The workflow produces read-quality and transcript-quantification outputs. It does not define a differential-expression contrast, clinical endpoint or patient-to-model pairing. No mapping rates, expression findings or sample exclusions have been claimed before the real run.
