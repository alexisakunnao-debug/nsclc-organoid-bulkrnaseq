# First sample FastQC check

I ran FastQC on the paired reads for one sample from the public NSCLC organoid study, SRP158722. The sample is GSM3356477, run SRR7749979. The two compressed read files matched the size and MD5 values in the sample sheet before analysis.

## What I ran

FastQC 0.12.1 analyzed Read 1 and Read 2 separately. Each file contained 31,912,606 reads, with a read length of 75 bases. The output is available here:

- [Read 1 report](../reports/fastqc/SRR7749979_1_fastqc.html)
- [Read 2 report](../reports/fastqc/SRR7749979_2_fastqc.html)

## What stood out

Both mates passed the per-base quality, per-sequence quality, adapter-content, ambiguous-base, and sequence-length checks. FastQC reported zero reads flagged as poor quality.

Read 1 failed the per-base sequence-content and sequence-duplication checks. Read 2 had warnings for those checks. Both mates also had warnings for per-sequence GC distribution and overrepresented sequences.

These flags are prompts to investigate, not automatic reasons to discard the reads. RNA-seq libraries can have uneven sequence composition and repeated reads because transcript abundance is uneven. I would compare these patterns with the other samples and the library-preparation details before deciding whether they are expected here.

## Scope

This is a raw-read quality check for one sample. I have not run MultiQC, trimmed reads, aligned them to a genome, or quantified transcripts. FastQC results alone do not establish that the sample is biologically valid or suitable for downstream analysis.
