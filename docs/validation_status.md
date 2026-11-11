# Validation status

Recorded on 10 November 2026, America/Toronto.

| Component | Current evidence |
|---|---|
| Public run discovery | ENA metadata retrieved and joined to the saved GEO source |
| Cohort manifest | Nine organoid runs, three in the pilot; all paired and linked by experiment and BioSample |
| Compressed FASTQ inventory | 35.55 GiB for the cohort; 7.15 GiB for the pilot |
| Reference pinning | GENCODE 47 URLs and provider checksums verified from the release checksum list |
| Bash entry points | Syntax checks passed |
| Python tests | 31 tests passed with Python 3.14.7 with the local Python interpreter |
| Stage handoffs and resume | Synthetic orchestration test passed with explicit external-tool stubs |
| Conda environment solve | Pending on the target Linux system |
| FastQC on study reads | Pending |
| Full reference index | Pending |
| Salmon on study reads | Pending |
| MultiQC on study outputs | Pending |
| Slurm submission and resource measurements | Pending |

The tests validate workflow contracts and failure handling. They do not establish that the external tools have run successfully on this study or that the Slurm profile matches a particular cluster.

## Reproduce local checks

```bash
python3 -B -m unittest discover -s tests -v
bash -n run.sh submit-slurm.sh workflow/worker.sbatch
bash run.sh preflight
```

The standard-library tests need no bioinformatics executables and make no external network requests. Their small synthetic files are temporary and never enter study reports.

## First HPC run

1. Solve and activate the pinned environment on the target Linux system.
2. Run the tests and require the expected tool versions in preflight.
3. Confirm shared storage capacity, scratch retention, scheduler account and resource limits.
4. Submit the three-run pilot, then inspect all FastQC modules, inferred library types and mapping flags.
5. Review `sample_qc.tsv`, the abundance matrices, MultiQC and stage completion manifests.
6. Record elapsed time, peak memory, software builds and any interpretation decisions before expanding the cohort.

Keep this page in sync with actual completed runs. Add measured results only after reviewing the generated reports.
