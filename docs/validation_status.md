# Validation status

Last updated 10 November 2026 (America/Toronto).

| Component | What has been checked |
|---|---|
| Public run discovery | ENA metadata retrieved and matched to the saved GEO records |
| Cohort manifest | Nine organoid runs, three in the pilot; paired reads linked by experiment and BioSample |
| Compressed FASTQ inventory | 35.55 GiB for the cohort; 7.15 GiB for the pilot |
| Reference files | GENCODE 47 URLs and provider checksums checked against the release checksum list |
| Bash entry points | Syntax checks passed |
| Python tests | 31 tests passed with the local Python 3.14.7 interpreter |
| Stage handoffs and resume | Synthetic orchestration test passed using stand-in external tools |
| Conda environment solve | Still needs to be checked on the target Linux system |
| FastQC on study reads | Not run |
| Salmon reference index | Not built |
| Salmon on study reads | Not run |
| MultiQC on study outputs | Not run |
| Slurm submission and resource measurements | Not run |

The tests check workflow behavior and error handling. They do not show that FastQC or Salmon has run successfully on this study, or that the Slurm resource requests fit a particular cluster.

## Re-running the local checks

```bash
python3 -B -m unittest discover -s tests -v
bash -n run.sh submit-slurm.sh workflow/worker.sbatch
bash run.sh preflight
```

The standard-library tests do not need bioinformatics tools or network access. They use small synthetic files that are temporary and are not included in the study reports.

## First HPC run

1. Solve and activate the pinned environment on the target Linux system.
2. Run the tests and check the expected tool versions with preflight.
3. Confirm shared storage, scratch retention, scheduler account and resource limits.
4. Submit the three-sample pilot. Review all FastQC modules, inferred library types and mapping flags.
5. Check `sample_qc.tsv`, abundance tables, MultiQC and stage completion manifests.
6. Record elapsed time, peak memory, software builds and any decisions before running the rest of the cohort.

I will update this page after the real pilot run. Any measured results should come from reviewed output reports.
