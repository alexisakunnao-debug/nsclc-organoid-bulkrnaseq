#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo 'Usage: worker.sbatch STAGE PROJECT_ROOT SAMPLES_TSV WORK_DIR' >&2
    exit 2
fi
stage="$1"
root="$2"
samples="$3"
work_dir="$4"
threads="${SLURM_CPUS_PER_TASK:-8}"
args=("$stage" --samples "$samples" --work-dir "$work_dir" --threads "$threads")
if [[ "$stage" == quant ]]; then
    args+=(--task "${SLURM_ARRAY_TASK_ID:?Missing array task ID}")
fi
bash "$root/run.sh" "${args[@]}"
