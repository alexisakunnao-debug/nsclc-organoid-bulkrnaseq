#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
samples="${1:-$ROOT/config/pilot.tsv}"
work_dir="${2:-$ROOT/work}"
samples="$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$samples")"
work_dir="$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$work_dir")"
bash "$ROOT/run.sh" preflight --samples "$samples" --work-dir "$work_dir" --require-tools
command -v sbatch >/dev/null || { echo 'Slurm sbatch is unavailable.' >&2; exit 1; }
count="$(python3 "$ROOT/workflow/rnaseq.py" count --samples "$samples")"
mkdir -p "$ROOT/logs"
# Set SLURM_ACCOUNT or SLURM_PARTITION only if required by your cluster.
common=(--parsable --chdir="$ROOT" --output="$ROOT/logs/%x-%A_%a.log")
[[ -z "${SLURM_ACCOUNT:-}" ]] || common+=(--account="$SLURM_ACCOUNT")
[[ -z "${SLURM_PARTITION:-}" ]] || common+=(--partition="$SLURM_PARTITION")
ref_job="$(sbatch "${common[@]}" --job-name=nsclc-reference --cpus-per-task=8 --mem=48G --time=08:00:00 "$ROOT/workflow/worker.sbatch" reference "$ROOT" "$samples" "$work_dir")"
ref_job="${ref_job%%;*}"
quant_job="$(sbatch "${common[@]}" --job-name=nsclc-quant --dependency="afterok:$ref_job" --array="1-$count%2" --cpus-per-task=8 --mem=32G --time=12:00:00 "$ROOT/workflow/worker.sbatch" quant "$ROOT" "$samples" "$work_dir")"
quant_job="${quant_job%%;*}"
summary_job="$(sbatch "${common[@]}" --job-name=nsclc-report --dependency="afterok:$quant_job" --cpus-per-task=2 --mem=8G --time=02:00:00 "$ROOT/workflow/worker.sbatch" summarize "$ROOT" "$samples" "$work_dir")"
printf 'Reference job %s\nQuantification array %s\nReport job %s\n' "$ref_job" "$quant_job" "$summary_job"
