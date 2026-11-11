"""Run-level RNA-seq QC and Salmon quantification for SRP158722."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from itertools import zip_longest
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
import uuid
import zipfile

ROOT = Path(__file__).resolve().parents[1]
VERSIONS = {"fastqc": "0.12.1", "salmon": "1.10.2", "multiqc": "1.25.1"}
FIELDS = ["sample_id", "geo_accession", "experiment_accession", "run_accession",
          "read1_url", "read2_url", "read1_md5", "read2_md5", "read1_bytes", "read2_bytes"]


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def digest(path, algorithm="sha256"):
    checksum = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:20]


def read_samples(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != FIELDS:
            raise ValueError(f"Sample sheet columns must be {FIELDS}")
        rows = list(reader)
    if not rows:
        raise ValueError("The sample sheet is empty")
    runs, samples = set(), set()
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("Malformed sample sheet row")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row["sample_id"]):
            raise ValueError("Unsafe sample ID")
        if row["sample_id"] in samples or row["run_accession"] in runs:
            raise ValueError("Duplicate sample or run ID; use one row per run")
        samples.add(row["sample_id"])
        runs.add(row["run_accession"])
        for field, prefix in [("geo_accession", "GSM"), ("experiment_accession", "SRX"), ("run_accession", "SRR")]:
            if not re.fullmatch(prefix + r"\d+", row[field]):
                raise ValueError(f"Invalid {field}")
        for mate in (1, 2):
            url = urllib.parse.urlsplit(row[f"read{mate}_url"])
            if url.scheme != "https" or url.netloc != "ftp.sra.ebi.ac.uk":
                raise ValueError("FASTQ downloads must use the ENA HTTPS host")
            expected_name = f"{row['run_accession']}_{mate}.fastq.gz"
            if Path(url.path).name != expected_name:
                raise ValueError("Run accession or mate disagrees with FASTQ URL")
            if not re.fullmatch(r"[0-9a-f]{32}", row[f"read{mate}_md5"]):
                raise ValueError("Invalid MD5")
            if int(row[f"read{mate}_bytes"]) <= 0:
                raise ValueError("FASTQ byte count must be positive")
    return rows


def tool_versions():
    found = {}
    for tool, expected in VERSIONS.items():
        if shutil.which(tool) is None:
            raise RuntimeError(f"Missing {tool}; activate the environment from environment.yml")
        result = subprocess.run([tool, "--version"], check=True, capture_output=True, text=True)
        observed = (result.stdout + result.stderr).strip()
        if re.search(r"(?<![\d.])" + re.escape(expected) + r"(?![\d.])", observed) is None:
            raise RuntimeError(f"Expected {tool} {expected}, observed {observed}")
        found[tool] = observed
    return found


def nearest_existing(path):
    path = Path(path).resolve()
    while not path.exists():
        path = path.parent
    return path


def require_space(path, needed_bytes):
    free = shutil.disk_usage(nearest_existing(path)).free
    if free < needed_bytes:
        raise RuntimeError(f"Insufficient storage: need {needed_bytes / 2**30:.1f} GiB free, have {free / 2**30:.1f} GiB")


@contextmanager
def lock(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.mkdir()
    except FileExistsError as exc:
        raise RuntimeError(f"Stage is locked: {path}. Check the owning job before removing a stale lock.") from exc
    try:
        write_json(path / "owner.json", {"pid": os.getpid(), "started_utc": datetime.now(timezone.utc).isoformat()})
        yield
    finally:
        (path / "owner.json").unlink(missing_ok=True)
        path.rmdir()


def verify_file(path, expected_md5, expected_bytes=None):
    path = Path(path)
    if expected_bytes is not None and path.stat().st_size != expected_bytes:
        raise ValueError(f"Wrong file size: {path.name}")
    if digest(path, "md5") != expected_md5:
        raise ValueError(f"Checksum mismatch: {path.name}")


def download(url, target, expected_md5, expected_bytes=None):
    """Publish a download only after its checksum matches the metadata."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        verify_file(target, expected_md5, expected_bytes)
        return
    temporary = target.with_name(target.name + ".part-" + uuid.uuid4().hex)
    try:
        with urllib.request.urlopen(url, timeout=120) as response, temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
        verify_file(temporary, expected_md5, expected_bytes)
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)


def run_command(argv, log):
    with Path(log).open("w") as handle:
        handle.write(json.dumps(argv) + "\n")
        handle.flush()
        subprocess.run([str(part) for part in argv], stdout=handle, stderr=subprocess.STDOUT, check=True)


def checked_stage(path, signature):
    path = Path(path)
    if not path.exists():
        return False
    manifest = read_json(path / "complete.json")
    if manifest["signature"] != signature:
        raise ValueError(f"Checkpoint signature mismatch: {path}")
    for name, checksum in manifest["output_sha256"].items():
        if digest(path / name) != checksum:
            raise ValueError(f"Completed output was modified: {path / name}")
    return True


def finish_stage(build, final, signature, provenance):
    files = {str(path.relative_to(build)): digest(path) for path in sorted(build.rglob("*")) if path.is_file()}
    if not files:
        raise ValueError("Cannot complete an empty stage")
    write_json(build / "complete.json", {
        "signature": signature, "completed_utc": datetime.now(timezone.utc).isoformat(),
        "output_sha256": files, "provenance": provenance,
    })
    build.rename(final)


def reference_signature(root):
    return fingerprint({"reference": read_json(root / "config/reference.json"), "pipeline_sha256": digest(root / "workflow/rnaseq.py")})


def build_reference(root, work, threads):
    config = read_json(root / "config/reference.json")
    versions = tool_versions()
    signature = reference_signature(root)
    final = work / "reference" / signature
    with lock(work / "locks/reference"):
        if checked_stage(final, signature):
            print(f"Verified reference checkpoint {signature}")
            return
        require_space(work, 45 * 2**30)
        build = final.with_name(signature + ".building-" + uuid.uuid4().hex)
        build.mkdir(parents=True)
        inputs = []
        for key in ("transcripts", "genome"):
            source = config[key]
            target = work / "downloads" / Path(urllib.parse.urlsplit(source["url"]).path).name
            download(source["url"], target, source["md5"])
            inputs.append(target)
        make_gentrome(inputs[0], inputs[1], build / "gentrome.fa", build / "decoys.txt")
        command = ["salmon", "index", "--transcripts", str(build / "gentrome.fa"),
                   "--decoys", str(build / "decoys.txt"), "--gencode", "--index", str(build / "index"),
                   "--kmerLen", str(config["index_k"]), "--threads", str(threads)]
        run_command(command, build / "index.log")
        if not (build / "index/versionInfo.json").is_file():
            raise ValueError("Salmon did not produce an index version file")
        finish_stage(build, final, signature, {"reference": config, "tools": versions, "command": command})


def make_gentrome(transcripts, genome, output, decoys):
    """Write transcript targets followed by genome decoys, checking sequence IDs."""
    seen = set()
    with Path(output).open("wt") as target, Path(decoys).open("wt") as names:
        for source, is_decoy in [(transcripts, False), (genome, True)]:
            count = 0
            with gzip.open(source, "rt") as handle:
                for line in handle:
                    if line.startswith(">"):
                        name = line[1:].split()[0]
                        identifier = name if is_decoy else name.split("|")[0]
                        if identifier in seen:
                            raise ValueError(f"Duplicate FASTA identifier: {identifier}")
                        seen.add(identifier)
                        count += 1
                        if is_decoy:
                            names.write(identifier + "\n")
                    elif count == 0 and line.strip():
                        raise ValueError("Sequence encountered before a FASTA header")
                    target.write(line.rstrip("\r\n") + "\n")
            if not count:
                raise ValueError("Reference FASTA is empty")


def fastq_record(handle):
    header = handle.readline()
    if not header:
        return None
    sequence, plus, quality = (handle.readline().rstrip("\r\n") for _ in range(3))
    if not header.startswith("@") or not plus.startswith("+") or not sequence:
        raise ValueError("Malformed or truncated FASTQ record")
    if len(sequence) != len(quality) or any(ord(char) < 33 or ord(char) > 126 for char in quality):
        raise ValueError("Invalid FASTQ quality string")
    if re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence) is None:
        raise ValueError("Invalid nucleotide sequence")
    name = re.sub(r"/[12]$", "", header[1:].split()[0])
    return name, len(sequence)


def check_pair(read1, read2):
    count, shortest, longest = 0, math.inf, 0
    with gzip.open(read1, "rt") as first, gzip.open(read2, "rt") as second:
        while True:
            one, two = fastq_record(first), fastq_record(second)
            if one is None and two is None:
                break
            if one is None or two is None or one[0] != two[0]:
                raise ValueError("FASTQ mates have different counts or read identifiers")
            count += 1
            shortest = min(shortest, one[1], two[1])
            longest = max(longest, one[1], two[1])
    if count == 0:
        raise ValueError("FASTQ pair is empty")
    return {"fragments": count, "shortest_read": shortest, "longest_read": longest}


def sample_signature(root, row, ref_signature):
    return fingerprint({"sample": row, "reference": ref_signature,
                        "qc": read_json(root / "config/qc.json"), "versions": VERSIONS,
                        "pipeline_sha256": digest(root / "workflow/rnaseq.py")})


def quantify(root, work, row, threads):
    versions = tool_versions()
    reference_id = reference_signature(root)
    reference = work / "reference" / reference_id
    if not checked_stage(reference, reference_id):
        raise ValueError("Build the reference before quantification")
    qc = read_json(root / "config/qc.json")
    signature = sample_signature(root, row, reference_id)
    final = work / "quant" / row["sample_id"] / signature
    with lock(work / "locks" / row["sample_id"]):
        if checked_stage(final, signature):
            print(f"Verified quantification checkpoint {row['sample_id']}")
            return
        raw_bytes = sum(int(row[f"read{mate}_bytes"]) for mate in (1, 2))
        require_space(work, raw_bytes + 5 * 2**30)
        reads = []
        for mate in (1, 2):
            path = work / "reads" / f"{row['run_accession']}_{mate}.fastq.gz"
            download(row[f"read{mate}_url"], path, row[f"read{mate}_md5"], int(row[f"read{mate}_bytes"]))
            reads.append(path)
        pair_stats = check_pair(*reads)
        build = final.with_name(signature + ".building-" + uuid.uuid4().hex)
        (build / "fastqc").mkdir(parents=True)
        write_json(build / "read_pair_qc.json", pair_stats)
        fastqc = ["fastqc", "--threads", str(min(threads, 2)), "--outdir", str(build / "fastqc"), *map(str, reads)]
        run_command(fastqc, build / "fastqc.log")
        for mate in (1, 2):
            for suffix in (".html", ".zip"):
                expected = build / "fastqc" / f"{row['run_accession']}_{mate}_fastqc{suffix}"
                if not expected.is_file():
                    raise ValueError(f"FastQC output missing: {expected.name}")
        command = ["salmon", "quant", "--index", str(reference / "index"), "--libType", qc["library_type"],
                   "--mates1", str(reads[0]), "--mates2", str(reads[1]), "--threads", str(threads),
                   "--validateMappings", "--numBootstraps", str(qc["bootstrap_replicates"]),
                   "--output", str(build / "salmon")]
        if qc["sequence_bias_correction"]:
            command.append("--seqBias")
        if qc["gc_bias_correction"]:
            command.append("--gcBias")
        run_command(command, build / "salmon.log")
        metrics = salmon_metrics(build / "salmon", pair_stats["fragments"], qc)
        write_json(build / "metrics.json", metrics)
        finish_stage(build, final, signature, {"sample": row, "tools": versions, "reference_id": reference_id,
                                             "commands": [fastqc, command], "qc": qc})


def salmon_metrics(directory, expected_fragments, qc):
    directory = Path(directory)
    meta = read_json(directory / "aux_info/meta_info.json")
    processed, mapped, rate = meta["num_processed"], meta["num_mapped"], meta["percent_mapped"]
    if processed != expected_fragments or not 0 <= mapped <= processed or processed <= 0:
        raise ValueError("Salmon fragment counts disagree with input FASTQ pairs")
    if not math.isfinite(rate) or abs(rate - mapped * 100 / processed) > 0.1:
        raise ValueError("Salmon mapping percentage disagrees with fragment counts")
    tpm_sum, count_sum, transcripts = 0.0, 0.0, set()
    with (directory / "quant.sf").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
            raise ValueError("Unexpected quant.sf columns")
        for row in reader:
            if row["Name"] in transcripts:
                raise ValueError("Duplicate transcript in quant.sf")
            transcripts.add(row["Name"])
            values = [float(row[key]) for key in ("Length", "EffectiveLength", "TPM", "NumReads")]
            if not all(math.isfinite(value) and value >= 0 for value in values) or values[0] <= 0:
                raise ValueError("Invalid numerical transcript estimate")
            tpm_sum += values[2]
            count_sum += values[3]
    if not transcripts or (mapped > 0 and abs(tpm_sum - 1e6) > 100):
        raise ValueError("Empty transcript table or TPM sum outside numerical tolerance")
    review = []
    if rate < qc["mapping_rate_review_below_percent"]:
        review.append("LOW_MAPPING_RATE")
    if processed < qc["minimum_fragments_review_below"]:
        review.append("LOW_FRAGMENT_COUNT")
    return {"num_processed": processed, "num_mapped": mapped, "percent_mapped": rate,
            "transcripts": len(transcripts), "tpm_sum": tpm_sum, "estimated_assigned_fragments": count_sum,
            "library_types": meta.get("library_types", []), "review_flags": review}


def fastqc_statuses(path):
    with zipfile.ZipFile(path) as archive:
        candidates = [name for name in archive.namelist() if name.endswith("/summary.txt")]
        if len(candidates) != 1:
            raise ValueError("FastQC archive must contain one summary.txt")
        lines = archive.read(candidates[0]).decode().splitlines()
    result = {}
    for line in lines:
        status, module, _ = line.split("\t")
        if status not in {"PASS", "WARN", "FAIL"} or module in result:
            raise ValueError("Invalid FastQC module status")
        result[module] = status
    if not result:
        raise ValueError("FastQC summary is empty")
    return result


def abundance_matrices(directories, sample_ids, output):
    """Merge a shared transcript index without dropping versioned identifiers."""
    from contextlib import ExitStack

    with ExitStack() as stack:
        readers = [csv.DictReader(stack.enter_context((path / "salmon/quant.sf").open()), delimiter="\t")
                   for path in directories]
        tpm = csv.writer(stack.enter_context((output / "transcript_tpm.tsv").open("w", newline="")), delimiter="\t")
        counts = csv.writer(stack.enter_context((output / "transcript_estimated_counts.tsv").open("w", newline="")), delimiter="\t")
        tpm.writerow(["transcript_id", *sample_ids])
        counts.writerow(["transcript_id", *sample_ids])
        for rows in zip_longest(*readers):
            if any(row is None for row in rows) or len({row["Name"] for row in rows}) != 1:
                raise ValueError("Transcript identifiers or ordering differ between samples")
            tpm.writerow([rows[0]["Name"], *(row["TPM"] for row in rows)])
            counts.writerow([rows[0]["Name"], *(row["NumReads"] for row in rows)])


def summarize(root, work, rows):
    versions = tool_versions()
    reference_id = reference_signature(root)
    signature = fingerprint({"samples": rows, "reference": reference_id, "qc": read_json(root / "config/qc.json")})
    final = work / "reports" / signature
    with lock(work / "locks/reports"):
        if checked_stage(final, signature):
            print(f"Verified report checkpoint {final}")
            return
        records, paths, quant_directories = [], [], []
        for row in rows:
            sample_id = sample_signature(root, row, reference_id)
            path = work / "quant" / row["sample_id"] / sample_id
            if not checked_stage(path, sample_id):
                raise ValueError(f"Missing completed sample: {row['sample_id']}")
            metrics = read_json(path / "metrics.json")
            record = {"sample_id": row["sample_id"], "run_accession": row["run_accession"], **metrics}
            record["review_flags"] = ";".join(metrics["review_flags"])
            record["library_types"] = ";".join(metrics["library_types"])
            for mate in (1, 2):
                statuses = fastqc_statuses(path / "fastqc" / f"{row['run_accession']}_{mate}_fastqc.zip")
                record[f"read{mate}_fastqc_warn"] = sum(status == "WARN" for status in statuses.values())
                record[f"read{mate}_fastqc_fail"] = sum(status == "FAIL" for status in statuses.values())
            records.append(record)
            quant_directories.append(path)
            paths.extend([str(path / "fastqc"), str(path / "salmon")])
        build = final.with_name(signature + ".building-" + uuid.uuid4().hex)
        build.mkdir(parents=True)
        with (build / "sample_qc.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
        abundance_matrices(quant_directories, [row["sample_id"] for row in rows], build)
        command = ["multiqc", "--force", "--module", "fastqc", "--module", "salmon", "--outdir", str(build),
                   "--filename", "multiqc_report.html", *paths]
        run_command(command, build / "multiqc.log")
        if not (build / "multiqc_report.html").is_file():
            raise ValueError("MultiQC did not produce an HTML report")
        finish_stage(build, final, signature, {"sample_sheet": rows, "tools": versions, "command": command})
        print(f"Completed cohort report: {final}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["preflight", "count", "reference", "quant", "summarize"])
    parser.add_argument("--samples", type=Path, default=ROOT / "config/pilot.tsv")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "work")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--task", type=int, help="One-based sample sheet row for a Slurm array")
    parser.add_argument("--require-tools", action="store_true")
    args = parser.parse_args()
    rows = read_samples(args.samples)
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.stage == "count":
        print(len(rows))
    elif args.stage == "preflight":
        raw = sum(int(row[f"read{mate}_bytes"]) for row in rows for mate in (1, 2))
        print(json.dumps({"runs": len(rows), "compressed_fastq_gib": round(raw / 2**30, 2),
                          "available_work_volume_gib": round(shutil.disk_usage(nearest_existing(args.work_dir)).free / 2**30, 2),
                          "reference_reserve_gib": 45, "suggested_work_volume_gib": math.ceil(raw / 2**30) + 60,
                          "tools": {tool: shutil.which(tool) for tool in VERSIONS}}, indent=2))
        if args.require_tools:
            tool_versions()
    elif args.stage == "reference":
        build_reference(ROOT, args.work_dir.resolve(), args.threads)
    elif args.stage == "quant":
        if args.task is None or not 1 <= args.task <= len(rows):
            parser.error("quant requires --task between 1 and the sample count")
        quantify(ROOT, args.work_dir.resolve(), rows[args.task - 1], args.threads)
    else:
        summarize(ROOT, args.work_dir.resolve(), rows)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
