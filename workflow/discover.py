"""Join GEO organoid annotations to ENA run files by experiment and BioSample."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.request

from rnaseq import FIELDS, ROOT, read_samples

ENA_FIELDS = "run_accession,experiment_accession,sample_accession,study_accession,library_layout,library_strategy,scientific_name,fastq_ftp,fastq_bytes,fastq_md5"
ENA_URL = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP158722&result=read_run&fields={ENA_FIELDS}&format=tsv"
SOFT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119004/soft/GSE119004_family.soft.gz"
PILOT_GEO = {"GSM3356477", "GSM4063806", "GSM4063812"}


def organoids(soft):
    records, current = [], None
    for line in soft.splitlines():
        if line.startswith("^"):
            if current is not None:
                records.append(current)
            current = {"geo_accession": line.split(" = ", 1)[1]} if line.startswith("^SAMPLE = ") else None
        elif current is not None and line.startswith("!Sample_"):
            key, value = line.split(" = ", 1)
            current.setdefault(key, []).append(value)
    if current is not None:
        records.append(current)
    selected = []
    for record in records:
        characteristics = record.get("!Sample_characteristics_ch1", [])
        if "sample type: organoid" not in [value.lower().strip() for value in characteristics]:
            continue
        relations = "\n".join(record.get("!Sample_relation", []))
        experiments = re.findall(r"term=(SRX\d+)", relations)
        biosamples = re.findall(r"biosample/(SAMN\d+)", relations)
        if len(experiments) != 1 or len(biosamples) != 1:
            raise ValueError("Organoid needs exactly one experiment and BioSample relation")
        selected.append({"geo_accession": record["geo_accession"], "title": record["!Sample_title"][0],
                         "experiment_accession": experiments[0], "biosample": biosamples[0]})
    if not selected or len({item["geo_accession"] for item in selected}) != len(selected):
        raise ValueError("Missing or duplicate organoid records")
    return selected


def join_runs(selected, ena_text):
    ena = list(csv.DictReader(io.StringIO(ena_text), delimiter="\t"))
    rows = []
    for sample in selected:
        matches = [row for row in ena if row["experiment_accession"] == sample["experiment_accession"]]
        if len(matches) != 1:
            raise ValueError("This cohort expects one ENA run per GEO organoid")
        run = matches[0]
        if (run["sample_accession"] != sample["biosample"] or run["study_accession"] != "PRJNA487654"
                or run["library_layout"] != "PAIRED" or run["library_strategy"] != "RNA-Seq"
                or run["scientific_name"] != "Homo sapiens"):
            raise ValueError("ENA study, sample, assay or layout disagrees with cohort scope")
        columns = [run[field].split(";") for field in ("fastq_ftp", "fastq_md5", "fastq_bytes")]
        if any(len(values) != 2 for values in columns):
            raise ValueError("Expected exactly two FASTQ mates")
        row = {"sample_id": sample["geo_accession"], "geo_accession": sample["geo_accession"],
               "experiment_accession": sample["experiment_accession"], "run_accession": run["run_accession"]}
        for mate in (1, 2):
            row[f"read{mate}_url"] = "https://" + columns[0][mate - 1]
            row[f"read{mate}_md5"] = columns[1][mate - 1]
            row[f"read{mate}_bytes"] = columns[2][mate - 1]
        rows.append(row)
    return sorted(rows, key=lambda row: row["geo_accession"])


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--soft", type=Path, help="Reuse a local compressed GEO SOFT snapshot")
    parser.add_argument("--ena", type=Path, help="Reuse a local ENA file report")
    parser.add_argument("--out", type=Path, default=ROOT)
    args = parser.parse_args()
    if args.soft:
        source = args.soft.read_bytes()
    else:
        with urllib.request.urlopen(SOFT_URL, timeout=60) as response:
            source = response.read()
    if args.ena:
        ena = args.ena.read_bytes()
    else:
        with urllib.request.urlopen(ENA_URL, timeout=60) as response:
            ena = response.read()
    selected = organoids(gzip.decompress(source).decode())
    rows = join_runs(selected, ena.decode())
    pilot = [row for row in rows if row["geo_accession"] in PILOT_GEO]
    if len(rows) != 9 or len(pilot) != 3:
        raise ValueError("Cohort membership changed; review before replacing the pinned sample sheets")
    (args.out / "config").mkdir(parents=True, exist_ok=True)
    metadata = args.out / "data/metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    for name, cohort in [("cohort.tsv", rows), ("pilot.tsv", pilot)]:
        write_tsv(args.out / "config" / name, cohort, FIELDS)
        read_samples(args.out / "config" / name)
    write_tsv(metadata / "organoid_annotations.tsv", selected, list(selected[0]))
    (metadata / "ena_runs.tsv").write_bytes(ena)
    provenance = {"retrieved_utc": datetime.now(timezone.utc).isoformat(), "geo_source": SOFT_URL,
                  "ena_source": ENA_URL, "soft_sha256": hashlib.sha256(source).hexdigest(),
                  "ena_sha256": hashlib.sha256(ena).hexdigest(), "organoid_runs": len(rows),
                  "pilot_geo_accessions": sorted(PILOT_GEO), "selection": "Documented GEO sample type, joined by SRX and confirmed BioSample"}
    (metadata / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Validated {len(rows)} organoid runs; selected {len(pilot)} pilot runs")


if __name__ == "__main__":
    main()
