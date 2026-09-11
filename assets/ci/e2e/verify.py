#!/usr/bin/env python3
"""Verify deterministic biological invariants of an xasm E2E profile run."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


# lives at assets/ci/e2e/verify.py, three levels below the repository root
ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "assets/test/test_data/e2e"
GOLDEN = json.loads((Path(__file__).with_name("golden.json")).read_text())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def verify_fixture() -> None:
    manifest = json.loads((FIXTURE / "manifest.json").read_text())
    for relative, expected in manifest["sha256"].items():
        path = FIXTURE / relative
        require(path.is_file(), f"fixture file missing: {relative}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        require(actual == expected, f"fixture checksum changed: {relative}")

    for sample, designed in manifest["reads"].items():
        observed: Counter[str] = Counter()
        with gzip.open(FIXTURE / "reads" / f"{sample}_1.fastq.gz", "rt") as handle:
            for line_number, line in enumerate(handle):
                if line_number % 4 == 0:
                    observed[line.rstrip().split(":", 2)[1]] += 1
        for category, expected in designed.items():
            name = category.removesuffix("_pairs")
            require(observed[name] == expected, f"{sample} {name}: {observed[name]} != {expected}")


def line_count(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for _ in handle)


def chromosome_counts(path: Path, gtf: bool = False) -> tuple[Counter[str], list[str]]:
    counts: Counter[str] = Counter()
    transcript_ids: list[str] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if gtf and (len(fields) < 9 or fields[2] != "transcript"):
                continue
            counts[fields[0]] += 1
            if gtf:
                match = re.search(r'transcript_id "([^"]+)"', fields[8])
                require(match is not None, f"transcript lacks transcript_id: {line.rstrip()}")
                transcript_ids.append(match.group(1))
            else:
                require(len(fields) >= 4, f"BED record has fewer than four fields: {line.rstrip()}")
                transcript_ids.append(fields[3])
    return counts, transcript_ids


def latest_trace(output: Path) -> Path:
    traces = list((output / "pipeline_info").glob("execution_trace_*.txt"))
    require(bool(traces), f"no execution trace under {output}")
    return max(traces, key=lambda path: path.stat().st_mtime_ns)


def completed_tasks(output: Path) -> Counter[str]:
    task_counts: Counter[str] = Counter()
    with latest_trace(output).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            require(row["status"] == "COMPLETED", f"task did not complete: {row['name']} ({row['status']})")
            process = row["name"].split(" (")[0].rsplit(":", 1)[-1]
            task_counts[process] += 1
    return task_counts


def verify_profile(profile: str) -> None:
    expected = GOLDEN["profiles"][profile]
    output = ROOT / expected["output"]
    require(output.is_dir(), f"profile output missing: {output}")

    samples: dict[str, list[str]] = {}
    with (output / "samplesheets/samplesheet.csv").open(newline="") as handle:
        for row in csv.reader(handle):
            require(len(row) == 10, f"unexpected samplesheet row: {row}")
            require(
                row[9] in {"forward", "reverse", "unstranded"},
                f"unexpected strandedness: {row}",
            )
            samples[row[0]] = row

    task_counts = completed_tasks(output)

    if expected.get("smoke"):
        # No golden counts exist for this profile; assert the CBQ path ran end to end
        # and that the STAR/FASTP path it replaces did not.
        require(set(samples) == set(GOLDEN["common_samples"]), f"unexpected samples: {sorted(samples)}")
        for process in expected["tasks"]:
            require(task_counts[process] >= 1, f"{process} did not run")
        for process in expected.get("forbidden_tasks", []):
            require(task_counts[process] == 0, f"task that should be absent ran: {process}")
        print(f"verified {profile}: smoke OK")
        return

    require(set(samples) == set(GOLDEN["common_samples"]), f"unexpected samples: {sorted(samples)}")
    for sample, common in GOLDEN["common_samples"].items():
        row = samples[sample]
        require(int(row[3]) == common["reads_after_trim"], f"{sample}: trimmed-read count changed")
        require(abs(float(row[4]) - common["trim_percent"]) < 1e-7, f"{sample}: trim percent changed")
        require(abs(float(row[5]) - common["deacon_retained_percent"]) < 1e-4, f"{sample}: Deacon count changed")
        require(abs(float(row[6]) - common["mapped_percent"]) < 1e-4, f"{sample}: STAR mapping changed")
        require(int(row[8]) == expected["assembled"][sample], f"{sample}: assembly count changed")

    gtf_counts, transcript_ids = chromosome_counts(output / expected["metassembly"], gtf=True)
    require(dict(gtf_counts) == expected["gtf_transcripts"], f"metassembly counts changed: {dict(gtf_counts)}")
    require(len(transcript_ids) == len(set(transcript_ids)), "metassembly transcript IDs are not globally unique")

    final_total = None
    if expected.get("final_glob"):
        final_files = list((output / "10_final").glob(expected["final_glob"]))
        require(len(final_files) == 1, f"expected one final BED, found {len(final_files)}")
        final_counts, final_ids = chromosome_counts(final_files[0])
        require(dict(final_counts) == expected["final_records"], f"final predictions changed: {dict(final_counts)}")
        require(len(final_ids) == len(set(final_ids)), "final transcript IDs are not unique")
        with final_files[0].open() as handle:
            coordinates = [
                (fields[0], int(fields[1]), int(fields[2]))
                for line in handle if line.strip()
                for fields in [line.rstrip().split("\t")]
            ]
        require(coordinates == sorted(coordinates), "final BED is not coordinate sorted")
        for relative, count in expected["polish_line_counts"].items():
            # '*' marks filenames that are order-dependent (xorf merged outputs are
            # named after whichever seleno group wins the submodule collect); require
            # exactly one match so the assertion stays strict.
            if "*" in relative:
                matches = list(output.glob(relative))
                require(len(matches) == 1, f"polishing checkpoint glob ambiguous: {relative}: {matches}")
                path = matches[0]
            else:
                path = output / relative
            require(path.is_file(), f"polishing checkpoint missing: {relative}")
            require(line_count(path) == count, f"polishing checkpoint changed: {relative}")
        for relative in expected.get("absent_paths", []):
            require(not (output / relative).exists(), f"unexpected polishing output: {relative}")
        final_total = sum(final_counts.values())

    task_counts = completed_tasks(output)
    for process, count in expected["tasks"].items():
        require(task_counts[process] == count, f"{process} fan-out changed: {task_counts[process]} != {count}")
    for process in expected.get("forbidden_tasks", []):
        require(task_counts[process] == 0, f"disabled two-pass task ran unexpectedly: {process}")

    if final_total is None:
        print(f"verified {profile}: {len(transcript_ids)} metassembly transcript(s)")
    else:
        print(f"verified {profile}: {len(transcript_ids)} metassembly and {final_total} final transcript(s)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("profile", choices=tuple(GOLDEN["profiles"]) + ("all",))
    args = parser.parse_args()
    verify_fixture()
    profiles = GOLDEN["profiles"] if args.profile == "all" else (args.profile,)
    for profile in profiles:
        verify_profile(profile)


if __name__ == "__main__":
    main()
