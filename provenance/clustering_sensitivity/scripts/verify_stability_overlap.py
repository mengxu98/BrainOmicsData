"""Recompute full-cohort partition metrics from retained overlap counts."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import csv
import hashlib
import json
import math
import sys


def read_table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify(directory):
    parameters, = read_table(directory / "parameters.tsv")
    summary, = read_table(directory / "summary.tsv")
    completion = json.loads((directory / "COMPLETE.json").read_text())
    overlaps = read_table(directory / "overlap.tsv")
    cluster_rows = read_table(directory / "cluster_stability.tsv")
    expected = {
        row["Cluster"]: int(row["Cells"])
        for row in read_table(
            Path(__file__).resolve().parents[1]
            / "annotations"
            / "final_cluster_annotation.tsv"
        )
    }

    baseline, comparison, pairs = Counter(), Counter(), {}
    for row in overlaps:
        pair = row["Baseline"], row["Perturbed"]
        count = int(row["Freq"])
        assert count > 0 and pair not in pairs
        pairs[pair] = count
        baseline[pair[0]] += count
        comparison[pair[1]] += count

    assert dict(baseline) == expected
    cells = sum(baseline.values())
    assert cells == 2602031 == int(summary["Cells"]) == int(parameters["Cells"])
    assert int(summary["Clusters"]) == len(comparison)
    assert int(parameters["Task"]) == int(summary["Task"])
    assert parameters["Group_Singletons"] == "TRUE"
    assert parameters["Subsampling"] == "FALSE"
    assert int(parameters["Seed"]) == int(summary["Seed"])
    assert float(parameters["Resolution"]) == float(summary["Resolution"])
    assert float(summary["Seconds"]) > 0
    assert completion["state"] == "COMPLETE"
    assert completion["task"] == int(summary["Task"])
    assert completion["cells"] == int(summary["Cells"])

    choose2 = lambda value: value * (value - 1) // 2
    row_pairs = sum(map(choose2, baseline.values()))
    column_pairs = sum(map(choose2, comparison.values()))
    observed = sum(map(choose2, pairs.values()))
    expected_pairs = Fraction(row_pairs * column_pairs, choose2(cells))
    ari = float(
        (observed - expected_pairs)
        / (Fraction(row_pairs + column_pairs, 2) - expected_pairs)
    )
    assert math.isclose(ari, float(summary["ARI"]), abs_tol=1e-12)
    assert math.isclose(completion["ARI"], ari, abs_tol=5e-5)
    assert len(cluster_rows) == 75
    assert {row["Baseline"] for row in cluster_rows} == set(baseline)
    for row in cluster_rows:
        key = row["Baseline"]
        fragments = {
            right: count
            for (left, right), count in pairs.items()
            if left == key
        }
        assert int(row["Cells"]) == baseline[key]
        assert int(row["Fragments"]) == len(fragments)
        largest = max(fragments.values()) / baseline[key]
        jaccard = max(
            count / (baseline[key] + comparison[right] - count)
            for right, count in fragments.items()
        )
        assert math.isclose(
            largest, float(row["Largest_Fragment_Fraction"]), abs_tol=1e-12
        )
        assert math.isclose(
            jaccard, float(row["Best_Jaccard"]), abs_tol=1e-12
        )

    retained_files = [
        "parameters.tsv", "summary.tsv", "COMPLETE.json", "overlap.tsv",
        "cluster_stability.tsv"
    ]
    return {
        "state": "verified",
        "task": int(parameters["Task"]),
        "cells": cells,
        "clusters": len(comparison),
        "ARI": ari,
        "graph_md5": parameters["Graph_MD5"],
        "reference_assignment_md5": parameters["Reference_Assignment_MD5"],
        "input_sha256": {
            name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
            for name in retained_files
        },
    }


if __name__ == "__main__":
    print(json.dumps(verify(Path(sys.argv[1])), indent=2))
