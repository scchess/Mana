import os
import numpy as np
import pandas as pd


def parse(file):
    assert(os.path.exists(file))
    df = pd.read_csv(file, sep="\t")
    cols = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
    df = df.drop(df.columns[cols], axis=1)

    df["Match per base"] = df["matches"] / df["reads_all"]
    df["Mismatch per base"] = df["mismatches"] / df["reads_all"]
    df["Deletion per base"] = df["deletions"] / df["reads_all"]
    df["Insertion per base"] = df["insertions"] / df["reads_all"]

    match_mean = np.mean(df["Match per base"])
    match_min = np.min(df["Match per base"])
    match_max = np.max(df["Match per base"])

    mismatches_mean = np.mean(df["Mismatch per base"])
    mismatches_min = np.min(df["Mismatch per base"])
    mismatches_max = np.max(df["Mismatch per base"])

    deletions_mean = np.mean(df["Deletion per base"])
    deletions_min = np.min(df["Deletion per base"])
    deletions_max = np.max(df["Deletion per base"])

    insertions_mean = np.mean(df["Insertion per base"])
    insertions_min = np.min(df["Insertion per base"])
    insertions_max = np.max(df["Insertion per base"])

    total_mean = np.mean([mismatches_mean, deletions_mean, insertions_mean])
    total_min = np.min([mismatches_min, deletions_min, insertions_min])
    total_max = np.max([mismatches_max, deletions_max, insertions_max])

    mean_coverage = np.mean(df["reads_all"])
    min_coverage = np.min(df["reads_all"])
    max_coverage = np.max(df["reads_all"])

    return {"match_mean": match_mean,
            "match_min": match_min,
            "match_max": match_max,
            "mean_coverage": mean_coverage,
            "min_coverage": min_coverage,
            "max_coverage": max_coverage,
            "mismatches_mean": mismatches_mean,
            "mismatches_min": mismatches_min,
            "mismatches_max": mismatches_max,
            "deletions_mean": deletions_mean,
            "deletions_min": deletions_min,
            "deletions_max": deletions_max,
            "insertions_mean": insertions_mean,
            "insertions_min": insertions_min,
            "insertions_max": insertions_max,
            "total_mean": total_mean,
            "total_min": total_min,
            "total_max": total_max}
