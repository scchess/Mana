import os
import numpy as np
import pandas as pd


def parse(file):
    assert(os.path.exists(file))
    df = pd.read_csv(file, sep="\t")
    cols = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
    df = df.drop(df.columns[cols], axis=1)

    match_mean = np.mean(df["matches"])
    match_min = np.min(df["matches"])
    match_max = np.max(df["matches"])

    mismatches_mean = np.mean(df["mismatches"])
    mismatches_min = np.min(df["mismatches"])
    mismatches_max = np.max(df["mismatches"])

    deletions_mean = np.mean(df["deletions"])
    deletions_min = np.min(df["deletions"])
    deletions_max = np.max(df["deletions"])

    insertions_mean = np.mean(df["insertions"])
    insertions_min = np.min(df["insertions"])
    insertions_max = np.max(df["insertions"])

    total_mean = np.mean([mismatches_mean, deletions_mean, insertions_mean])
    total_min = np.min([mismatches_min, deletions_min, insertions_min])
    total_max = np.min([mismatches_max, deletions_max, insertions_max])

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
