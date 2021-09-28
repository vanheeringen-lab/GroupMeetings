# 1) Y-axis normalization

import pandas as pd
import numpy as np
from qnorm import quantile_normalize


def tpm_normalization(
        tpms: pd.DataFrame,
        column_order: list,
        minimum_value: int = None,
) -> pd.DataFrame:
    """filter and order a tpm table, then quantile normalize and log transform"""
    bc = tpms[column_order]                       # filter & order samples
    if minimum_value:
        b4 = bc.shape[0]
        bc = bc[bc.max(axis=1) >= minimum_value]  # filter genes
        aft = b4 - bc.shape[0]
        print(f"Genes with TPM below {minimum_value}: {aft} of {b4} ({round(100*aft/b4,0)}%)")
    bc = quantile_normalize(bc, axis=1)           # normalize
    bc = np.log2(bc+1)                            # transform
    return bc


# preprocessing
template_samples = pd.read_csv("/path/to/samples_file", sep="\t", index_col=0)
sample_order, time2samples = get_sample_info(template_samples)  # function not in this code review
# head(sample_order, 4): ['1-cell-1', '1-cell-2', '2-cell-1', '2-cell-2']

tpms = pd.read_csv("/path/to/tpm_file", sep="\t", index_col=0)
template_tpms = tpm_normalization(tpms, sample_order, minimum_value=5)


# 2) X-axis normalization


def _str_inf_ts(timepoints, n):
    t_extended = []
    for t in timepoints:
        for dec in range(n):
            if dec == 0:
                t_extended.append(f"{t}")
            else:
                t_extended.append(f"{t}+{dec}/{n}")
    return t_extended


def _flt_inf_ts(timepoints, n):
    t_extended = []
    timepoints = [round(float(e), 2) for e in timepoints]
    for t in range(len(timepoints)):
        curr_timepoint = timepoints[t]
        diff_time = 0
        if t < len(timepoints)-1:
            next_timepoint = timepoints[t+1]
            diff_time = next_timepoint-curr_timepoint
        for dec in range(n):
            t_extended.append(round(curr_timepoint + diff_time*(dec/n), 2))
    return t_extended


def inference_timeseries(timepoints: list, n: int = 10) -> list:
    """
    Extends the given list of timepoints to n points between for each original point, up to the final point.
    Returns the extended list of named timepoints with added time.
    e.g. func([1, 2, 3], n=3) -> [1, 1.33, 1.66, 2, 2.33, 2.66, 3]
    e.g. func(["a","b","c"], n=3) -> ["a", "a+1/3", "a+2/3", "b", "b+1/3", "b+2/3", "c"]
    """
    if all_numeric(timepoints):  # function not in this code review
        t_extended = _flt_inf_ts(timepoints, n)
    else:
        t_extended = _str_inf_ts(timepoints, n)

    # remove the points after the last real timepoint
    total_timepoints = len(timepoints) * n - n + 1
    t_extended = t_extended[:total_timepoints]

    return t_extended


extended_timepoints = inference_timeseries(list(time2samples), n=10)
# print(len(extended_timepoints)): 171
# print(extended_timepoints[0:4]): [0.0, 4.5, 9.0, 13.5]


# 3) Distance metric

from scipy.spatial.distance import cdist


def get_cost_matrix(template_tpms, query_tpms, metric='euclidean'):
    template = template_tpms.to_numpy(dtype=np.float64)
    query = query_tpms.to_numpy(dtype=np.float64)
    cost_matrix = cdist(template.T, query.T, metric=metric).T  # pairwise distance matrix
    return cost_matrix

cost_matrix = get_cost_matrix(template_tpms, query_tpms, metric)
# metrics = ["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "dice", "euclidean", "hamming", "jaccard", "jensenshannon", "kulsinski", "mahalanobis", "matching", "minkowski", "rogerstanimoto", "russellrao", "seuclidean", "sokalmichener", "sokalsneath", "sqeuclidean", "wminkowski", "yule"]
