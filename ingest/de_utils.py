"""Utilities for differential expression
"""

import re


def get_size(metric):
    # Scanpy: logfoldchanges; Seurat: avg_log2FC
    SIZE_RE = r"(logfoldchange|log2FC|log2foldchange)"
    match = re.search(SIZE_RE, metric, re.IGNORECASE)
    size = match.group() if match else None
    return size


### Significance parsers ###
### TODO (SCP-): Custom tooltips for custom metrics
def get_pval_adj(metric):
    # Scanpy: pvals_adj; Seurat: p_val_adj
    ADJUSTED_P_VALUE_RE = r"(pvals_adj|p_val_adj)"
    match = re.search(ADJUSTED_P_VALUE_RE, metric, re.IGNORECASE)
    pval_adj = match.group() if match else None
    return pval_adj


def get_pval(metric):
    # Scanpy: pvals; Seurat: p_val
    P_VALUE_RE = r"(pval|p_val)"
    match = re.search(P_VALUE_RE, metric, re.IGNORECASE)
    pval = match.group() if match else None
    pval_adj = get_pval_adj(metric)
    if not pval_adj:
        return pval


def get_qval(metric):
    # Scanpy: qvals; Seurat: q_val (?)
    Q_VALUE_RE = r"(qval|q_val)"
    match = re.search(Q_VALUE_RE, metric, re.IGNORECASE)
    qval = match.group() if match else None
    return qval


def get_significance(metric):
    pval_adj = get_pval_adj(metric)
    if pval_adj:
        return pval_adj
    else:
        pval = get_pval(metric)
        if pval:
            return pval
        else:
            qval = get_qval(metric)
            if qval:
                return qval
            else:
                return None
### End significance parsers ###


def get_size_and_significance(metrics):
    sig_list = list(filter(get_size, metrics))
    significance = sig_list[0] if len(sig_list) > 0 else None
    size_list = list(filter(get_significance, metrics))
    size = size_list[0] if len(size_list) > 0 else None

    return {size, significance}
