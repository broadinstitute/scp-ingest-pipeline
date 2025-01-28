"""Ingest differential expression uploaded by study owner or editor

EXAMPLE:
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 ingest_differential_expression --annotation-name General_Celltype --annotation-type group --annotation-scope study --cluster-name cluster_umap_txt --study-accession SCPdev --ingest-differential-expression --differential-expression-file gs://fc-febd4c65-881d-497f-b101-01a7ec427e6a/author_de_test/lfc_qval_scanpy-like.csv --method wilcoxon
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import csv
import logging

from monitor import setup_logger, log_exception
from de import DifferentialExpression
from ingest_files import IngestFiles

sanitize_string = DifferentialExpression.sanitize_string


def sort_comparison(groups):
    """Naturally sort groups in a pairwise comparison; specially handle one-vs-rest

    https://en.wikipedia.org/wiki/Natural_sort_order

    :param groups (list<str>) A list of groups, e.g. ["B cells", "CSN1S1 macrophages"]
    """

    if any(i.isdigit() for i in groups):
        sorted_arr = sorted(
            groups, key=lambda x: int("".join([i for i in x if i.isdigit()]))
        )
        return sorted_arr
    elif "rest" == groups[1]:
        return groups
    elif "rest" == groups[0]:
        return [groups[1], groups[0]]
    else:
        return sorted(groups)


def canonicalize_name_and_order(data, header_refmap):
    """Ensure dataframe has column names and ordering as expected by DE UI

    Canonical order:
    gene,group,<comparison group,>size,significance,<other>
    """

    is_one_vs_rest_only = header_refmap["comparison_group"] == "None"
    if is_one_vs_rest_only:
        data = data.assign(comparison_group="rest")

    # Update column names in DF to expected header names
    rename_map = {}
    for ref_header in header_refmap:
        if ref_header in ["size", "significance"]:
            continue
        actual_header = header_refmap[ref_header]
        rename_map[actual_header] = ref_header
    data = data.rename(columns=rename_map)

    data = data.astype({"group": "string"})
    data['order'] = data.index.astype(int)

    # Update headers to expected order
    unsorted_headers = list(data.columns)
    size = header_refmap["size"]
    significance = header_refmap["significance"]
    sorted_headers = sort_headers(unsorted_headers, size, significance)
    canonical_long_df = data[sorted_headers]

    return canonical_long_df


def convert_long_to_wide(data):
    """Convert from canonical DE long format to SCP DE wide format

    (Long format is typical uploaded, but this module internally uses wide.)
    """
    metrics = list(data.columns[3:])
    data["combined"] = data["group"] + "--" + data["comparison_group"]
    frames = []
    # Example headers from a hypothetical author DE file:
    # genes  groups   comparison_groups    log2foldchanges  pvals_adj   qvals   means    cats dogs
    # metrics = ["log2foldchanges", "pvals_adj", "qvals", "mean", "cat", "dog"]
    for metric in metrics:
        wide_metric = pd.pivot(data, index="gene", columns="combined", values=metric)
        wide_metric = wide_metric.add_suffix(f"--{metric}")
        frames.append(wide_metric)

    result = pd.concat(frames, axis=1, join="inner")
    result.columns.name = " "
    result = result.reset_index()

    return result


def get_data_by_column(data):
    """
    outputs:
        - full list of column names
        - full list of genes
        - "rest" which is columns isolated from genes

    col_values outputs dict with the comparison characteristic and all the
    numeric data in an array, e.g. 'type_0'_'type_1'_qval: [0, 1, 2]
    """
    genes = data.iloc[:, 0].tolist()
    rest = data[data.columns[1:]]
    columns = list(rest.columns)

    # split col into two lists: pairwise, one_vs_rest
    split_values = {"one_vs_rest": [], "pairwise": []}
    for column in columns:
        if "rest" in column:
            split_values["one_vs_rest"].append(column)
        else:
            split_values["pairwise"].append(column)

    return columns, genes, rest, split_values


def generate_manifest(stem, clean_val, clean_val_p, qual):
    """Create manifest file of each comparison in the initial data
    if the comparison is with rest, rest is omitted and just the type is written
    """
    file_names_one_vs_rest = []
    for i in clean_val:
        clean_comparison = '--'.join([x for x in i if x != "rest" and x not in qual])
        if clean_comparison not in file_names_one_vs_rest:
            file_names_one_vs_rest.append(clean_comparison)

    file_names_pairwise = []
    for i in clean_val_p:
        values = [i[0], i[1]]
        if values not in file_names_pairwise:
            file_names_pairwise.append(values)

    with open(f"{stem}--manifest.tsv", "w", newline="") as f:
        tsv_output = csv.writer(f, delimiter="\t")
        if len(file_names_one_vs_rest) != 0:
            for value in range(len(file_names_one_vs_rest)):
                tsv_output.writerow([file_names_one_vs_rest[value]])

        if len(file_names_pairwise) != 0:
            for value in range(len(file_names_pairwise)):
                tsv_output.writerow(
                    [
                        file_names_pairwise[value][0],
                        file_names_pairwise[value][1],
                    ]
                )


def sort_all_group(all_group):
    """Filter and sort all_group so it can be later rearranged by a stride range

    Each inner_array in all_group can have raw, unsorted values like
    ['A--B--logfoldchanges', 'A--C--logfoldchanges', 'A--B--qval', 'A--C--qval', 'A--B--mean', 'A--C--mean']
    this sorts those like:
    ['A--B--logfoldchanges', 'A--B--mean', 'A--B--qval', 'A--C--logfoldchanges', 'A--C--mean', 'A--C--qval']

    This way, elements are sorted by 1st ***and 2nd group*** names,
    enabling grouped_comparison to rearrange with a simple stride length.
    """
    all_group_fin = []
    for inner_array in all_group:
        if inner_array == []:
            continue  # Filter out / skip empty arrays
        sorted_column_names = sorted(inner_array)
        all_group_fin.append(sorted_column_names)

    return all_group_fin


def sort_comparison_metrics(comparison_metrics, size, significance):
    """Ensure comparison_metrics has the order expected in the UI

    E.g., sort raw input:
    ['A--B--logfoldchanges', 'A--B--mean', 'A--B--qval']

    to output:
    ['A--B--logfoldchanges', 'A--B--qval', 'A--B--mean']

    (Gene, log2(fold change), q-value are the columns in the UI's DE table.)

    TODO: Generalize to output
    ['A--B--<size_metric>', 'A--B--<significance_metric>', <other metric identifiers>]

    Here `logfoldchanges` and `qval` are _particular_ size and significance
    metrics, but you can imagine how the author might provide `pvals_adj`
    instead of `qvals`, and we'd want to cleanly display `pvals_adj` (with a
    polished label, of course) in the UI.
    """

    # Sort alphabetically
    comparison_metrics = sorted(comparison_metrics)

    # Arrange significance in expected order (ultimately ranked 3rd)
    comparison_metrics = sorted(
        comparison_metrics, key=lambda x: x.split('--')[-1] == significance
    )

    # Arrange size in expected order (ultimately ranked 2nd)
    comparison_metrics = sorted(
        comparison_metrics, key=lambda x: x.split('--')[-1] == size
    )

    # Rank 1st with "gene", "group", then (if present) "comparison_group"
    comparison_metrics = sorted(
        comparison_metrics, key=lambda x: x.split('--')[-1] == "comparison_group"
    )
    comparison_metrics = sorted(
        comparison_metrics, key=lambda x: x.split('--')[-1] == "group"
    )
    comparison_metrics = sorted(
        comparison_metrics, key=lambda x: x.split('--')[-1] == "gene"
    )

    comparison_metrics.reverse()

    return comparison_metrics


def sort_headers(headers, size, significance):
    """Like `sort_comparison_metrics`, but for bare headers / metrics"""

    # Sort alphabetically
    headers = sorted(headers)

    # Rank significance 1st (ultimately ranked 5th)
    headers = sorted(headers, key=lambda x: x == significance)

    # Rank size 1st (ultimately ranked 4th)
    headers = sorted(headers, key=lambda x: x == size)

    # Rank 1st with "gene", "group", then (if present) "comparison_group"
    headers = sorted(headers, key=lambda x: x == "comparison_group")
    headers = sorted(headers, key=lambda x: x == "group")
    headers = sorted(headers, key=lambda x: x == "gene")

    headers.reverse()

    return headers


# note: my initial files had pval, qval, logfoldchanges.
# David's files have qval, mean, logfoldchanges.
# For the purposes of this validation I will be using his column values/formatting.


def validate_size_and_significance(metrics, size, significance, logger):
    """Locally log whether size and/or significance are detected among metrics

    TODO:
        - Log to Sentry / Mixpanel
    """
    has_size = any([metric.split('--')[-1] == size for metric in metrics])
    has_significance = any(
        [metric.split('--')[-1] == significance for metric in metrics]
    )

    in_headers = f"in headers: {metrics}"

    if not has_size or not has_significance:
        instruction = f'Column headers must include "{size}" and "{significance}".'
        if not has_size and not has_significance:
            msg = f"{instruction}  No such size or significance metrics found {in_headers}"
        elif not size:
            msg = f"{instruction}  No such size metric found {in_headers}"
        elif not has_significance:
            msg = f"{instruction}  No such significance metric found {in_headers}"
        logger.error(msg)
        raise ValueError(msg)
    elif has_size and has_significance:
        logger.info(
            f'Found size ("{size}") and significance ("{significance}") metrics {in_headers}'
        )


def get_groups_and_metrics(raw_column_names, size, significance, logger):
    """Cleans column names, splits them into compared groups and metrics

    A "metric" here is a statistical measure, like "logfoldchanges", "qval",
    "mean", or others provided by authors.

    A "group" either of the groups being compared, e.g. "B cells",
    "macrophages", "rest", or others provided by authors.  Note that "rest" is
    a special group that means "all other groups in this annotation".

    ---

    :param raw_column_names (list<str>) A list of raw column names, where
        each of the 3 items / parts of the raw column name elements is
        delimited by a double-hyphen.
        Element format: "<group>--<comparison group>--<metric>"
        Element example: "B cells--rest--logfoldchanges'"

    :return list<groups, split_headers, metrics>
    """
    split_headers = []

    for raw_column_name in raw_column_names:
        column_items = raw_column_name.split("--")
        split_header = []
        for item in column_items:
            item = item.replace(
                "'", ""
            )  # Remove quotes in e.g. 'type_0'--'type_1'--qval
            if (item != "") and (item != "_"):
                split_header.append(item.strip("_"))
        split_headers.append(split_header)

    groups = []
    metrics = []

    # isolate the groups and values in submitted file
    # expected format:
    # groups: ['group_0', 'group_1', 'group_2', 'group_3']
    # Example metrics: ['logfoldchanges', 'qval', 'mean']
    for split_header in split_headers:
        [group, comparison_group, metric] = split_header
        if group not in groups:
            groups.append(group)
        if comparison_group not in groups:
            groups.append(comparison_group)
        if metric not in metrics:
            metrics.append(metric)

    validate_size_and_significance(metrics, size, significance, logger)

    return groups, split_headers, metrics


def detect_seurat_findallmarkers(headers):
    """Reports whether headers likely derive from FindAllMarkers() in Seurat

    E.g.: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#finding-differentially-expressed-features-cluster-biomarkers

    These headers were observed in a real user-uploaded DE file.
    """
    findallmarkers_headers = [
        'p_val',
        'avg_log2FC',
        'pct.1',
        'pct.2',
        'p_val_adj',
        'cluster',
        'gene',
    ]
    is_seurat_findallmarkers = len(headers) == len(findallmarkers_headers) and all(
        headers == findallmarkers_headers
    )
    return is_seurat_findallmarkers


def order_not_significant(array_1, array_2):
    correlation, pval = spearmanr(array_1, array_2)
    if correlation > 0.95:
        return False
    else:
        return True


def organize_results(df):
    # processing turned values into strings, convert to numeric for sorting
    df["order"] = df["order"].astype(float)
    df["order"] = df["order"].astype(int)
    # sort dataframe by input row order
    df = df.sort_values(by="order")
    df = df.set_index('order')
    # maintain unnamed index column in DE results file
    df.index.name = None
    # processing ensures the significance metric is the 3rd column of the df
    df[df.columns[2]] = df[df.columns[2]].astype(float)
    input_order = df[df.columns[2]].to_numpy()
    sig_sorted = df[df.columns[2]].sort_values()
    sig_array = sig_sorted.to_numpy()

    if order_not_significant(sig_array, input_order):
        # sort dataframe by significance metric
        df = df.sort_values(by=df.columns[2])
        return df
    else:
        # leave dataframe sorted by input row order
        return df


class AuthorDifferentialExpression:
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")
    author_de_logger = setup_logger(
        __name__ + ".author_de_logger",
        "author_de_log.txt",
        level=logging.INFO,
        format="support_configs",
    )

    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    def __init__(
        self,
        cluster_name,
        annotation_name,
        study_accession,
        annotation_scope,
        method,
        differential_expression_file,
        header_refmap,
    ):
        """
        :param cluster_name (string) Name of cluster, e.g. "All Cells UMAP"
        :param annotation_name (string) Name of annotation, e.g. "General_Celltype"
        :param study_accession (string) Study accession, e.g. "SCP123"
        :param annotation_scope (string) Scope of annotation; one of "cluster" or "study"
        :param method (string) Statistical method used for DE, e.g. "wilcoxon"
        :param differential_expression_file (string) Path to author DE file
        :param header_refmap (dict): Map of canonical DE headers to actual headers in DE file
        """
        # AuthorDifferentialExpression.de_logger.info(
        #    "Initializing DifferentialExpression instance"
        # )
        self.cluster_name = sanitize_string(cluster_name)
        self.annotation = sanitize_string(annotation_name)
        self.accession = study_accession
        self.annot_scope = annotation_scope
        self.method = method
        self.header_refmap = header_refmap
        self.size_metric = header_refmap["size"]
        self.significance_metric = header_refmap["significance"]
        self.stem = f"{self.cluster_name}--{self.annotation}"

        author_de_file_gcs_url = differential_expression_file
        allowed_file_types = AuthorDifferentialExpression.ALLOWED_FILE_TYPES
        raw_de_file_obj = IngestFiles(author_de_file_gcs_url, allowed_file_types)
        author_file_handle, local_file_path = IngestFiles.resolve_path(
            raw_de_file_obj, author_de_file_gcs_url
        )
        self.author_de_file = local_file_path

    def execute(self):
        logger = AuthorDifferentialExpression.dev_logger

        clean_val = []
        clean_val_p = []
        metrics = []
        file_path = self.author_de_file

        # sep=None invokes detecting separator via csv.Sniffer, per
        # https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
        data = pd.read_csv(file_path, sep=None)

        headers = data.columns

        is_wide_format = '--' in headers[1]
        if is_wide_format:
            data_by_col = get_data_by_column(data)
        else:
            data = canonicalize_name_and_order(data, self.header_refmap)
            wide_df = convert_long_to_wide(data)
            data_by_col = get_data_by_column(wide_df)

        col, genes, rest, split_values = data_by_col
        pairwise = split_values["pairwise"]
        one_vs_rest = split_values["one_vs_rest"]

        if len(one_vs_rest) != 0:
            groups, clean_val, metrics = get_groups_and_metrics(
                one_vs_rest, self.size_metric, self.significance_metric, logger
            )
            self.generate_result_files(
                one_vs_rest, genes, rest, groups, clean_val, metrics
            )

        if len(pairwise) != 0:
            groups_p, clean_val_p, metrics = get_groups_and_metrics(
                pairwise, self.size_metric, self.significance_metric, logger
            )
            self.generate_result_files(
                pairwise, genes, rest, groups_p, clean_val_p, metrics
            )
        generate_manifest(self.stem, clean_val, clean_val_p, metrics)

        print("Author DE transformation succeeded")

    def generate_result_files(self, col, genes, rest, groups, clean_val, metrics):
        """
        Create an individual DE result file for each comparison, pairwise or rest,
        with all the metrics being used (e.g. logfoldchanges, qval, mean)

        For the desired format, e.g. if we have:
            type_0--type_1--logfoldchanges
            type_0--type_1--qval
            type_0--type_1--mean

        Then final format should have type 0 type 1 in the title, and genes, logfoldchanges, qval, and mean as columns
        """

        comparisons_dict = {}
        all_group = []
        for i in range(len(groups)):
            curr_group = groups[i]
            type_group = []
            for j in range(len(clean_val)):
                group = clean_val[j][0]
                comparison_group = clean_val[j][1]
                comparison = f"{group}--{comparison_group}"
                comparisons_dict[comparison] = []
                real_title = col[j]
                if curr_group == group:
                    type_group.append(real_title)

            all_group.append(type_group)

        # An array of arrays of comparison-metrics, where each inner array has
        # elements with the same 1st-group name
        all_group_fin = sort_all_group(all_group)

        grouped_comparison_metrics = []
        num_metrics = len(metrics)  # Stride length
        for i in all_group_fin:
            for j in range(0, len(i), num_metrics):
                x = j
                comparison_metrics = i[x : x + num_metrics]
                sorted_comparison_metrics = sort_comparison_metrics(
                    comparison_metrics, self.size_metric, self.significance_metric
                )
                grouped_comparison_metrics.append(sorted_comparison_metrics)

        for comparison in comparisons_dict:
            for comparison_metrics in grouped_comparison_metrics:
                for comparison_metric in comparison_metrics:
                    if comparison in comparison_metric:
                        comparisons_dict[comparison].append(comparison_metric)

        # Now we have all the columns grouped in lists by comparison name, with logfoldchanges, qval, mean
        # have to pair with corresponding gene for that row
        # dictionary format:
        # comparison name: [[gene, logfoldchanges, qval, mean], [gene, logfoldchanges, qval, mean], ...]
        comparisons = comparisons_dict.keys()
        rows_by_comparison = dict.fromkeys(comparisons, [])

        for comparison_metrics in grouped_comparison_metrics:
            metric_values = []
            for comparison_metric in comparison_metrics:

                # E.g. "B cells--CSN1S1 macrophages"
                comparison = '--'.join(comparison_metric.split('--')[0:2])

                # Numerical values for metric in this comparison
                values = rest[comparison_metric].tolist()
                metric_values.append(values)

            rows = metric_values
            rows.insert(0, genes)
            rows_by_comparison[comparison] = rows

        headers = sort_headers(metrics, self.size_metric, self.significance_metric)
        headers.insert(0, "gene")

        for comparison_metric in rows_by_comparison:

            if "rest" in comparison_metric:
                comparison = comparison_metric.split("--")[0]

            else:
                group = comparison_metric.split("--")[0]
                comparison_group = comparison_metric.split("--")[1]
                sorted_list = sort_comparison([group, comparison_group])
                comparison = f'{sorted_list[0]}--{sorted_list[1]}'

            clean_comparison_metric = '--'.join(
                [sanitize_string(group) for group in comparison.split('--')]
            )

            tsv_name = f'{self.stem}--{clean_comparison_metric}--{self.annot_scope}--{self.method}.tsv'

            rows = rows_by_comparison[comparison_metric]

            arr = np.array(rows)

            t_arr = arr.transpose()

            if len(t_arr) == 0:
                print(
                    f"No data to output for TSV, skip preparation to write {tsv_name}"
                )
                continue

            # Drop rows that are all "nan", as seen sometimes in Seurat FindAllMarkers()
            t_arr = t_arr[~(t_arr == 'nan').any(axis=1)]

            inner_df = pd.DataFrame(data=t_arr, columns=headers)
            inner_df = organize_results(inner_df)
            inner_df.to_csv(tsv_name, sep='\t')

            print(f"Wrote TSV: {tsv_name}")
