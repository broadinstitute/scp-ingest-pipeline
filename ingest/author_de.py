"""Ingest differential expression uploaded by study owner or editor

EXAMPLE:
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 ingest_differential_expression --annotation-name General_Celltype --annotation-type group --annotation-scope study --cluster-name cluster_umap_txt --study-accession SCPdev --ingest-differential-expression --differential-expression-file gs://fc-febd4c65-881d-497f-b101-01a7ec427e6a/author_de_test/lfc_qval_scanpy-like.csv --method wilcoxon
"""

import pandas as pd
import numpy as np
import csv
import logging

from monitor import setup_logger, log_exception
from de import DifferentialExpression
from ingest_files import IngestFiles

sanitize_string = DifferentialExpression.sanitize_string


def check_group(names_dict, name):
    """Helper function to return the comparison based on the comparison with the value

    Takes in dictionary that stores all names, and given name

    E.g. input type_2_type_3_mean returns type_2_type_3
    """
    for i in names_dict:
        if name in names_dict[i]:
            return i


def sort_comparison(groups):
    """Naturally sort groups in a pairwise comparison; specially handle one-vs-rest

    https://en.wikipedia.org/wiki/Natural_sort_order

    :param groups (list<str>) A list of groups, e.g. ["B cells", "CSN1S1 macrophages"]
    """

    if any(i.isdigit() for i in groups):
        sorted_arr = sorted(groups, key=lambda x: int("".join([i for i in x if i.isdigit()])))
        return sorted_arr
    elif "rest" == groups[1]:
        return groups
    elif "rest" == groups[0]:
        return [groups[1], groups[0]]
    else:
        return sorted(groups)


def convert_seurat_findallmarkers_to_wide(data):
    """Convert from Seurat FindAllMarkers() format to SCP DE wide format
    """
    print('0')
    data = data.rename(columns={"cluster": "group", "gene": "genes"})
    data = data.astype({"group": "string"})
    print('1')
    data = data.assign(comparison_group="rest")
    expected_header_order = [
        "genes", "group", "comparison_group", "avg_log2FC", "p_val_adj", "p_val", "pct.1", "pct.2"
    ]
    data = data[expected_header_order]
    print('reindexed data')
    print(data)
    # print("data.dtypes")
    # print(data.dtypes)
    wide_df = convert_long_to_wide(data)
    print("wide_df data.dtypes")
    print(data.dtypes)
    print('3')
    return wide_df


def convert_long_to_wide(data):
    """Convert from long format to SCP DE wide format

    (Long format is typical uploaded, but this module internally uses wide.)
    """
    print("data data.dtypes")
    print(data.dtypes)
    print('A 0')
    metrics = list(data.columns[3:])
    print('A 1')
    data["combined"] = data["group"] + "--" + data["comparison_group"]
    frames = []
    print('A 2')
    # E.g.
    # gene  group   comparison_group    log2foldchange  pvals_adj   qvals   mean    cat dog
    # metrics = ["log2foldchange", "pvals_adj", "qvals", "mean", "cat", "dog"]
    for metric in metrics:
        wide_metric = pd.pivot(data, index="genes", columns="combined", values=metric)
        wide_metric = wide_metric.add_suffix(f"--{metric}")
        frames.append(wide_metric)

    print('A 3')
    result = pd.concat(frames, axis=1, join="inner")
    print('A 4')
    result.columns.name = " "
    print('result before', result)
    result = result.reset_index()
    print('result', result)
    print('A 5')
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
                tsv_output.writerow([file_names_pairwise[value][0], file_names_pairwise[value][1],])


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

    # Rank significance 1st (ultimately ranked 2nd)
    comparison_metrics = sorted(
        comparison_metrics,
        key=lambda x: x.split('--')[-1] == significance
    )

    # Rank size 1st (ultimately ranked 1st)
    comparison_metrics = sorted(
        comparison_metrics,
        key=lambda x: x.split('--')[-1] == size
    )

    comparison_metrics.reverse()

    return comparison_metrics


def sort_metrics(metrics, size, significance):
    """Like `sort_comparison_metrics`, but for bare metrics
    """

    # Sort alphabetically
    metrics = sorted(metrics)

    # Rank significance 1st (ultimately ranked 2nd)
    metrics = sorted(
        metrics,
        key=lambda x: x == significance
    )

    # Rank size 1st (ultimately ranked 1st)
    metrics = sorted(
        metrics,
        key=lambda x: x == size
    )

    metrics.reverse()

    return metrics

# note: my initial files had pval, qval, logfoldchanges.
# David's files have qval, mean, logfoldchanges.
# For the purposes of this validation I will be using his column values/formatting.


def validate_size_and_significance(metrics, size, significance, logger):
    """Locally log whether size and/or significance are detected among metrics

    TODO:
        - Log to Sentry / Mixpanel
    """
    has_size = any([metric.split('--')[-1] == size for metric in metrics])
    has_significance = any([metric.split('--')[-1] == significance for metric in metrics])

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
        logger.info(f'Found size ("{size}") and significance ("{significance}") metrics {in_headers}')


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
            item = item.replace("'", "")  # Remove quotes in e.g. 'type_0'--'type_1'--qval
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
    findallmarkers_headers = ['p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene']
    is_seurat_findallmarkers = (
        len(headers) == len(findallmarkers_headers) and all(headers == findallmarkers_headers)
    )
    return is_seurat_findallmarkers


class AuthorDifferentialExpression:
    # TODO: reorder author's columns in input file so output is logfoldchanges qval mean
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
        size_metric,
        significance_metric
    ):
        # AuthorDifferentialExpression.de_logger.info(
        #    "Initializing DifferentialExpression instance"
        # )
        self.cluster_name = sanitize_string(cluster_name)
        self.annotation = sanitize_string(annotation_name)
        self.accession = study_accession
        self.annot_scope = annotation_scope
        self.method = method
        self.size_metric = size_metric
        self.significance_metric = significance_metric
        self.stem = f"{self.cluster_name}--{self.annotation}"

        author_de_file_gcs_url = differential_expression_file
        allowed_file_types = AuthorDifferentialExpression.ALLOWED_FILE_TYPES
        raw_de_file_obj = IngestFiles(author_de_file_gcs_url, allowed_file_types)
        author_file_handle, local_file_path = IngestFiles.resolve_path(
            raw_de_file_obj, author_de_file_gcs_url
        )
        self.author_de_file = local_file_path

    def execute(self):
        clean_val = []
        clean_val_p = []
        metrics = []
        file_path = self.author_de_file

        # TODO:
        #   - Throw well-formatted error if file type not in ALLOWED_FILE_TYPES
        #   - Consider if / how this merges with centralized handling in ingest_files.py
        file_type = IngestFiles.get_file_type(file_path)[0]
        if file_type == 'text/tab-separated-values':
            delimiter = '\t'
        else:
            delimiter = ','
        data = pd.read_csv(file_path, delimiter)

        first_cols = data.columns

        is_seurat_findallmarkers = detect_seurat_findallmarkers(first_cols)

        if first_cols[0] == "genes" and first_cols[1] == "group" and first_cols[2] == "comparison_group":
            # Raw data is in long format
            wide_format = convert_long_to_wide(data)
            data_by_col = get_data_by_column(wide_format)
        elif '--' in first_cols[1]:
            # Raw data is in wide format
            data_by_col = get_data_by_column(data)
        elif is_seurat_findallmarkers:
            print('Data is in Seurat FindAllMarkers() format')
            # Raw data is a simple plaintext export from Seurat
            wide_format = convert_seurat_findallmarkers_to_wide(data)
            print('wide_format', wide_format)
            data_by_col = get_data_by_column(wide_format)

        col, genes, rest, split_values = data_by_col
        pairwise = split_values["pairwise"]
        one_vs_rest = split_values["one_vs_rest"]

        logger = AuthorDifferentialExpression.dev_logger
        if len(one_vs_rest) != 0:
            groups, clean_val, metrics = get_groups_and_metrics(
                one_vs_rest, self.size_metric, self.significance_metric, logger
            )
            self.generate_result_files(one_vs_rest, genes, rest, groups, clean_val, metrics)

        if len(pairwise) != 0:
            groups_p, clean_val_p, metrics = get_groups_and_metrics(
                pairwise, self.size_metric, self.significance_metric, logger
            )
            self.generate_result_files(pairwise, genes, rest, groups_p, clean_val_p, metrics)
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
                comparison_metrics = i[x: x + num_metrics]
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
                comparison = check_group(comparisons_dict, comparison_metric)

                # Numerical values for metric in this comparison
                values = rest[comparison_metric].tolist()

                metric_values.append(values)

            rows = metric_values
            rows.insert(0, genes)
            rows_by_comparison[comparison] = rows

        headers = sort_metrics(metrics, self.size_metric, self.significance_metric)
        headers.insert(0, "genes")

        for comparison in rows_by_comparison:
            arr = np.array(rows_by_comparison[comparison])
            t_arr = arr.transpose()
            inner_df = pd.DataFrame(data=t_arr, columns=headers)

            if "rest" in comparison:
                comparison = comparison.split("--")[0]

            else:
                group = comparison.split("--")[0]
                comparison_group = comparison.split("--")[1]
                sorted_list = sort_comparison([group, comparison_group])
                comparison = f'{sorted_list[0]}--{sorted_list[1]}'

            comparison = '--'.join([sanitize_string(group) for group in comparison.split('--')])

            tsv_name = f'{self.stem}--{comparison}--{self.annot_scope}--{self.method}.tsv'
            inner_df.to_csv(tsv_name, sep='\t')

            print(f"Wrote TSV: {tsv_name}")
