"""Ingest differential expression uploaded by study owner or editor

EXAMPLE:
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 ingest_differential_expression --annotation-name General_Celltype --annotation-type group --annotation-scope study --cluster-name cluster_umap_txt --study-accession SCPdev --ingest-differential-expression --differential-expression-file gs://fc-febd4c65-881d-497f-b101-01a7ec427e6a/author_de_test_data_human_milk_All_Cells_UMAP_General_celltype.csv --method wilcoxon
"""

import pandas as pd
import numpy as np
import csv
import logging

from monitor import setup_logger, log_exception
from de import DifferentialExpression
from ingest_files import IngestFiles


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


def convert_long_to_wide(data):
    """Convert from long format to wide format

    (Long format is typical uploaded, but this module internally uses wide.)
    """
    data["combined"] = data['group'] + "--" + data["comparison_group"]
    frames = []
    metrics = ["logfoldchanges", "qval", "mean"]
    for metric in metrics:
        wide_metric = pd.pivot(data, index="genes", columns="combined", values=metric)
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
    """
    create manifest file of each comparison in the initial data
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


# note: my initial files had pval, qval, logfoldchanges.
# David's files have qval, mean, logfoldchanges.
# For the purposes of this validation I will be using his column values/formatting.


def get_groups_and_metrics(raw_column_names):
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
    # metrics: ['qval', 'logfoldchanges', 'mean']
    for split_header in split_headers:
        [group, comparison_group, metric] = split_header
        if group not in groups:
            groups.append(group)
        if comparison_group not in groups:
            groups.append(comparison_group)
        if metric not in metrics:
            metrics.append(metric)

    # TODO: Report this error to Sentry
    if ("logfoldchanges" not in metrics) or ("qval" not in metrics):
        raise Exception("Comparisons must include at least logfoldchanges and qval to be valid")

    return groups, split_headers, metrics


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
        **kwargs
    ):
        # AuthorDifferentialExpression.de_logger.info(
        #    "Initializing DifferentialExpression instance"
        # )
        self.cluster_name = DifferentialExpression.sanitize_strings(cluster_name)
        self.annotation = DifferentialExpression.sanitize_strings(annotation_name)
        self.kwargs = kwargs
        self.accession = self.kwargs["study_accession"]
        self.annot_scope = self.kwargs["annotation_scope"]
        self.method = self.kwargs["method"]
        self.stem = f"{self.cluster_name}--{self.annotation}"

        author_de_file_gcs_url = self.kwargs["differential_expression_file"]
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

        data = pd.read_csv(file_path)
        first_cols = data.columns
        if first_cols[0] == "genes" and first_cols[1] == "group" and first_cols[2] == "comparison_group":
            wide_format = convert_long_to_wide(data)
            data_by_col = get_data_by_column(wide_format)
        else:
            data_by_col = get_data_by_column(data)

        col, genes, rest, split_values = data_by_col
        pairwise = split_values["pairwise"]
        one_vs_rest = split_values["one_vs_rest"]

        if len(one_vs_rest) != 0:
            groups, clean_val, metrics = get_groups_and_metrics(one_vs_rest)
            self.generate_result_files(one_vs_rest, genes, rest, groups, clean_val, metrics)

        if len(pairwise) != 0:
            groups_p, clean_val_p, metrics = get_groups_and_metrics(pairwise)
            self.generate_result_files(pairwise, genes, rest, groups_p, clean_val_p, metrics)

        generate_manifest(self.stem, clean_val, clean_val_p, metrics)

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

        names_dict = {}
        all_group = []
        for i in range(len(groups)):
            curr_group = groups[i]
            type_group = []

            for j in range(len(clean_val)):
                curr_val = clean_val[j][0]
                comp_val = clean_val[j][1]
                file_naming = f"{curr_val}--{comp_val}"
                names_dict[file_naming] = []
                real_title = col[j]
                if curr_group == curr_val:
                    type_group.append(real_title)

            all_group.append(type_group)
        all_group_fin = [ele for ele in all_group if ele != []]
        grouped_lists = []

        # TODO: fix sorting error here. if you have comparison set 1 with foo and bar,
        # then you have comparison set 2 with bar and baz, error is triggered. adjust sorting method
        for i in all_group_fin:
            for j in range(0, len(i), 3):
                x = j
                grouped_lists.append(i[x: x + 3])

        for i in names_dict:
            for j in grouped_lists:
                for k in j:
                    if i in k:
                        names_dict[i].append(k)

        # Now we have all the columns grouped in lists by pairwise comparison, with qval, logfoldchanges, mean
        # have to pair with corresponding gene for that row
        # dictionary format:
        # comparison name: [[gene, qval, logfoldchanges, mean] [gene, qval, logfoldchanges, mean] etc...]
        keys = names_dict.keys()
        file_d = dict.fromkeys(keys, [])
        for i in grouped_lists:
            list_i = []
            for j in i:
                f_name = check_group(names_dict, j)
                col_v = rest[j].tolist()
                list_i.append(col_v)
            file_d[f_name] = genes, list_i[0], list_i[1], list_i[2]

        headers = metrics
        headers.insert(0, "genes")

        final_files_to_find = []
        for i in file_d:
            arr = np.array(file_d[i])
            t_arr = arr.transpose()
            inner_df = pd.DataFrame(data=t_arr, columns=headers)

            if "rest" in i:
                i = i.split("--")[0]

            else:
                first_group = i.split("--")[0]
                second_group = i.split("--")[1]
                sorted_list = sort_comparison([first_group, second_group])
                i = f'{sorted_list[0]}--{sorted_list[1]}'

            comparison = '--'.join([DifferentialExpression.sanitize_strings(group) for group in i.split('--')])

            tsv_name = f'{self.stem}--{comparison}--{self.annot_scope}--{self.method}.tsv'
            final_files_to_find.append(tsv_name)
            inner_df.to_csv(tsv_name, sep='\t')
        return final_files_to_find
