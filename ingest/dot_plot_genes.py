import pandas as pd
import json
import gzip
import os
import glob
import re
import multiprocessing
import sys
import datetime
from dateutil.relativedelta import relativedelta
from functools import partial
from bson.objectid import ObjectId

try:
    from annotations import Annotations
    from expression_writer import ExpressionWriter
    from writer_functions import get_cluster_cells
    from monitor import setup_logger, log_exception
    from mongo_connection import MongoConnection, graceful_auto_reconnect
    import config

except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .annotations import Annotations
    from .expression_writer import ExpressionWriter
    from .writer_functions import get_cluster_cells
    from .monitor import setup_logger, log_exception
    from .mongo_connection import MongoConnection, graceful_auto_reconnect
    from . import config

class DotPlotGenes:
    COLLECTION_NAME = "dot_plot_genes"
    BATCH_SIZE = 100
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    EXP_WRITER_SETTINGS = {"output_format": "dict", "sparse": True}
    denominator = 2 if re.match('darwin', sys.platform) else 1
    num_cores = int(multiprocessing.cpu_count() / denominator) - 1
    # Logger provides more details
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")

    def __init__(
            self,
            study_id,
            study_file_id, # expression matrix file
            cluster_group_id,
            cluster_file,
            annotation_file,
            matrix_file_path,
            matrix_file_type,
            **kwargs,
    ):
        self.study_id = study_id
        self.study_file_id = study_file_id
        self.cluster_group_id = cluster_group_id
        self.cluster_file = cluster_file
        self.annotation_file = annotation_file
        self.matrix_file_path = matrix_file_path
        self.matrix_file_type = matrix_file_type
        self.kwargs = kwargs

        if matrix_file_type == "mtx":
            self.genes_path = self.kwargs["gene_file"]
            self.barcodes_path = self.kwargs["barcode_file"]
        else:
            self.genes_path = None
            self.barcodes_path = None

        self.mongo_connection = MongoConnection()

        # the cluster name here is not important, it is only used for directory names
        self.cluster_name = f"cluster_{self.cluster_group_id}"
        self.exp_writer = ExpressionWriter(
            self.matrix_file_path, self.matrix_file_type, self.cluster_file, self.cluster_name, self.genes_path,
            self.barcodes_path
        )

        self.annotation_map = {}
        self.cluster_cells = []

    def set_annotation_map(self):
        """
        Preprocess all associated annotation data to generate the following:
        - list of cluster cells
        - map of all qualifying annotations with list of clusters cells in each annotation label
        """
        self.dev_logger.info(f"getting cluster cells from {self.cluster_file}")
        self.cluster_cells = get_cluster_cells(self.cluster_file)
        self.dev_logger.info(f"preprocessing annotation data from {self.annotation_file}")
        raw_annotations = Annotations(self.annotation_file, self.ALLOWED_FILE_TYPES)
        raw_annotations.preprocess(False)
        valid_annot_names = [
            column[0] for column in raw_annotations.file.columns if
            column[1] == 'group' and
            len(raw_annotations.file[column[0]]['group'].unique()) > 1
        ]
        for annot in valid_annot_names:
            annotation_name = f"{annot}--group--study"
            groups = dict(raw_annotations.file['NAME']['TYPE'].groupby(raw_annotations.file[annot]['group']).unique())
            self.dev_logger.info(f"reading {annotation_name}, found {len(groups)} labels")
            self.annotation_map[annotation_name] = {}
            for group in groups:
                filtered_cells = set(self.cluster_cells).intersection(set(groups[group]))
                self.annotation_map[annotation_name][group] = filtered_cells

        self.dev_logger.info(f"Annotation data preprocessing for {self.annotation_file} complete")

    def render_gene_expression(self):
        """
        Render gene-level expression for all cells in the given cluster, filtering for non-zero expression
        :return:
        """
        self.dev_logger.info(f"Rendering cluster-filtered gene expression from {self.matrix_file_path}")
        self.exp_writer.render_artifacts(**self.EXP_WRITER_SETTINGS)
        self.dev_logger.info(f"Gene expression rendering for {self.matrix_file_path} complete")

    def preprocess(self):
        """
        Preprocess all data in preparation of creating DotPlotGene entries
        :return:
        """
        self.set_annotation_map()
        self.render_gene_expression()
        self.dev_logger.info("All data preprocessing complete")

    @staticmethod
    def process_gene(gene_file, output_path, dot_plot_gene, annotation_map, cluster_cells):
        """
        Read gene-level document and compute both mean and percent cells expressing for all applicable annotations
        :param gene_file: (str) name of gene-level JSON file
        :param output_path (str) path to write output files to
        :param dot_plot_gene (dict) empty DotPlotGene with IDs already populated
        :param annotation_map (dict) class-level map of all annotations/labels and cells in each
        :param cluster_cells (list) list of all cells from cluster
        :return: (dict) fully processed DotPlotGene
        """
        gene_name = DotPlotGenes.get_gene_name(gene_file)
        dot_plot_gene["_id"] = ObjectId()
        dot_plot_gene["gene_symbol"] = gene_name
        dot_plot_gene["searchable_gene"] = gene_name.lower()
        gene_dict = DotPlotGenes.get_gene_dict(gene_file)
        exp_scores = DotPlotGenes.get_expression_metrics(gene_dict, annotation_map, cluster_cells)
        dot_plot_gene['exp_scores'] = exp_scores
        with gzip.open(f"{output_path}/{gene_name}.json", "wt") as file:
            json.dump(dot_plot_gene, file, separators=(',', ':'))

    @staticmethod
    def get_expression_metrics(gene_doc, annotation_map, cluster_cells):
        """
        Set the mean expression and percent cells expressing for all available annotations/labels
        :param gene_doc: (dict) gene-level expression dict
        :param annotation_map (dict) class-level map of all annotations/labels and cells in each
        :param cluster_cells (list) list of all cells from cluster
        :return: (dict)
        """
        expression_metrics = {}
        for annotation in annotation_map:
            expression_metrics[annotation] = {}
            for label in annotation_map[annotation]:
                label_cells = annotation_map[annotation][label]
                filtered_expression = DotPlotGenes.filter_expression_for_label(gene_doc, label_cells)
                mean = DotPlotGenes.mean_expression(filtered_expression)
                pct_exp = DotPlotGenes.pct_expression(filtered_expression, cluster_cells)
                expression_metrics[annotation][label] = [mean, pct_exp]
        return expression_metrics

    @staticmethod
    def filter_expression_for_label(gene_doc, filter_cells):
        """
        Filter gene expression for cells present in a given label
        :param gene_doc: (dict) gene-level expression dict
        :param filter_cells: (list) list of cells to filter on
        :return: (dict) original gene doc filtered by cells from annotation label
        """
        return {cell: exp for cell, exp in gene_doc.items() if cell in filter_cells}

    @staticmethod
    def get_gene_name(gene_file_path):
        """
        Extract gene symbol from filepath
        :param gene_file_path: (str) path to gene JSON file
        :return: (str)
        """
        return re.sub(r'\.json', '', gene_file_path.split('/')[1])

    @staticmethod
    def get_gene_dict(gene_path):
        """
        Read a gene document and process as a dict
        :param gene_path: (str) path to gzipped gene doc
        :return: (dict)
        """
        return json.load(gzip.open(gene_path, 'rt'))

    @staticmethod
    def mean_expression(gene_doc):
        """
        Get the mean expression of cells for a given gene
        :param gene_doc: (dict) gene-level significant expression values
        :return: (float)
        """
        exp_values = pd.DataFrame(gene_doc.values())
        return 0.0 if exp_values.empty else round(exp_values.mean()[0], 3)

    @staticmethod
    def pct_expression(gene_doc, cluster_cells):
        """
        Get the percentage of cells expressing for a given gene relative to the cells in the cluster
        :param gene_doc: (dict) gene-level significant expression values
        :return: (float)
        """
        observed_cells = gene_doc.keys()
        return round(len(observed_cells) / len(cluster_cells), 3)

    def process_all_genes(self):
        """
        Parallel function to process all files and render out DotPlotGene dicts
        :return:
        """
        output_path = f"{self.cluster_name}/dot_plot_genes"
        os.mkdir(output_path)
        gene_files = glob.glob(f"{self.cluster_name}/*.json")
        blank_dot_plot_gene = {
            "study_id": self.study_id,
            "study_file_id": self.study_file_id,
            "cluster_group_id": self.cluster_group_id,
            "exp_scores": {}
        }
        output_path = f"{self.cluster_name}/dot_plot_genes"
        self.dev_logger.info(f"beginning parallel rendering of {len(gene_files)} DotPlotGene entries")
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(
            DotPlotGenes.process_gene,
            dot_plot_gene=blank_dot_plot_gene,
            output_path=output_path,            annotation_map=self.annotation_map,
            cluster_cells=self.cluster_cells
        )
        pool.map(processor, gene_files)

    def run_processor(self):
        """
        Main handler to process all data and render DotPlotGenes
        :return:
        """
        start_time = datetime.datetime.now()
        self.dev_logger.info(f"beginning rendering of {self.matrix_file_path} into DotPlotGene entries")
        self.preprocess()
        self.process_all_genes()
        end_time = datetime.datetime.now()
        time_diff = relativedelta(end_time, start_time)
        self.dev_logger.info(
            f" completed, total runtime: {time_diff.hours}h, {time_diff.minutes}m, {time_diff.seconds}s"
        )
