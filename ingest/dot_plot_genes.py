import pandas as pd
import json
import gzip
import os
import glob
import re
import multiprocessing
import sys
import datetime
import string
import random
from dateutil.relativedelta import relativedelta
from functools import partial
from bson.objectid import ObjectId

from annotations import Annotations
from expression_writer import ExpressionWriter
from writer_functions import get_cluster_cells
from monitor import setup_logger, bypass_mongo_writes
from mongo_connection import MongoConnection, graceful_auto_reconnect


class DotPlotGenes:
    COLLECTION_NAME = "dot_plot_genes"
    BATCH_SIZE = 50
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    EXP_WRITER_SETTINGS = {"output_format": "dict", "sparse": True, "delocalize": False}
    denominator = 2 if re.match('darwin', sys.platform) else 1
    num_cores = int(multiprocessing.cpu_count() / denominator) - 1
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")

    def __init__(
            self,
            study_id,
            study_file_id,  # expression matrix file
            cluster_group_id,
            cluster_file,
            cell_metadata_file,
            matrix_file_path,
            matrix_file_type,
            **kwargs,
    ):
        self.study_id = study_id
        self.study_file_id = study_file_id
        self.cluster_group_id = cluster_group_id
        self.cluster_file = cluster_file
        self.cell_metadata_file = cell_metadata_file
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
        # a random 6-letter slug is appended to the end to avoid directory collisions when running in CI
        random_slug = ''.join(random.sample(string.ascii_letters, 6))
        self.cluster_name = f"cluster_entry_{self.cluster_group_id}_{random_slug}"
        self.output_path = f"{self.cluster_name}/dot_plot_genes"
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
        self.dev_logger.info(f"preprocessing annotation data from {self.cell_metadata_file}")
        cell_metadata = Annotations(self.cell_metadata_file, self.ALLOWED_FILE_TYPES)
        cell_metadata.preprocess(False)
        valid_metadata = [
            [column[0], 'study', cell_metadata] for column in cell_metadata.file.columns if
            self.annotation_is_valid(column, cell_metadata)
        ]
        # check for annotations in cluster file
        cluster_annots = Annotations(self.cluster_file, self.ALLOWED_FILE_TYPES)
        cluster_annots.preprocess(False)
        valid_cluster_annots = [
            [column[0], 'cluster', cluster_annots] for column in cluster_annots.file.columns if
            self.annotation_is_valid(column, cluster_annots)
        ]
        all_annotations = valid_metadata + valid_cluster_annots
        for annotation_name, annotation_scope, source_data in all_annotations:
            self.add_annotation_to_map(annotation_name, annotation_scope, source_data)
        self.dev_logger.info(f"Annotation data preprocessing for {self.cell_metadata_file} complete")

    def add_annotation_to_map(self, annotation_name, annotation_scope, source_data):
        """
        Take an individual annotation, filter cells and add it to the annotation_map dictionary
        :param annotation_name: (str) name of annotation
        :param annotation_scope: (str) scope of annotation, either study or cluster
        :param source_data: (DataFrame) pandas dataframe of source data
        """
        annotation_id = f"{annotation_name}--group--{annotation_scope}"
        groups = dict(source_data.file['NAME']['TYPE'].groupby(source_data.file[annotation_name]['group']).unique())
        self.dev_logger.info(f"reading {annotation_name}, found {len(groups)} labels")
        self.annotation_map[annotation_id] = {}
        for group in groups:
            filtered_cells = set(self.cluster_cells).intersection(set(groups[group]))
            self.annotation_map[annotation_id][group] = filtered_cells

    def render_gene_expression(self):
        """
        Render gene-level expression for all cells in the given cluster, filtering for non-zero expression
        """
        self.dev_logger.info(f"Rendering cluster-filtered gene expression from {self.matrix_file_path}")
        self.exp_writer.render_artifacts(**self.EXP_WRITER_SETTINGS)
        self.dev_logger.info(f"Gene expression rendering for {self.matrix_file_path} complete")

    def preprocess(self):
        """
        Preprocess all data in preparation of creating DotPlotGene entries
        """
        self.set_annotation_map()
        self.render_gene_expression()
        self.dev_logger.info("All data preprocessing complete")

    @staticmethod
    def annotation_is_valid(column, source_data):
        """
        Determine if a given column in an annotation file is valid
        must be group-based and have between 2 and 250 values
        :param column: (str) name of column
        :param source_data: (DataFrame) pandas dataframe of source data
        :return: (bool)
        """
        viz_range = range(2, 200, 1)
        column_name, annotation_type = column
        return annotation_type == 'group' and len(source_data.file[column_name][annotation_type].unique()) in viz_range

    @staticmethod
    def process_gene(gene_file, output_path, dot_plot_gene, annotation_map):
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
        dot_plot_gene["gene_symbol"] = gene_name
        dot_plot_gene["searchable_gene"] = gene_name.lower()
        gene_dict = DotPlotGenes.get_gene_dict(gene_file)
        exp_scores = DotPlotGenes.get_expression_metrics(gene_dict, annotation_map)
        dot_plot_gene['exp_scores'] = exp_scores
        with gzip.open(f"{output_path}/{gene_name}.json.gz", "wt") as file:
            json.dump(dot_plot_gene, file, separators=(',', ':'))

    @staticmethod
    def get_expression_metrics(gene_doc, annotation_map):
        """
        Set the mean expression and percent cells expressing for all available annotations/labels
        :param gene_doc: (dict) gene-level expression dict
        :param annotation_map (dict) class-level map of all annotations/labels and cells in each
        :return: (dict)
        """
        expression_metrics = {}
        for annotation in annotation_map:
            expression_metrics[annotation] = {}
            for label in annotation_map[annotation]:
                label_cells = annotation_map[annotation][label]
                filtered_expression = DotPlotGenes.filter_expression_for_label(gene_doc, label_cells)
                pct_exp = DotPlotGenes.pct_expression(filtered_expression, label_cells)
                scaled_mean = DotPlotGenes.scaled_mean_expression(filtered_expression, pct_exp)
                expression_metrics[annotation][label] = [scaled_mean, pct_exp]
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
        return re.sub(r'\.json\.gz', '', gene_file_path.split('/')[1])

    @staticmethod
    def get_gene_dict(gene_path):
        """
        Read a gene document and process as a dict
        :param gene_path: (str) path to gzipped gene doc
        :return: (dict)
        """
        return json.load(gzip.open(gene_path, 'rt'))

    @staticmethod
    def to_model(gene_dict):
        """
        Convert a raw dict into a document that can be inserted into MongoDB
        :param gene_dict: (dict) raw processed dot plot gene entry
        :return: (dict) transformed dict with ObjectId entries
        """
        model_dict = gene_dict.copy()
        model_dict['study_id'] = ObjectId(gene_dict['study_id'])
        model_dict['study_file_id'] = ObjectId(gene_dict['study_file_id'])
        model_dict['cluster_group_id'] = ObjectId(gene_dict['cluster_group_id'])
        return model_dict

    @staticmethod
    def scaled_mean_expression(gene_doc, pct_exp):
        """
        Get the scaled mean expression of cells for a given gene
        :param gene_doc: (dict) gene-level significant expression values
        :param pct_exp: (float) percentage of cells expressing for gene to scale mean by
        :return: (float)
        """
        exp_values = pd.DataFrame(gene_doc.values())
        if exp_values.empty:
            return 0.0
        else:
            raw_exp = exp_values.mean()[0]
            return round(raw_exp * pct_exp, 3)

    @staticmethod
    def pct_expression(gene_doc, cells):
        """
        Get the percentage of cells expressing for a given gene relative to the cells in the cluster
        :param gene_doc: (dict) gene-level significant expression values
        :param cells: (list) list of cells for given annotation label
        :return: (float)
        """
        if len(cells) == 0:
            return 0.0
        observed_cells = gene_doc.keys()
        return round(len(observed_cells) / len(cells), 4)

    def process_all_genes(self):
        """
        Parallel function to process all files and render out DotPlotGene dicts
        """
        os.mkdir(self.output_path)
        gene_files = glob.glob(f"{self.cluster_name}/*.json.gz")
        blank_dot_plot_gene = {
            "study_id": self.study_id,
            "study_file_id": self.study_file_id,
            "cluster_group_id": self.cluster_group_id,
            "exp_scores": {}
        }
        self.dev_logger.info(f"beginning parallel rendering of {len(gene_files)} DotPlotGene entries")
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(
            DotPlotGenes.process_gene,
            dot_plot_gene=blank_dot_plot_gene,
            output_path=self.output_path,
            annotation_map=self.annotation_map
        )
        pool.map(processor, gene_files)

    @graceful_auto_reconnect
    def load(self, collection, documents):
        """
        Insert batch of documents into MongoDB
        :param collection: (String) name of db collection (ensure method signature matches @graceful_auto_reconnect)
        :param documents: (list) list of rendered documents to insert
        """
        if not bypass_mongo_writes():
            self.mongo_connection._client[collection].insert_many(documents, ordered=False)
        else:
            dev_msg = f"Extracted {len(documents)} DotPlotGenes for {self.matrix_file_path}"
            self.dev_logger.info(dev_msg)

    def transform(self):
        """
        Main handler to process all data and render/insert DotPlotGenes
        """
        start_time = datetime.datetime.now()
        self.dev_logger.info(f"beginning rendering of {self.matrix_file_path} into DotPlotGene entries")
        self.preprocess()
        self.process_all_genes()
        self.dev_logger.info(f"rendering of {self.matrix_file_path} complete, beginning load")
        gene_docs = []
        collection = self.COLLECTION_NAME
        for gene_path in glob.glob(f"{self.output_path}/*.json.gz"):
            rendered_gene = DotPlotGenes.get_gene_dict(gene_path)
            model_dict = DotPlotGenes.to_model(rendered_gene)
            gene_docs.append(model_dict)
            if len(gene_docs) == self.BATCH_SIZE:
                self.load(collection, gene_docs)
                gene_docs.clear()
        if len(gene_docs) > 0:
            self.load(collection, gene_docs)
            gene_docs.clear()
        end_time = datetime.datetime.now()
        time_diff = relativedelta(end_time, start_time)
        self.dev_logger.info(
            f" completed, total runtime: {time_diff.hours}h, {time_diff.minutes}m, {time_diff.seconds}s"
        )
