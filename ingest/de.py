import numpy as np
import pandas as pd
import scanpy as sc

try:
    from annotations import Annotations
    from clusters import Clusters
    from cell_metadata import CellMetadata
    from ingest_files import IngestFiles
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .annotations import Annotations
    from .clusters import Clusters
    from .cell_metadata import CellMetadata
    from .ingest_files import IngestFiles


class DifferentialExpression:
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    def __init__(
        self,
        cluster,
        cell_metadata,
        matrix_file_path,
        matrix_file_type,
        annotation,
        **kwargs,
    ):
        self.cluster = cluster
        self.metadata = cell_metadata
        self.annotation = annotation
        self.matrix_file_path = matrix_file_path
        self.matrix_file_type = matrix_file_type

        if matrix_file_type == "mtx":
            self.genes_path = kwargs.pop("gene_file")
            genes_ingest_file = IngestFiles(genes_path, self.ALLOWED_FILE_TYPES)
            self.genes_file = genes_ingest_file.resolve_path(self.genes_path)[0]
            self.genes: List[str] = [
                g.strip().strip('"') for g in self.genes_file.readlines()
            ]

            self.barcodes_path = kwargs.pop("barcode_file")
            barcodes_ingest_file = IngestFiles(barcodes_path, self.ALLOWED_FILE_TYPES)
            self.barcodes_file = barcodes_ingest_file.resolve_path(self.barcodes_path)[
                0
            ]
            self.cells: List[str] = [
                c.strip().strip('"') for c in self.barcodes_file.readlines()
            ]

    @staticmethod
    def get_cluster_cells(cluster):
        """ ID cells in cluster file """
        cluster_cell_values = cluster.file['NAME'].values.tolist()
        cluster_cell_list = []
        for value in cluster_cell_values:
            cluster_cell_list.extend(value)
        return cluster_cell_list

    @staticmethod
    def determine_dtypes(metadata):
        """ use SCP TYPE data to coerce data to proper dtypes:
                numeric-like group annotations
                missing values in group annotations (to avoid NaN)
        """
        dtype_info = dict(zip(metadata.headers, metadata.annot_types))
        dtypes = {}
        for header, type_info in dtype_info.items():
            if type_info == "group":
                dtypes[header] = "string"
        return dtypes

    @staticmethod
    def load_raw_annots(metadata, dtypes):
        """ using SCP metadata header lines
            create properly coerced pandas dataframe of all study metadata
        """
        annot_redux = IngestFiles(metadata.file_path, metadata.ALLOWED_FILE_TYPES)
        annot_file_type = annot_redux.get_file_type(metadata.file_path)[0]
        annot_file_handle = annot_redux.open_file(metadata.file_path)[1]
        raw_annots = annot_redux.open_pandas(
            metadata.file_path,
            annot_file_type,
            open_file_object=annot_file_handle,
            names=metadata.headers,
            skiprows=2,
            index_col=0,
            dtype=dtypes,
        )
        return raw_annots

    @staticmethod
    def prepare_annots(metadata, de_cells):
        """ subset metadata based on cells in cluster
        """
        dtypes = DifferentialExpression.determine_dtypes(metadata)
        raw_annots = DifferentialExpression.load_raw_annots(metadata, dtypes)
        cluster_annots = raw_annots[raw_annots.index.isin(de_cells)]
        return cluster_annots

    def execute_de(self):
        self.prepare_h5ad(
            self.cluster, self.metadata, self.matrix_file_path, self.annotation
        )

    @staticmethod
    def prepare_h5ad(cluster, metadata, matrix_file_path, annotation):
        """
        """
        de_cells = DifferentialExpression.get_cluster_cells(cluster)
        de_annots = DifferentialExpression.prepare_annots(metadata, de_cells)

        # dense matrix
        data = sc.read(matrix_file_path)
        adata = data.transpose()

        # subset matrix based on cells in cluster
        matrix_subset_list = np.in1d(adata.obs_names, de_cells)
        adata = adata[matrix_subset_list]

        adata.obs = de_annots

        # ideally include cluster file name in either filename or as directory
        # have rails provide name as input?
        file_name = f'/Volumes/jlc2T/active/SCP1677/DE/{metadata.study_accession}_raw_to_DE.h5ad'
        adata.write_h5ad(file_name)

        adata.raw = adata
        adata.write_h5ad(file_name)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # adata.write_h5ad(file_name)
        rank_method = 'wilcoxon'
        rank_key = "rank." + annotation + "." + rank_method
        sc.tl.rank_genes_groups(
            adata,
            annotation,
            key_added=rank_key,
            use_raw=False,
            method=rank_method,
            pts=True,
        )
        adata.write_h5ad(file_name)
        groups = np.unique(adata.obs[annotation]).tolist()
        for group in groups:
            rank = sc.get.rank_genes_groups_df(adata, key=rank_key, group=str(group))
            out_file = f'/Volumes/jlc2T/active/SCP1677/DE/{annotation}-{str(group)}-{rank_method}.tsv'
            # when ready, add compression='gzip'
            rank.to_csv(out_file, sep='\t')

        rank_method = 't-test'
        rank_key = "rank." + annotation + "." + rank_method
        sc.tl.rank_genes_groups(
            adata,
            annotation,
            key_added=rank_key,
            use_raw=False,
            method=rank_method,
            pts=True,
        )
        adata.write_h5ad(file_name)
        print("bar")

