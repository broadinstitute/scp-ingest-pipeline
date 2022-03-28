import logging
import numpy as np
import pandas as pd
import scanpy as sc

try:
    from monitor import setup_logger, log_exception
    from annotations import Annotations
    from clusters import Clusters
    from cell_metadata import CellMetadata
    from ingest_files import IngestFiles
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .monitor import setup_logger, log_exception
    from .annotations import Annotations
    from .clusters import Clusters
    from .cell_metadata import CellMetadata
    from .ingest_files import IngestFiles


class DifferentialExpression:
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")
    de_logger = setup_logger(
        __name__ + ".de_logger",
        "de_log.txt",
        level=logging.INFO,
        format="support_configs",
    )

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
        self.kwargs = kwargs
        self.accession = self.kwargs["study_accession"]
        # only used in output filename, removing spaces
        self.cluster_name = self.kwargs["name"].replace(" ", "_")
        self.method = self.kwargs["method"]

        if matrix_file_type == "mtx":
            self.genes_path = self.kwargs["gene_file"]
            genes_ingest_file = IngestFiles(self.genes_path, self.ALLOWED_FILE_TYPES)
            self.genes_file = genes_ingest_file.resolve_path(self.genes_path)[0]
            self.genes: List[str] = [
                g.strip().strip('"') for g in self.genes_file.readlines()
            ]

            self.barcodes_path = self.kwargs["barcode_file"]
            barcodes_ingest_file = IngestFiles(
                self.barcodes_path, self.ALLOWED_FILE_TYPES
            )
            self.barcodes_file = barcodes_ingest_file.resolve_path(self.barcodes_path)[
                0
            ]
            self.barcodes: List[str] = [
                c.strip().strip('"') for c in self.barcodes_file.readlines()
            ]
        DifferentialExpression.de_logger.info(f"DifferentialExpression initialized")

    @staticmethod
    def get_cluster_cells(cluster_cells):
        """ ID cells in cluster file """
        # cluster_cells.tolist() yields a list of lists that needs to be flattened
        # using extend converts a single-value list to a plain value
        cluster_cell_values = cluster_cells.tolist()
        cluster_cell_list = []
        for value in cluster_cell_values:
            cluster_cell_list.extend(value)
        return cluster_cell_list

    @staticmethod
    def determine_dtypes(headers, annot_types):
        """ use SCP TYPE data to coerce data to proper dtypes:
                numeric-like group annotations
                missing values in group annotations (to avoid NaN)
        """
        dtype_info = dict(zip(headers, annot_types))
        dtypes = {}
        for header, type_info in dtype_info.items():
            if type_info == "group":
                dtypes[header] = "string"
        return dtypes

    @staticmethod
    def load_raw_annots(metadata_file_path, allowed_file_types, headers, dtypes):
        """ using SCP metadata header lines
            create properly coerced pandas dataframe of all study metadata
        """
        annot_redux = IngestFiles(metadata_file_path, allowed_file_types)
        annot_file_type = annot_redux.get_file_type(metadata_file_path)[0]
        annot_file_handle = annot_redux.open_file(metadata_file_path)[1]
        raw_annots = annot_redux.open_pandas(
            metadata_file_path,
            annot_file_type,
            open_file_object=annot_file_handle,
            names=headers,
            skiprows=2,
            index_col=0,
            dtype=dtypes,
        )
        return raw_annots

    @staticmethod
    def subset_annots(metadata, de_cells):
        """ subset metadata based on cells in cluster
        """
        DifferentialExpression.de_logger.info(
            f"subsetting metadata on cells in clustering"
        )
        dtypes = DifferentialExpression.determine_dtypes(
            metadata.headers, metadata.annot_types
        )
        raw_annots = DifferentialExpression.load_raw_annots(
            metadata.file_path, metadata.ALLOWED_FILE_TYPES, metadata.headers, dtypes
        )
        cluster_annots = raw_annots[raw_annots.index.isin(de_cells)]
        return cluster_annots

    @staticmethod
    def order_annots(metadata, adata_cells):
        """ order metadata based on cells order in matrix
        """
        matrix_cell_order = adata_cells.tolist()
        return metadata.reindex(matrix_cell_order)

    @staticmethod
    def subset_adata(adata, de_cells):
        """ subset adata object based on cells in cluster
        """
        DifferentialExpression.de_logger.info(
            f"subsetting matrix on cells in clustering"
        )
        matrix_subset_list = np.in1d(adata.obs_names, de_cells)
        adata = adata[matrix_subset_list]
        return adata

    def execute_de(self):
        if self.matrix_file_type == "mtx":
            self.prepare_h5ad(
                self.cluster,
                self.metadata,
                self.matrix_file_path,
                self.matrix_file_type,
                self.annotation,
                self.accession,
                self.genes,
                self.barcodes,
            )
            DifferentialExpression.de_logger.info("preparing DE on sparse matrix")
        else:
            self.run_h5ad(
                self.cluster,
                self.metadata,
                self.matrix_file_path,
                self.matrix_file_type,
                self.annotation,
                self.accession,
                self.cluster_name,
                self.method,
            )
            DifferentialExpression.de_logger.info("preparing DE on dense matrix")

    @staticmethod
    def run_h5ad(
        cluster,
        metadata,
        matrix_file_path,
        matrix_file_type,
        annotation,
        study_accession,
        cluster_name,
        method,
        genes=None,
        barcodes=None,
    ):

        de_cells = DifferentialExpression.get_cluster_cells(cluster.file['NAME'].values)
        de_annots = DifferentialExpression.subset_annots(metadata, de_cells)

        if matrix_file_type == "dense":
            # will need try/except (SCP-4205)
            adata = sc.read(matrix_file_path)
        else:
            # MTX DE UNTESTED (SCP-4203)
            # will want try/except here to catch failed data object composition
            adata = sc.read_mtx(matrix_file_path)
            # For AnnData, obs are cells and vars are genes
            # BUT transpose needed for both dense and sparse
            # so transpose step is after this data object composition step
            # therefore the assignements below are the reverse of expected
            adata.obs_names = genes
            adata.var_names = barcodes

        adata = adata.transpose()

        adata = DifferentialExpression.subset_adata(adata, de_cells)

        # will need try/except (SCP-4205)
        adata.obs = DifferentialExpression.order_annots(de_annots, adata.obs_names)

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        rank_key = "rank." + annotation + "." + method
        DifferentialExpression.de_logger.info(f"calculating DE")
        try:
            sc.tl.rank_genes_groups(
                adata,
                annotation,
                key_added=rank_key,
                use_raw=False,
                method=method,
                pts=True,
            )
        except KeyError:
            msg = f"Missing expected annotation in metadata: {annotation}, unable to calculate DE."
            log_exception(
                DifferentialExpression.dev_logger, DifferentialExpression.de_logger, msg
            )
            raise KeyError(msg)

        groups = np.unique(adata.obs[annotation]).tolist()
        for group in groups:
            group_filename = group.replace(" ", "_")
            DifferentialExpression.de_logger.info(f"Writing DE output for {str(group)}")
            rank = sc.get.rank_genes_groups_df(adata, key=rank_key, group=str(group))

            out_file = (
                f'{cluster_name}--{annotation}--{str(group_filename)}--{method}.tsv'
            )
            # Round numbers to 4 significant digits while respecting fixed point
            # and scientific notation (note: trailing zeros are removed)
            rank.to_csv(out_file, sep='\t', float_format='%.4g')

        # Provide h5ad of DE analysis as reference computable object
        # DifferentialExpression.de_logger.info(f"Writing DE h5ad file")
        # file_name = f'{study_accession}_{cluster_name}_to_DE.h5ad'
        # adata.write_h5ad(file_name)

        DifferentialExpression.de_logger.info(f"DE processing complete")

