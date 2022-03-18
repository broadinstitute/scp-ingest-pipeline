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
        # self.matrix_file = IngestFiles(
        #     matrix_file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        # )
        self.annotation = annotation

        if matrix_file_type == "mtx":
            genes_path = kwargs.pop("gene_file")
            genes_ingest_file = IngestFiles(genes_path, self.ALLOWED_FILE_TYPES)
            self.genes_file = genes_ingest_file.resolve_path(genes_path)[0]
            self.genes: List[str] = [
                g.strip().strip('"') for g in self.genes_file.readlines()
            ]

            barcodes_path = kwargs.pop("barcode_file")
            barcodes_ingest_file = IngestFiles(barcodes_path, self.ALLOWED_FILE_TYPES)
            self.barcodes_file = barcodes_ingest_file.resolve_path(barcodes_path)[0]
            self.cells: List[str] = [
                c.strip().strip('"') for c in self.barcodes_file.readlines()
            ]

        # ID cells in cluster
        cluster_cell_values = self.cluster.file['NAME'].values.tolist()
        cluster_cell_list = []
        for value in cluster_cell_values:
            cluster_cell_list.extend(value)

        # coerce numeric-like group annotations ?and groups with missing values?
        dtype_info = dict(zip(self.metadata.headers, self.metadata.annot_types))
        dtypes = {}
        for header, type_info in dtype_info.items():
            if type_info == "group":
                dtypes[header] = "string"

        annot_redux = IngestFiles(self.metadata.file_path, self.ALLOWED_FILE_TYPES)
        annot_file_type = annot_redux.get_file_type(self.metadata.file_path)[0]
        annot_file_handle = annot_redux.open_file(self.metadata.file_path)[1]
        raw_annots = annot_redux.open_pandas(
            self.metadata.file_path,
            annot_file_type,
            open_file_object=annot_file_handle,
            names=self.metadata.headers,
            skiprows=2,
            index_col=0,
            dtype=dtypes,
        )

        # subset metadata based on cells in cluster
        cluster_annots = raw_annots[raw_annots.index.isin(cluster_cell_list)]

        # dense matrix
        data = sc.read(matrix_file_path)
        adata = data.transpose()

        # subset matrix based on cells in cluster
        matrix_subset_list = np.in1d(adata.obs_names, cluster_cell_list)
        adata = adata[matrix_subset_list]

        adata.obs = cluster_annots

        # ideally include cluster file name in either filename or as directory
        # have rails provide name as input?
        file_name = f'/Volumes/jlc2T/active/SCP1671/DE/{self.metadata.study_accession}_raw_to_DE.h5ad'
        adata.write_h5ad(file_name)

        adata.raw = adata
        adata.write_h5ad(file_name)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # adata.write_h5ad(file_name)
        rank_method = 'wilcoxon'
        rank_key = "rank." + self.annotation + "." + rank_method
        sc.tl.rank_genes_groups(
            adata,
            self.annotation,
            key_added=rank_key,
            use_raw=False,
            method=rank_method,
            pts=True,
        )
        adata.write_h5ad(file_name)
        groups = np.unique(adata.obs[self.annotation]).tolist()
        for group in groups:
            rank = sc.get.rank_genes_groups_df(adata, key=rank_key, group=str(group))
            out_file = f'/Volumes/jlc2T/active/SCP1671/DE/{self.annotation}-{str(group)}-{rank_method}.tsv'
            # when ready, add compression='gzip'
            rank.to_csv(out_file, sep='\t')

        rank_method = 't-test'
        rank_key = "rank." + self.annotation + "." + rank_method
        sc.tl.rank_genes_groups(
            adata,
            self.annotation,
            key_added=rank_key,
            use_raw=False,
            method=rank_method,
            pts=True,
        )
        adata.write_h5ad(file_name)
        print("bar")

        return 0
