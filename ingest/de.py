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
        **kwargs
    ):
        self.cluster = cluster
        self.metadata = cell_metadata
        self.matrix_file = IngestFiles(
            matrix_file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        )
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

        print("foo")
