try:
    from ingest_files import IngestFiles
    from monitor import log_exception
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import IngestFiles
    from .monitor import log_exception


class H5adIngestor(IngestFiles):
    ALLOWED_FILE_TYPES = ['application/x-hdf5']

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        IngestFiles.__init__(
            self, file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        )
        pass

    def obtain_adata(self):
        try:
            self.adata = self.open_file(self.file_path)[0]
            print(self.adata)
            IngestFiles.dev_logger.info(str(self.adata))
        except ValueError as e:
            raise ValueError(e)

    def validate(self):
        """
        Currently, file passes "basic validation" if file
        can be opened by scanpy
        """
        try:
            self.adata = self.obtain_adata()
            return True
        except ValueError:
            return False

