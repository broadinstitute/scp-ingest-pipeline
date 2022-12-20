import pandas as pd  # NOqa: F821

try:
    from ingest_files import IngestFiles
    from monitor import log_exception
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import IngestFiles
    from .monitor import log_exception


class AnnDataIngestor(IngestFiles):
    ALLOWED_FILE_TYPES = ['application/x-hdf5']

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        IngestFiles.__init__(
            self, file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        )

    def obtain_adata(self):
        try:
            adata = self.open_file(self.file_path)[0]
            # for faster dev, print adata info to screen, may want to remove in future
            print(adata)
            IngestFiles.dev_logger.info(str(adata))
            return adata
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

    @staticmethod
    def generate_cluster_header(adata, clustering_name, study_file_id):
        """
        Based on clustering dimensions, write clustering NAME line to file
        """
        dim = ['NAME', 'X', 'Y']
        clustering_dimension = adata.obsm[clustering_name].shape[1]
        if clustering_dimension == 3:
            headers = dim.append('Z')
        elif clustering_dimension == 2:
            headers = dim
        elif clustering_dimension > 3:
            msg = f"Too many dimensions for visualization in obsm \"{clustering_name}\", found {clustering_dimension}, expected 2 or 3."
            raise ValueError(msg)
        else:
            msg = f"Too few dimensions for visualization in obsm \"{clustering_name}\", found {clustering_dimension}, expected 2 or 3."
            raise ValueError(msg)
        file_prefix = f"{study_file_id}.{clustering_name}"
        with open(AnnDataIngestor.set_clustering_filename(file_prefix), "w") as f:
            f.write('\t'.join(headers) + '\n')

    @staticmethod
    def generate_cluster_type_declaration(adata, clustering_name, study_file_id):
        """
        Based on clustering dimensions, write clustering TYPE line to file
        """
        clustering_dimension = adata.obsm[clustering_name].shape[1]
        types = ["TYPE", *["numeric"] * clustering_dimension]
        file_prefix = f"{study_file_id}.{clustering_name}"
        with open(AnnDataIngestor.set_clustering_filename(file_prefix), "a") as f:
            f.write('\t'.join(types) + '\n')

    @staticmethod
    def generate_cluster_body(adata, clustering_name, study_file_id):
        """
        Append clustering data to clustering file
        """
        cluster_cells = pd.DataFrame(adata.obs_names)
        cluster_body = pd.concat(
            [cluster_cells, pd.DataFrame(adata.obsm[clustering_name])], axis=1
        )
        file_prefix = f"{study_file_id}.{clustering_name}"
        pd.DataFrame(cluster_body).to_csv(
            AnnDataIngestor.set_clustering_filename(file_prefix),
            sep="\t",
            mode="a",
            header=None,
            index=False,
        )

    @staticmethod
    def set_clustering_filename(name):
        return f"{name}.cluster.anndata_segment.tsv"

    @staticmethod
    def generate_metadata_file(adata, output_name):
        """
        Generate metadata NAME and TYPE lines
        """
        headers = adata.obs.columns.tolist()
        types = []
        for header in headers:
            if pd.api.types.is_number(adata.obs[header]):
                types.append("NUMERIC")
            else:
                types.append("GROUP")
        headers.insert(0, "NAME")
        types.insert(0, "TYPE")
        with open(output_name, "w") as f:
            f.write('\t'.join(headers) + '\n')
            f.write('\t'.join(types) + '\n')
        adata.obs.to_csv(output_name, sep="\t", mode="a", header=None, index=True)

    @staticmethod
    def clusterings_to_delocalize(arguments, study_file_id):
        # ToDo - check if names using obsm_keys need sanitization
        cluster_file_names = []
        for name in arguments["obsm_keys"]:
            file_prefix = f"{study_file_id}.{name}"
            cluster_file_names.append(
                AnnDataIngestor.set_clustering_filename(file_prefix)
            )
        return cluster_file_names

    @staticmethod
    def delocalize_extracted_files(file_path, study_file_id, files_to_delocalize):
        """ Copy extracted files to study bucket
        """

        for file in files_to_delocalize:
            IngestFiles.delocalize_file(
                study_file_id,
                None,
                file_path,
                file,
                f"_scp_internal/anndata_ingest/{file}",
            )
