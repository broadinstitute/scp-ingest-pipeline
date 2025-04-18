import pandas as pd  # NOqa: F821
import os
import gzip
import shutil
import scanpy as sc
import scipy
from scipy.io.mmio import MMFile


# scipy.io.mmwrite uses scientific notation by default
# https://stackoverflow.com/questions/64748513
class MMFileFixedFormat(MMFile):
    def _field_template(self, field, precision):
        # Override MMFile._field_template.
        return f'%.{precision}f\n'


try:
    from ingest_files import DataArray, IngestFiles
    from expression_files.expression_files import GeneExpression
    from monitor import log_exception, bypass_mongo_writes
    from validation.validate_metadata import list_duplicates
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import DataArray, IngestFiles
    from .expression_files.expression_files import GeneExpression
    from .monitor import log_exception, bypass_mongo_writes
    from .validation.validate_metadata import list_duplicates


class AnnDataIngestor(GeneExpression, IngestFiles, DataArray):
    ALLOWED_FILE_TYPES = ['application/x-hdf5']

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        GeneExpression.__init__(self, file_path, study_file_id, study_id)
        IngestFiles.__init__(
            self, file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        )
        self.kwargs = kwargs

    def obtain_adata(self):
        try:
            adata = self.open_file(self.file_path)[0]
            # for faster dev, print adata info to screen, may want to remove in future
            print(adata)
            IngestFiles.dev_logger.info(str(adata))
            return adata
        except ValueError as e:
            raise ValueError(e)

    def basic_validation(self):
        """
        Currently, file passes "basic validation" if file
        can be opened by scanpy
        """
        try:
            self.adata = self.obtain_adata()
            return True
        except ValueError:
            return False

    def validate_raw_location(self):
        """
        Confirm file has data at raw_location
        """
        adata = self.obtain_adata()
        raw_location = self.kwargs.get("raw_location")
        if raw_location is not None:
            try:
                if raw_location == ".raw":
                    if adata.raw is None:
                        msg = f'No data found in .raw slot'
                        log_exception(
                            IngestFiles.dev_logger, IngestFiles.user_logger, msg
                        )
                        raise ValueError(msg)
                else:
                    if raw_location not in adata.layers.keys():
                        msg = f'No data found at adata.layers["{raw_location}"]'
                        log_exception(
                            IngestFiles.dev_logger, IngestFiles.user_logger, msg
                        )
                        raise ValueError(msg)
                return True
            except ValueError:
                return False
        else:
            msg = 'Must specify location of raw counts in AnnData object'
            log_exception(IngestFiles.dev_logger, IngestFiles.user_logger, msg)
            return False

    def create_cell_data_arrays(self):
        """Extract cell name DataArray documents for raw data"""
        adata = self.obtain_adata()
        cells = list(adata.obs_names)
        # use filename denoting a raw 'fragment' to allow successful ingest and downstream queries
        raw_filename = "h5ad_frag.matrix.raw.mtx.gz"
        data_arrays = []
        for data_array in GeneExpression.create_data_arrays(
            name=f"{raw_filename} Cells",
            array_type="cells",
            values=cells,
            linear_data_type="Study",
            linear_data_id=self.study_file_id,
            cluster_name=raw_filename,
            study_file_id=self.study_file_id,
            study_id=self.study_id,
        ):
            data_arrays.append(data_array)

        return data_arrays

    def ingest_raw_cells(self):
        """Insert raw count cells into MongoDB"""
        arrays = self.create_cell_data_arrays()
        if not bypass_mongo_writes():
            self.load(arrays, DataArray.COLLECTION_NAME)
        else:
            dev_msg = f"Extracted {len(arrays)} DataArray for {self.study_file_id}:{arrays[0]['name']}"
            IngestFiles.dev_logger.info(dev_msg)

    @staticmethod
    def generate_cluster_header(adata, clustering_name):
        """
        Based on clustering dimensions, write clustering NAME line to file
        """
        dim = ['NAME', 'X', 'Y']
        clustering_dimension = adata.obsm[clustering_name].shape[1]
        if clustering_dimension == 3:
            dim.append('Z')
            headers = dim
        elif clustering_dimension == 2:
            headers = dim
        elif clustering_dimension > 3:
            msg = f'Too many dimensions for visualization in obsm "{clustering_name}", found {clustering_dimension}, expected 2 or 3.'
            raise ValueError(msg)
        else:
            msg = f'Too few dimensions for visualization in obsm "{clustering_name}", found {clustering_dimension}, expected 2 or 3.'
            raise ValueError(msg)
        filename = AnnDataIngestor.set_clustering_filename(clustering_name)
        with open(filename, "w") as f:
            f.write('\t'.join(headers) + '\n')

    @staticmethod
    def generate_cluster_type_declaration(adata, clustering_name):
        """
        Based on clustering dimensions, write clustering TYPE line to file
        """
        clustering_dimension = adata.obsm[clustering_name].shape[1]
        types = ["TYPE", *["numeric"] * clustering_dimension]
        filename = AnnDataIngestor.set_clustering_filename(clustering_name)
        with open(filename, "a") as f:
            f.write('\t'.join(types) + '\n')

    @staticmethod
    def generate_cluster_body(adata, clustering_name):
        """
        Append clustering data to clustering file
        """
        cluster_cells = pd.DataFrame(adata.obs_names)
        cluster_body = pd.concat(
            [cluster_cells, pd.DataFrame(adata.obsm[clustering_name])], axis=1
        )
        filename = AnnDataIngestor.set_clustering_filename(clustering_name)
        pd.DataFrame(cluster_body).to_csv(
            filename, sep="\t", mode="a", header=None, index=False
        )
        AnnDataIngestor.compress_file(filename)

    @staticmethod
    def set_clustering_filename(name):
        return f"h5ad_frag.cluster.{name}.tsv"

    @staticmethod
    def generate_metadata_file(adata, output_name):
        """
        Generate metadata NAME and TYPE lines
        """
        headers = adata.obs.columns.tolist()
        types = []
        for header in headers:
            if pd.api.types.is_bool_dtype(adata.obs[header]):
                types.append("GROUP")
            elif pd.api.types.is_numeric_dtype(adata.obs[header]):
                types.append("NUMERIC")
            else:
                types.append("GROUP")
        headers.insert(0, "NAME")
        types.insert(0, "TYPE")
        with open(output_name, "w") as f:
            f.write('\t'.join(headers) + '\n')
            f.write('\t'.join(types) + '\n')
        adata.obs.to_csv(output_name, sep="\t", mode="a", header=None, index=True)
        AnnDataIngestor.compress_file(output_name)

    @staticmethod
    def clusterings_to_delocalize(arguments):
        # ToDo - check if names using obsm_keys need sanitization
        cluster_file_names = []
        for name in arguments["obsm_keys"]:
            compressed_file = AnnDataIngestor.set_clustering_filename(name) + ".gz"
            cluster_file_names.append(compressed_file)
        return cluster_file_names

    @staticmethod
    def compress_file(filename):
        with open(filename, 'rb') as file_in:
            compressed_file = filename + '.gz'
            with gzip.open(compressed_file, 'wb') as file_gz:
                shutil.copyfileobj(file_in, file_gz)
        os.remove(filename)

    @staticmethod
    def generate_processed_matrix(adata):
        """
        Generate matrix files with the following file names:
        h5ad_frag.matrix.processed.mtx
        h5ad_frag.barcodes.processed.tsv
        h5ad_frag.features.processed.tsv
        Gzip files for faster delocalization
        """
        if AnnDataIngestor.check_ensembl_index(adata):
            # CELLxGENE indexes by Ensembl gene ID, not gene name (i.e. symbol).
            # Gene name is encoded in feature_name, which is needed for gene search.
            feature_frame = adata.var.feature_name
            index = True
        else:
            # Default case
            feature_frame = adata.var.index
            index = False
        pd.DataFrame(feature_frame).to_csv(
            "h5ad_frag.features.processed.tsv.gz",
            sep="\t",
            index=index,
            header=False,
            compression="gzip",
        )
        pd.DataFrame(adata.obs.index).to_csv(
            "h5ad_frag.barcodes.processed.tsv.gz",
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )
        mtx_filename = "h5ad_frag.matrix.processed.mtx"
        MMFileFixedFormat().write(
            mtx_filename, a=scipy.sparse.csr_matrix(adata.X.T), precision=3
        )
        AnnDataIngestor.compress_file(mtx_filename)

    @staticmethod
    def check_ensembl_index(adata):
        """Check if an AnnData file is indexed on Ensembl gene IDs (e.g. ENSG00000243485) instead of gene symbols"""
        if adata.var.index.name == 'gene_ids':
            return True
        else:
            prefixes = list(set(gene_id[:3] for gene_id in adata.var_names))
            return len(prefixes) == 1 and prefixes[0] == 'ENS'

    @staticmethod
    def delocalize_extracted_files(
        file_path, study_file_id, accession, files_to_delocalize
    ):
        """Copy extracted files to study bucket"""

        for file in files_to_delocalize:
            IngestFiles.delocalize_file(
                file_path,
                file,
                f"_scp_internal/anndata_ingest/{accession}_{study_file_id}/{file}",
            )
