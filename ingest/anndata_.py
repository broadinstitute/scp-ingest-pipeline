import pandas as pd  # NOqa: F821
import os
import datetime
import scanpy as sc

try:
    from ingest_files import IngestFiles
    from expression_files.expression_files import GeneExpression
    from monitor import log_exception
    from validation.validate_metadata import list_duplicates
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import IngestFiles
    from .expression_files.expression_files import GeneExpression
    from .monitor import log_exception
    from .validation.validate_metadata import list_duplicates


class AnnDataIngestor(GeneExpression, IngestFiles):
    ALLOWED_FILE_TYPES = ['application/x-hdf5']

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        GeneExpression.__init__(self, file_path, study_file_id, study_id)
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

    @staticmethod
    def generate_cluster_header(adata, clustering_name):
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
    def clusterings_to_delocalize(arguments):
        # ToDo - check if names using obsm_keys need sanitization
        cluster_file_names = []
        for name in arguments["obsm_keys"]:
            cluster_file_names.append(AnnDataIngestor.set_clustering_filename(name))
        return cluster_file_names

    @staticmethod
    def delocalize_extracted_files(file_path, study_file_id, files_to_delocalize):
        """Copy extracted files to study bucket"""

        for file in files_to_delocalize:
            IngestFiles.delocalize_file(
                study_file_id,
                None,
                file_path,
                file,
                f"_scp_internal/anndata_ingest/{study_file_id}/{file}",
            )

    @staticmethod
    def check_valid(adata):
        error_messages = []

        try:
            AnnDataIngestor.check_names_unique(adata.var_names, "Feature")
        except ValueError as v:
            error_messages.append(str(v))
        try:
            AnnDataIngestor.check_names_unique(adata.obs_names, "Obs")
        except ValueError as v:
            error_messages.append(str(v))
        if len(error_messages) > 0:
            raise ValueError("; ".join(error_messages))

        return True

    def process_matrix(self):
        """Perform matrix processing"""
        if self.check_valid(self.adata):
            self.transform()

    @staticmethod
    def check_names_unique(names, name_type):
        """Return True if names are unique, else false
        Expected name_types: ["Feature", "Obs"]
        """
        # check feature_name and obs names, feature_id logic not included
        # TODO (SCP-5105) non-happy path - add feature_id assessment
        if len(names) == len(names.unique()):
            return True
        else:
            dups = list_duplicates(names)
            features_for_msg = 2
            end = features_for_msg if len(dups) > features_for_msg else len(dups)
            dup_list = dups[:end]
            dup_string = " ".join(dup_list)

            msg = (
                f"{name_type} names must be unique within a file. "
                f"{len(dups)} duplicates found, including: {dup_string}"
            )
            GeneExpression.log_for_mixpanel(
                "error", "content:duplicate:values-within-file", msg
            )
            raise ValueError(msg)

    def transform(self):
        """Transforms matrix into gene data model."""
        # initialize settings for mock data loads in tests
        self.test_models = None
        self.models_processed = 0

        # derive file name from file path
        file_name = os.path.basename(self.file_path)
        start_time = datetime.datetime.now()
        GeneExpression.dev_logger.info("Starting run at " + str(start_time))
        num_processed = 0
        gene_models = []
        data_arrays = []
        for all_cell_model in GeneExpression.create_data_arrays(
            name=f"{file_name} Cells",
            array_type="cells",
            values=self.adata.obs.index.tolist(),
            linear_data_type="Study",
            linear_data_id=self.study_file_id,
            **self.data_array_kwargs,
        ):
            data_arrays.append(all_cell_model)

        # ASSUMPTION all_cell_model same for raw_count and processed_expression
        # TODO (SCP-5103): if raw counts is indicated check that .raw slot is populated

        # Iterate over feature names (for happy path)
        for feature in self.adata.var_names.tolist():
            print(f"processing feature: {feature}")
            feature_expression_series = sc.get.obs_df(self.adata, keys=feature)
            if feature_expression_series.hasnans:
                msg = (
                    f'Expected numeric expression score - '
                    f'expression data has NaN values for feature "{feature}"'
                )
                GeneExpression.log_for_mixpanel(
                    "error", "content:type:not-numeric", msg
                )
                raise ValueError(msg)
            # capture sparse (only non zero values and their cell IDs)
            # check mtx.py for all zero gene handling
            filtered_expression_series = feature_expression_series[
                feature_expression_series.values > 0
            ]

            exp_cells = filtered_expression_series.index.tolist()

            untrimmed_exp_scores = filtered_expression_series.values.tolist()

            # trim expression data to three significant digits
            exp_scores = [round(float(value), 3) for value in untrimmed_exp_scores]
            # TODO (SCP-5105) for None value below, replace with feature ID (string)
            data_arrays, gene_models, num_processed = self.create_models(
                exp_cells,
                exp_scores,
                feature,
                None,
                gene_models,
                data_arrays,
                num_processed,
                False,
            )
        # Load any remaining models. This is necessary because the amount of
        # models may be less than the batch size.
        if len(gene_models) > 0 or len(data_arrays) > 0:
            self.create_models(
                [], [], None, None, gene_models, data_arrays, num_processed, True
            )
