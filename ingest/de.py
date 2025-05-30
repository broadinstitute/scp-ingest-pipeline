import logging
import numpy as np
import pandas as pd
import scanpy as sc
import re
import glob
from anndata import AnnData

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
        annotation_name,
        **kwargs,
    ):
        DifferentialExpression.de_logger.info(
            "Initializing DifferentialExpression instance"
        )
        self.cluster = cluster
        self.metadata = cell_metadata
        self.annotation = annotation_name
        self.matrix_file_path = matrix_file_path
        self.matrix_file_type = matrix_file_type
        self.kwargs = kwargs
        self.accession = self.kwargs["study_accession"]
        self.annot_scope = self.kwargs["annotation_scope"]
        self.method = self.kwargs["method"]

        if matrix_file_type == "mtx":
            self.genes_path = self.kwargs["gene_file"]
            self.barcodes_path = self.kwargs["barcode_file"]

    @staticmethod
    def assess_annotation(annotation, metadata, extra_params):
        """Check that annotation for DE is not of TYPE numeric"""
        dtype_info = dict(zip(metadata.headers, metadata.annot_types))
        annotation_info = dtype_info.get(annotation, None)
        if annotation_info == "numeric":
            msg = f"DE analysis infeasible for numeric annotation \"{annotation}\"."
            print(msg)
            log_exception(
                DifferentialExpression.dev_logger, DifferentialExpression.de_logger, msg
            )
            raise TypeError(msg)
        elif annotation_info is None:
            msg = f"Provided annotation \"{annotation}\" not found in metadata file."
            print(msg)
            log_exception(
                DifferentialExpression.dev_logger, DifferentialExpression.de_logger, msg
            )
            raise KeyError(msg)
        elif annotation_info == "group":
            de_type = extra_params.get("de_type")
            group1 = extra_params.get("group1")
            group2 = extra_params.get("group2")
            if de_type == "pairwise":
                metadata.preprocess(False)
                values = metadata.file[annotation]['group'].unique()
                if group1 and group2:
                    if group1 not in values or group2 not in values:
                        msg = (
                            f'Provided annotation value(s) group1, \"{group1}\", '
                            f'or group2, \"{group2}\", not found in metadata file'
                        )
                        log_exception(
                            DifferentialExpression.dev_logger,
                            DifferentialExpression.de_logger,
                            msg,
                        )
                        raise ValueError(msg)
                else:
                    msg = (
                        f'Provided annotation value(s) group1, \"{group1}\", or '
                        f'group2, \"{group2}\", were false-y'
                    )
                    log_exception(
                        DifferentialExpression.dev_logger,
                        DifferentialExpression.de_logger,
                        msg,
                    )
                    raise KeyError(msg)
            # de_type is "rest" by default, DE annotation TYPE expected to be group
            else:
                # assessment does not raise error
                return None
        else:
            msg = f"Error: \"{annotation}\" has unexpected type \"{annotation_info}\"."
            print(msg)
            log_exception(
                DifferentialExpression.dev_logger, DifferentialExpression.de_logger, msg
            )
            raise ValueError(msg)

    @staticmethod
    def get_cluster_cells(cluster_cells):
        """ID cells in cluster file"""
        # cluster_cells.tolist() yields a list of lists that needs to be flattened
        # using extend converts a single-value list to a plain value
        cluster_cell_values = cluster_cells.tolist()
        cluster_cell_list = []
        for value in cluster_cell_values:
            cluster_cell_list.extend(value)
        return cluster_cell_list

    @staticmethod
    def determine_dtypes(headers, annot_types):
        """use SCP TYPE data to coerce data to proper dtypes:
        numeric-like group annotations
        missing values in group annotations (to avoid NaN)
        """
        dtype_info = dict(zip(headers, annot_types))
        dtypes = {}
        for header, type_info in dtype_info.items():
            if type_info == "group":
                # using dtype str to avoid StringArray objects when "string" dtype used
                # StringArray is a Experimental extension array, behavior is not yet set
                dtypes[header] = str
        return dtypes

    @staticmethod
    def process_annots(metadata_file_path, allowed_file_types, headers, dtypes):
        """using SCP metadata header lines create pandas dataframe where
        numeric-seeming group annotations are properly set to dtype of str
        and NaN in group columns are converted to "__Unspecified__"
        """
        annot_redux = IngestFiles(metadata_file_path, allowed_file_types)
        annot_file_type = annot_redux.get_file_type(metadata_file_path)[0]
        annot_file_handle, local_file_path = IngestFiles.resolve_path(
            annot_redux, metadata_file_path
        )
        annots = annot_redux.open_pandas(
            local_file_path,
            annot_file_type,
            open_file_object=annot_file_handle,
            names=headers,
            skiprows=2,
            index_col=0,
            dtype=dtypes,
            keep_default_na=False,
            na_values=[""],
        )
        group_annots = [k for k, v in dtypes.items() if v == str]
        # Where group metadata is missing values (eg. optional or nonconventional metadata)
        # replace NaN with the string '__Unspecified__'
        # intent is to reflect the '--Unspecified--' annotation label added to
        # SCP plot legends in scatter plot visualizations for unannotated points
        annots[group_annots] = annots[group_annots].fillna('__Unspecified__')
        return annots

    @staticmethod
    def subset_annots(metadata, de_cells):
        """subset metadata based on cells in cluster"""
        DifferentialExpression.de_logger.info(
            "subsetting metadata on cells in clustering"
        )
        dtypes = DifferentialExpression.determine_dtypes(
            metadata.headers, metadata.annot_types
        )
        orig_annots = DifferentialExpression.process_annots(
            metadata.file_path, metadata.ALLOWED_FILE_TYPES, metadata.headers, dtypes
        )
        cluster_annots = orig_annots[orig_annots.index.isin(de_cells)]
        return cluster_annots

    @staticmethod
    def order_annots(metadata, adata_cells):
        """order metadata based on cells order in matrix"""
        matrix_cell_order = adata_cells.tolist()
        return metadata.reindex(matrix_cell_order)

    @staticmethod
    def subset_adata(adata, de_cells):
        """subset adata object based on cells in cluster"""
        DifferentialExpression.de_logger.info(
            "subsetting matrix on cells in clustering"
        )
        matrix_subset_list = np.in1d(adata.obs_names, de_cells)
        adata = adata[matrix_subset_list].copy()
        return adata

    def execute_de(self):
        DifferentialExpression.de_logger.info(
            f'dev_info: Starting DE for {self.accession}'
        )
        try:
            # only used in output filename, replacing non-alphanumeric with underscores
            # except '+' replaced with 'pos'
            self.cluster_name = DifferentialExpression.sanitize_string(
                self.kwargs["name"]
            )
            if self.matrix_file_type in ["dense", "mtx", "h5ad"]:
                DifferentialExpression.de_logger.info(
                    f"preparing DE on {self.matrix_file_type} matrix"
                )
                self.run_scanpy_de(
                    self.cluster,
                    self.metadata,
                    self.matrix_file_path,
                    self.matrix_file_type,
                    self.annotation,
                    self.cluster_name,
                    self.kwargs,
                )
            else:
                msg = f"Submitted matrix_file_type should be \"dense\", \"mtx\" or \"h5ad\" not \"{self.matrix_file_type}\""
                print(msg)
                log_exception(
                    DifferentialExpression.dev_logger,
                    DifferentialExpression.de_logger,
                    msg,
                )
                raise ValueError(msg)
        except (TypeError, KeyError, ValueError) as e:
            raise ValueError(e)
        return

    @staticmethod
    def get_genes(genes_path):
        """Genes file can have one or two columns of gene information
        Preferentially use gene names from second column.
        If duplicate gene names, check that 1st plus 2nd column provides uniqueness
        If unique when joined, join columns with pipe (|) for use as DE input
        """
        genes_object = IngestFiles(genes_path, None)
        local_genes_path = genes_object.resolve_path(genes_path)[1]

        genes_df = pd.read_csv(local_genes_path, sep="\t", header=None)
        if len(genes_df.columns) > 1:
            # if genes are not unique, try combining with gene_id (SCP-4283)
            # print so we're aware of dups during dev testing
            if genes_df[1].count() != genes_df[1].nunique():
                warning = (
                    "dev_info: Features file contains duplicate identifiers in column 2"
                )
                print(warning)
                genes_df['new_id'] = genes_df[[1, 0]].agg('|'.join, axis=1)
                if genes_df['new_id'].count() != genes_df['new_id'].nunique():
                    msg = "Duplicates in features file even after joining gene_id and gene_name"
                    log_exception(
                        DifferentialExpression.dev_logger,
                        DifferentialExpression.de_logger,
                        msg,
                    )
                    raise ValueError(msg)
                else:
                    return genes_df['new_id'].tolist()
            else:
                return genes_df[1].tolist()
        else:
            if genes_df[0].count() != genes_df[0].nunique():
                warning = (
                    "dev_info: Features file contains duplicate identifiers in column 1"
                )
                print(warning)
            return genes_df[0].tolist()

    @staticmethod
    def get_barcodes(barcodes_path):
        """Extract barcodes from file for mtx reconstitution"""
        barcodes_ingest_file = IngestFiles(barcodes_path, "text/tab-separated-values")
        barcodes_file = barcodes_ingest_file.resolve_path(barcodes_path)[0]
        barcodes = [c.strip().strip('"') for c in barcodes_file.readlines()]
        return barcodes

    @staticmethod
    def adata_from_mtx(matrix_file_path, genes_path, barcodes_path):
        """reconstitute AnnData object from matrix, genes, barcodes files"""
        # process smaller files before reading larger matrix file
        barcodes = DifferentialExpression.get_barcodes(barcodes_path)
        features = DifferentialExpression.get_genes(genes_path)
        matrix_object = IngestFiles(matrix_file_path, None)
        local_file_path = matrix_object.resolve_path(matrix_file_path)[1]
        adata = sc.read_mtx(local_file_path)
        # For AnnData, obs are cells and vars are genes
        # BUT transpose needed for both dense and sparse
        # so transpose step is after this data object composition step
        # therefore the assignements below are the reverse of expected
        adata.var_names = barcodes
        adata.obs_names = features
        return adata

    @staticmethod
    def remove_single_sample_data(adata, annotation):
        """identify and remove cells that would constitute an annotation label
        that has data with only a single sample
        """
        counts = adata.obs[annotation].value_counts(dropna=False)
        for label, count in counts.iteritems():
            if count == 1:
                adata = adata[adata.obs[annotation] != label]
        return adata

    @staticmethod
    def delimiter_in_gene_name(rank):
        """Check if pipe delimiter occurs in "names" column"""
        return rank['names'].str.contains('|', regex=False).any()

    @staticmethod
    def extract_gene_id_for_out_file(rank):
        """Separate out gene name from gene ID"""
        rank['feature_id'] = rank['names'].str.split('|').str[1]
        rank['names'] = rank['names'].str.split('|').str[0]

        return rank

    @staticmethod
    def all_match_ensembl_id_regex(adata, col_name):
        regex = r'ENS[A-Z]{0,4}\d{11}'
        return adata[col_name].astype(str).str.match(regex).all()

    @staticmethod
    def write_de_result(adata, group, annotation, rank_key, cluster_name, extra_params):
        de_type = extra_params.get("de_type")
        method = extra_params.get("method")
        annot_scope = extra_params.get("annotation_scope")
        clean_annotation = DifferentialExpression.sanitize_string(annotation)
        rank = sc.get.rank_genes_groups_df(adata, key=rank_key, group=group)
        if DifferentialExpression.all_match_ensembl_id_regex(rank, "names"):
            feature_name_rank = rank.merge(
                adata.var['feature_name'], left_on="names", right_index=True, how='left'
            )
            feature_name_rank.rename(columns={'names': 'feature_id'}, inplace=True)
            feature_name_rank.rename(columns={'feature_name': 'names'}, inplace=True)
            new_column_order = (
                ["names"]
                + [
                    col
                    for col in feature_name_rank.columns
                    if col not in {"names", "feature_id"}
                ]
                + ["feature_id"]
            )
            feature_name_rank = feature_name_rank[new_column_order]

            rank = feature_name_rank
        if DifferentialExpression.delimiter_in_gene_name(rank):
            DifferentialExpression.extract_gene_id_for_out_file(rank)
        if de_type == "rest":
            clean_group = DifferentialExpression.sanitize_string(group)
            out_file = f'{cluster_name}--{clean_annotation}--{clean_group}--{annot_scope}--{method}.tsv'
            DifferentialExpression.de_logger.info(
                f"Writing DE output for {clean_group} vs rest"
            )
        elif de_type == "pairwise":
            # rank_genes_groups accepts a list. For SCP pairwise, should be a list with one item
            # converting list to string for incorporation into result filename
            group1 = DifferentialExpression.sanitize_string(
                ''.join(extra_params["group1"])
            )
            group2 = DifferentialExpression.sanitize_string(extra_params["group2"])
            out_file = f'{cluster_name}--{clean_annotation}--{group1}--{group2}--{annot_scope}--{method}.tsv'
            DifferentialExpression.de_logger.info(
                f"Writing DE output for {group1} vs {group2} pairwise"
            )
        else:
            msg = f'Unknown de_type, {de_type}'
            print(msg)
            log_exception(
                DifferentialExpression.dev_logger, DifferentialExpression.de_logger, msg
            )
            raise ValueError(msg)
        # Round numbers to 4 significant digits while respecting fixed point
        # and scientific notation (note: trailing zeros are removed)
        rank.to_csv(out_file, sep='\t', float_format='%.4g')

    @staticmethod
    def run_scanpy_de(
        cluster,
        metadata,
        matrix_file_path,
        matrix_file_type,
        annotation,
        cluster_name,
        extra_params,
    ):
        method = extra_params.get("method")
        de_type = extra_params.get("de_type")
        raw_location = extra_params.get("raw_location")

        try:
            DifferentialExpression.assess_annotation(annotation, metadata, extra_params)
        except (TypeError, KeyError, ValueError) as e:
            raise ValueError(e)

        de_cells = DifferentialExpression.get_cluster_cells(cluster.file['NAME'].values)
        de_annots = DifferentialExpression.subset_annots(metadata, de_cells)

        if matrix_file_type == "dense":
            # will need try/except (SCP-4205)
            matrix_object = IngestFiles(matrix_file_path, None)
            local_file_path = matrix_object.resolve_path(matrix_file_path)[1]
            adata = sc.read(local_file_path)
        elif matrix_file_type == "h5ad":
            matrix_object = IngestFiles(matrix_file_path, None)
            local_file_path = matrix_object.resolve_path(matrix_file_path)[1]
            # orig_adata = matrix_object.open_anndata(local_file_path)
            orig_adata = sc.read_h5ad(local_file_path)
        else:
            # MTX reconstitution UNTESTED (SCP-4203)
            # will want try/except here to catch failed data object composition
            genes_path = extra_params.get("gene_file")
            barcodes_path = extra_params.get("barcode_file")
            adata = DifferentialExpression.adata_from_mtx(
                matrix_file_path, genes_path, barcodes_path
            )

        if matrix_file_type == "h5ad":
            if raw_location == ".raw":
                if orig_adata.raw is not None:
                    DifferentialExpression.de_logger.info(
                        f"Performing DE on {raw_location} data"
                    )
                    adata = AnnData(
                        # using .copy() for the AnnData components is good practice
                        # but we won't be using orig_adata for analyses
                        # choosing to avoid .copy() for memory efficiency
                        X=orig_adata.raw.X,
                        obs=orig_adata.obs,
                        var=orig_adata.var,
                    )
                else:
                    msg = f'{matrix_file_path} does not have a .raw attribute'
                    print(msg)
                    log_exception(
                        DifferentialExpression.dev_logger,
                        DifferentialExpression.de_logger,
                        msg,
                    )
                    raise ValueError(msg)
            else:
                if raw_location in orig_adata.layers.keys():
                    DifferentialExpression.de_logger.info(
                        f"Performing DE on adata.layers['{raw_location}'] data"
                    )
                    adata = AnnData(
                        # using .copy() for the AnnData components is good practice
                        # but we won't be using orig_adata for analyses
                        # choosing to avoid .copy() for memory efficiency
                        X=orig_adata.layers[raw_location],
                        obs=orig_adata.obs,
                        var=orig_adata.var,
                    )
                else:
                    msg = f'{matrix_file_path} does not have adata.layers["{raw_location}"]'
                    print(msg)
                    log_exception(
                        DifferentialExpression.dev_logger,
                        DifferentialExpression.de_logger,
                        msg,
                    )
                    raise ValueError(msg)
        # AnnData expects gene x cell so dense and mtx matrices require transposition
        else:
            adata = adata.transpose()

        adata = DifferentialExpression.subset_adata(adata, de_cells)

        # h5ad inputs will already have obs data, only non-h5ad need this step
        if not matrix_file_type == "h5ad":
            adata.obs = DifferentialExpression.order_annots(de_annots, adata.obs_names)

        adata = DifferentialExpression.remove_single_sample_data(adata, annotation)

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        DifferentialExpression.de_logger.info("calculating DE")
        rank_key = f"rank.{annotation}.{method}"
        if de_type == "rest":
            try:
                sc.tl.rank_genes_groups(
                    adata,
                    annotation,
                    key_added=rank_key,
                    method=method,
                    pts=True,
                )
            except KeyError:
                msg = f"Missing expected annotation in metadata: {annotation}, unable to calculate DE."
                print(msg)
                log_exception(
                    DifferentialExpression.dev_logger,
                    DifferentialExpression.de_logger,
                    msg,
                )
                raise KeyError(msg)
            # ToDo - detection and handling of annotations with only one sample (SCP-4282)
            except ValueError as e:
                print(e)
                log_exception(
                    DifferentialExpression.dev_logger,
                    DifferentialExpression.de_logger,
                    e,
                )
                raise KeyError(e)

            DifferentialExpression.de_logger.info("Gathering DE annotation labels")
            groups = np.unique(adata.obs[annotation]).tolist()
            for group in groups:
                DifferentialExpression.write_de_result(
                    adata, group, annotation, rank_key, cluster_name, extra_params
                )

        elif de_type == 'pairwise':
            group1 = [extra_params["group1"]]
            group2 = extra_params["group2"]
            try:
                sc.tl.rank_genes_groups(
                    adata,
                    groupby=annotation,
                    groups=group1,
                    reference=group2,
                    key_added=rank_key,
                    method=method,
                    pts=True,
                )
            except KeyError:
                msg = f"Missing expected annotation in metadata: {annotation}, unable to calculate DE."
                print(msg)
                log_exception(
                    DifferentialExpression.dev_logger,
                    DifferentialExpression.de_logger,
                    msg,
                )
                raise KeyError(msg)
            # ToDo - detection and handling of annotations with only one sample (SCP-4282)
            except ValueError as e:
                print(e)
                log_exception(
                    DifferentialExpression.dev_logger,
                    DifferentialExpression.de_logger,
                    e,
                )
                raise KeyError(e)
            stringified_group1 = ''.join(extra_params["group1"])
            DifferentialExpression.write_de_result(
                adata,
                stringified_group1,
                annotation,
                rank_key,
                cluster_name,
                extra_params,
            )

        DifferentialExpression.de_logger.info("DE processing complete")

    @staticmethod
    def sanitize_string(input_string):
        """
        Replace '+' with 'pos', then replace non-alphanumerics with underscore
        this allows distinct sanitization for "CD16+ monocyte" vs "CD16- monocyte"
        """
        plus_converted_string = re.sub('\+', 'pos', input_string)
        return re.sub(r'\W', '_', plus_converted_string)

    @staticmethod
    def string_for_output_match(arguments):
        cleaned_cluster_name = DifferentialExpression.sanitize_string(
            arguments["cluster_name"]
        )
        cleaned_annotation_name = DifferentialExpression.sanitize_string(
            arguments["annotation_name"]
        )
        files_to_match = f"{cleaned_cluster_name}--{cleaned_annotation_name}*.tsv"
        return files_to_match

    @staticmethod
    def delocalize_de_files(destination_file_path, study_file_id, files_to_match):
        """Copy DE output files to study bucket"""

        files = glob.glob(files_to_match)
        for file in files:
            IngestFiles.delocalize_file(
                destination_file_path,
                file,
                f"_scp_internal/differential_expression/{file}",
            )
