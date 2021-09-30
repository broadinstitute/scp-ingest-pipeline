"""
Module for ingesting MTX files

DESCRIPTION
This module provides extract and transforms function for an MTX file bundle. An
MTX file bundle consists of A) an .mtx file in Matrix
Market matrix coordinate format, B) a genes.tsv file, and a barcodes.tsv file.
These are commonly provided from 10x Genomics v2.

Unsorted MTX files can be parsed. Files are considered unsorted when file is not
sorted by gene. In this event, the file will automatically be sorted by gene.
Version of sort: (GNU coreutils) 8.25

"""

from typing import Dict, Generator, List, Tuple, IO  # noqa: F401

import os

import datetime
import subprocess

try:
    from expression_files import GeneExpression
    from ingest_files import IngestFiles
    from monitoring.mixpanel_log import custom_metric
    import config
except ImportError:
    # Used when importing as external package, e.g. imports in
    # single_cell_portal code
    from .expression_files import GeneExpression
    from ..ingest_files import IngestFiles
    from ..monitoring.mixpanel_log import custom_metric
    from .. import config


class MTXIngestor(GeneExpression, IngestFiles):
    ALLOWED_FILE_TYPES = ["text/tab-separated-values"]

    @staticmethod
    def matches_file_type(file_type):
        return "mtx" == file_type

    def __init__(self, mtx_path: str, study_file_id: str, study_id: str, **kwargs):
        GeneExpression.__init__(self, mtx_path, study_file_id, study_id)
        IngestFiles.__init__(self, mtx_path, self.ALLOWED_FILE_TYPES)
        self.mtx_file, self.mtx_path = self.resolve_path(mtx_path)

        genes_path = kwargs.pop("gene_file")
        genes_ingest_file = IngestFiles(genes_path, self.ALLOWED_FILE_TYPES)
        self.genes_file = genes_ingest_file.resolve_path(genes_path)[0]

        barcodes_path = kwargs.pop("barcode_file")
        barcodes_ingest_file = IngestFiles(barcodes_path, self.ALLOWED_FILE_TYPES)
        self.barcodes_file = barcodes_ingest_file.resolve_path(barcodes_path)[0]

        # A list ['N', 'K', 'M'] that represents a gene-barcode matrix where N
        # is the gene index, M is the barcode index, and K is the expression
        # score for the given gene index
        self.mtx_dimensions: List[int] = MTXIngestor.get_mtx_dimensions(self.mtx_file)

    @staticmethod
    def check_valid(
        barcodes: List[str], genes: List[str], mtx_dimensions, query_params
    ):
        error_messages = []

        try:
            MTXIngestor.check_bundle(barcodes, genes, mtx_dimensions)
        except ValueError as v:
            error_messages.append(str(v))
        try:
            MTXIngestor.check_duplicates(genes, "gene")
        except ValueError as v:
            error_messages.append(str(v))
        try:
            MTXIngestor.check_duplicates(barcodes, "barcode")
        except ValueError as v:
            error_messages.append(str(v))
        try:
            GeneExpression.check_unique_cells(barcodes, *query_params)
        except ValueError as v:
            error_messages.append(str(v))

        if len(error_messages) > 0:
            raise ValueError("; ".join(error_messages))
        return True

    @staticmethod
    def get_gene_expression_data(file_path: str, file_handler):
        """Pipes expression data out for zipped and unzipped files"""
        std_in = file_path
        file_size = os.path.getsize(file_path)
        if file_size == 0:
            raise ValueError(f"{file_path} is empty: " + str(file_size))

        if IngestFiles.is_gzipped(file_path):
            base = os.path.basename(file_path)
            base = os.path.splitext(base)[0]
            new_file_name = f"{base}_unzipped.mtx"
            # If uncompressed file does not exist, create it.
            if not os.path.isfile(new_file_name):
                # Uncompress file and write to new_file_name.
                with open(new_file_name, "w+") as f:
                    subprocess.run(
                        ["gunzip", "-c", file_path], stdout=f, universal_newlines=True
                    )
            # Set standard in to new uncompressed file
            std_in = new_file_name

        # Grab gene expression data which starts at 'n', or start_idx, lines from top of file (-n +{start_idx}).
        # The header and mtx dimension aren't included or else the file would always be considered unsorted due to the mtx
        # dimensions always being larger than the first row of data.
        start_idx: int = MTXIngestor.get_data_start_line_number(file_handler)
        p1_cmd = ["tail", "-n", f"+{start_idx}", std_in]
        p1 = subprocess.Popen(p1_cmd, stdout=subprocess.PIPE)
        return p1

    @staticmethod
    def is_sorted(file_path: str, file_handler: IO):
        """Checks if a file is sorted by gene index"""
        p1 = MTXIngestor.get_gene_expression_data(file_path, file_handler)
        try:
            # Check that the input file is sorted (-c). Check sorting by first and only first column (-k 1,1,). This value is
            # numeric(-n). Use stable sort (--stable) so that columns are only compared by the first column.
            # Without this argument line '1 3 4' and ' 1 5 4' would be considered unsorted.
            p2 = subprocess.run(
                ["sort", "-c", "--stable", "-n", "-k", "1,1"],
                stdin=p1.stdout,
                capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            raise e
        else:
            # When a file is unsorted, or 'disordered', the output includes the string 'disorder' which indicates the
            # first offending row that caused the 'disorder'
            if "disorder" in str(p2.stderr):
                return False
            else:
                return True

    @staticmethod
    def check_bundle(barcodes, genes, mtx_dimensions):
        """Confirms barcode and gene files have expected length of values"""
        expected_genes = mtx_dimensions[0]
        actual_genes = len(genes)

        expected_barcodes = mtx_dimensions[1]
        actual_barcodes = len(barcodes)

        if (actual_barcodes == expected_barcodes) and (actual_genes == expected_genes):

            return True
        else:
            msg = (
                f"Expected {expected_barcodes} cells and {expected_genes} genes. "
                f"Got {actual_barcodes} cells and {actual_genes} genes."
            )
            raise ValueError(msg)

    @staticmethod
    def check_duplicates(names: List, file_type: str):
        """Checks for duplicate values.
        Barcode and gene files cannot contain duplicate values within the file

        Parameters
        ----------
        names - Gene or cell values
        file_type - Barcode or gene files. Used in error message
        """
        unique_names: List[str] = set(names)
        if len(names) > len(unique_names):
            amount_of_duplicates = abs(len(unique_names) - len(names))
            msg = (
                "Duplicate values are not allowed. "
                f"There are {amount_of_duplicates} duplicates "
                f"in the {file_type} file"
            )
            raise ValueError(msg)
        return True

    @staticmethod
    def get_data_start_line_number(file_handler: IO) -> int:
        """ Determines what line number data starts.

        Parameters:
        ___________
            file_handler (IO): File handler of MTX file that contains headers.

         Returns
         ----------
            count (int): Line number where data starts
        """
        # Move file pointer to top of file
        file_handler.seek(0, 0)
        for count, line in enumerate(file_handler, start=1):
            if not line.startswith("%"):
                try:
                    line_values = line.strip().split()
                    float(line_values[0])  # Determines if value is numeric
                    # First line w/o '%' is mtx dimension. So skip this line (+1)
                    return count + 1
                except ValueError:
                    raise ValueError(
                        "Only header, comment lines starting with '%', and numeric data allowed in MTX file."
                    )
                except IndexError:
                    raise IndexError("MTX file cannot start with a space")
        raise ValueError(
            "MTX file did not contain expression data. Please check formatting and contents of file."
        )

    @staticmethod
    def get_mtx_dimensions(file_handler) -> List:
        for line in file_handler:
            if not line.startswith("%"):
                mtx_dimensions: List[str] = line.strip().split()
                try:
                    # Convert values in mtx_dimensions to int
                    dimensions = list(map(int, mtx_dimensions))
                    return dimensions
                except Exception as e:
                    raise e
        raise ValueError("MTX file did not contain data")

    @staticmethod
    def get_features(feature_row: str):
        """Determines gene id and gene name from a given row:str in a feature
            file
        """
        feature_data = feature_row.split("\t")
        gene_id = feature_data[0]
        gene_name = feature_data[0]
        if len(feature_data) >= 2:
            # gene_name field is present
            gene_name = feature_data[1]
        return gene_id, gene_name

    @staticmethod
    @custom_metric(config.get_metric_properties, "perfTime:sort_mtx", {"sorted": False})
    def sort_mtx(file_path, mtx_file_handler: IO) -> str:
        """
        Sorts MTX file by gene. File header, dimensions, and comments are not included in sort.

         Parameters:
            mtx_file_handler (IO): File handler that points to top of MTX file
            file_path (str): Path path to MTX file

         Returns
         ----------
            new_file_path (str) : Full path of newly sorted MTX file. This file does not contain original headers
                or MTX dimensions.
        """
        file_name = os.path.basename(file_path)
        file_name = os.path.splitext(file_name)[0]
        new_file_name = f"{file_name}_sorted_MTX.mtx"

        p1 = MTXIngestor.get_gene_expression_data(file_path, mtx_file_handler)
        with open(new_file_name, "w+") as f:
            GeneExpression.dev_logger.info("Starting to sort")
            start_time = datetime.datetime.now()
            # Sort output of p1 ( all gene expression data) by first and only first column (-k 1,1,)
            # using 20G for the memory buffer (-S 20G) with a max maximum number of 320 temporary files (--batch-size=320)
            # that can be merged at once (instead of default 16). The data being sorted is numeric (-n).
            # Use compress program gzip to compress temporary files (--compress-program=gzip).
            subprocess.run(
                [
                    "sort",
                    "--compress-program=gzip",
                    "-S",
                    "20G",
                    "--batch-size=320",
                    "-n",
                    "-k",
                    "1,1",
                ],
                stdin=p1.stdout,
                stdout=f,
            )
        GeneExpression.dev_logger.info(
            f"Time to sort {str(datetime.datetime.now() - start_time)} "
        )
        GeneExpression.dev_logger.info("Finished sorting")
        new_file_path = f"{os.getcwd()}/{new_file_name}"
        return new_file_path

    def execute_ingest(self):
        """Parses MTX files"""
        self.extract_feature_barcode_matrices()
        MTXIngestor.check_valid(
            self.cells,
            self.genes,
            self.mtx_dimensions,
            query_params=(
                self.study_id,
                self.study_file_id,
                self.mongo_connection._client,
            ),
        )
        if not GeneExpression.is_raw_count_file(
            self.study_id, self.study_file_id, self.mongo_connection._client
        ):
            self.is_raw_count = False
            if not MTXIngestor.is_sorted(self.mtx_path, self.mtx_file):
                new_mtx_file_path = MTXIngestor.sort_mtx(self.mtx_path, self.mtx_file)
                # Reset mtx variables to newly sorted file
                self.mtx_file, self.mtx_path = self.resolve_path(new_mtx_file_path)
        else:
            # Cell names are the only data stored for raw counts.
            # Therefore, no need to determine if file is sorted or sort file.
            self.is_raw_count = True
        self.transform()

    def extract_feature_barcode_matrices(self):
        """
        Sets relevant iterables for the gene and barcode file of the MTX bundle
        """
        self.genes: List[str] = [
            g.strip().strip('"') for g in self.genes_file.readlines()
        ]
        self.cells: List[str] = [
            c.strip().strip('"') for c in self.barcodes_file.readlines()
        ]

    def transform(self):
        """Transform data into data models.
        Transform() is expecting the file handle to be at the first line of data.
        """
        num_processed = 0
        prev_idx = 0
        gene_models = []
        data_arrays = []
        exp_cells = []
        exp_scores = []
        visited_expression_indices = {}
        current_gene = None
        current_gene_id = None

        # All observed cells
        for data_array in GeneExpression.create_data_arrays(
            name=f"{self.cluster_name} Cells",
            array_type="cells",
            values=self.cells,
            linear_data_type="Study",
            linear_data_id=self.study_file_id,
            **self.data_array_kwargs,
        ):
            data_arrays.append(data_array)
        # Create models for non-raw count files
        if not self.is_raw_count:
            for row in self.mtx_file:
                raw_gene_idx, raw_barcode_idx, raw_exp_score = row.split()
                current_idx = int(raw_gene_idx)
                if current_idx != prev_idx:
                    if not current_idx > prev_idx:
                        raise ValueError("MTX file must be sorted")
                    GeneExpression.dev_logger.debug(
                        f"Processing {self.genes[prev_idx - 1]}"
                    )
                    visited_expression_indices[current_idx] = True
                    if prev_idx != 0:
                        # Expressed cells and scores are associated with prior gene
                        prev_gene_id, prev_gene = MTXIngestor.get_features(
                            self.genes[prev_idx - 1]
                        )
                        # If the previous gene exists, load its models
                        data_arrays, gene_models, num_processed = self.create_models(
                            exp_cells,
                            exp_scores,
                            prev_gene,
                            prev_gene_id,
                            gene_models,
                            data_arrays,
                            num_processed,
                            False,
                        )
                        exp_cells = []
                        exp_scores = []
                    prev_idx = current_idx
                exp_cell = self.cells[int(raw_barcode_idx) - 1]
                exp_score = round(float(raw_exp_score), 3)
                exp_cells.append(exp_cell)
                exp_scores.append(exp_score)

            # create gene entries for genes with no positive expression values
            for idx, gene in enumerate(self.genes):
                if not visited_expression_indices.get(idx + 1):
                    current_gene_id, current_gene = MTXIngestor.get_features(gene)

                    data_arrays, gene_models, num_processed = self.create_models(
                        [],
                        [],
                        current_gene,
                        current_gene_id,
                        gene_models,
                        data_arrays,
                        num_processed,
                        False,
                    )

            # Create data array for last row
            current_gene_id, current_gene = MTXIngestor.get_features(
                self.genes[prev_idx - 1]
            )
        self.create_models(
            exp_cells,
            exp_scores,
            current_gene,
            current_gene_id,
            gene_models,
            data_arrays,
            num_processed,
            True,
        )
