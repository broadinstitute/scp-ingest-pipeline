"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

PREREQUISITES
Must have python 3.6 or higher.
"""
from typing import List  # noqa: F401

from gene_data_model import Gene
from ingest_files import IngestFiles


class Dense(IngestFiles):
    def __init__(self, file_path, file_id, study_accession, file_params):
        self.ALLOWED_FILE_TYPES = [
            "text/csv",
            "text/plain",
            "text/tab-separated-values",
        ]
        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
        self.file_id = file_id
        self.study_accession = study_accession
        self.cell_names = self.get_next_line(increase_line_count=False)[1:]
        self.file_params = file_params

    def transform_expression_data_by_gene(self, expression_scores: List[str]) -> Gene:
        """Transforms dense matrix into firestore data model for genes.

        Args:
            lines : List[str]
                Lines from dense matrix file

        Returns:
                transformed_data : List[Gene]
                A list of Gene objects
        """
        gene_model = Gene(
            name=expression_scores[0],
            source_file_type="Dense",
            expression_scores=expression_scores[1:],
            cell_names=self.cell_names,
            study_accession=self.study_accession,
            file_id=self.file_id,
            **self.file_params
        )
        print(gene_model.top_level_doc)
        return gene_model

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.file.close()
