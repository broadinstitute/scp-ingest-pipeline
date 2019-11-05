"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

PREREQUISITES
Must have python 3.6 or higher.
"""
from typing import List  # noqa: F401
from gene_data_model import GeneExpression


class Dense(GeneExpression):
    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        self.ALLOWED_FILE_TYPES = [
            "text/csv",
            "text/plain",
            "text/tab-separated-values",
        ]
        GeneExpression.__init__(
            self,
            file_path,
            study_file_id,
            study_id,
            self.ALLOWED_FILE_TYPES,
            open_as='dataframe',
        )
        self.study_file_id = study_file_id
        self.study_id = study_id
        # Remove from dictionary any keys that have value=None
        self.matrix_params = kwargs
        self.preproccess()

    def transform_expression_data_by_gene(self):
        """Transforms dense matrix into firestore data model for genes.

        Args:
            lines : List[str]
                Lines from dense matrix file

        Returns:
                transformed_data : List[Gene]
                A list of Gene objects
        """
        for gene in self.file['GENE']:
            gene = gene.replace('"', '').replace("'", '')
            yield self.Model(
                {
                    'name': gene,
                    'searchable_name': gene.lower(),
                    'study_file_id': self.study_file_id,
                    'study_id': self.study_id,
                    'gene_id': self.matrix_params['gene_id']
                    if 'gene_id' in self.matrix_params
                    else None,
                }
            )

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.file.close()
