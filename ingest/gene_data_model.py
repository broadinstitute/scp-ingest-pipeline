"""
Creates a gene data model.

DESCRIPTION
This module currently takes in parameters from a file and creates a data model
for genes.

PREREQUISITES
Must have python 3.6 or higher.
"""
from ingest_files import IngestFiles
import abc
from dataclasses import dataclass
from mypy_extensions import TypedDict
from typing import List  # noqa: F401


class GeneExpression(IngestFiles):
    __metaclass__ = abc.ABCMeta
    # SUBCOLLECTION_NAME = "gene_expression"
    # COLLECTION_NAME = "genes"

    # This model pertains to columns from cell metadata files
    @dataclass
    class Model(TypedDict):
        name: str
        # downcase version of 'name'
        searchable_name: str
        study_file_id: str
        study_id: str
        gene_id: str = None

    def __init__(
        self,
        file_path: str,
        study_id: str,
        study_file_id: str,
        allowed_file_types: str,
        source_file_type: str = None,
        expression_scores: List = None,
        *,
        gene_id: str = None,
        taxon_name: str = None,
        taxon_common_name: str = None,
        ncbi_taxid: int = None,
        genome_assembly_accession: str = None,
        genome_annotation: str = None,
        check_for_zero_values: bool = True,
        open_as=None,
    ) -> None:
        IngestFiles.__init__(self, file_path, allowed_file_types, open_as=open_as)
        # self.preproccess()
        # self.name = name.replace('"', "")
        # self.gene_id = gene_id
        # self.source_file_type = (source_file_type,)

        # if check_for_zero_values:
        #     self.cell_names, self.expression_scores = self.set_expression_scores_and_cell_names(
        #         expression_scores, cell_names
        #     )
        # else:
        #     self.cell_names = cell_names
        #     self.expression_scores = expression_scores
        # # Subdocument that contains all cell names and expression scores for a
        # # given gene
        # self.subdocument = {
        #     "cell_names": self.cell_names,
        #     "expression_scores": self.expression_scores,
        #     "source_file_type": source_file_type,
        # }
        # # This is the top level document for the gene data model
        # self.top_level_doc = {
        #     "file_id": file_id,
        #     "searchable_name": name.lower(),
        #     "name": self.name,
        #     "gene_id": self.gene_id,
        #     "study_accession": study_accession,
        #     "taxon_name": taxon_name,
        #     "taxon_common_name": taxon_common_name,
        #     "ncbi_taxid": ncbi_taxid,
        #     "genome_assembly_accession": genome_assembly_accession,
        #     "genome_annotation": genome_annotation,
        # }

    def has_subcollection_data(self):
        return self.cell_names is not None

    @abc.abstractmethod
    def transform(self):
        """Returns Gene data model"""

    def preproccess(self):
        # Remove trailing white spaces, and quotes from column names
        self.file.rename(
            columns=lambda x: x.strip().replace('"', '').replace("'", ''), inplace=True
        )
