"""Creates a gene data model for Firestore
Creates a gene data model for Firestore

DESCRIPTION
This module currently takes in parameters from a file and creates a gene data
model for Firesotre.

PREREQUISITES
Must have python 3.6 or higher.
"""
from typing import *


class Gene:
    def __init__(self, name: str, source_file_name: str, source_file_type: str, *,
                 gene_id: str = '', study_accession: str = "", taxon_name: str = "",
                 taxon_common_name: str = "", ncbi_taxid: str = "", genome_assembly_accession: str = "",
                 genome_annotation: str = "", cell_names: List[str] = [], expression_scores: List[float] = []) -> None:

        self.name = name
        self.gene_id = gene_id
        self.study_accession = study_accession
        self.source_file_name = source_file_name
        self.source_file_type = source_file_type
        self.taxon_name = taxon_name
        self.taxon_common_name = taxon_common_name
        self.ncbi_taxid = ncbi_taxid
        self.genome_assembly_accession = genome_assembly_accession
        self.genome_annotation = genome_annotation
        self.cell_names = cell_names
        self.expression_scores = expression_scores

        self.gene = {'Gene': {
            # TODO:  This part will need to be changed when we've decided on how to
            # idenitify  Document IDs
            self.name: {
                'name': self.name,
                'gene_id': self.gene_id,
                'study_accession': self.study_accession,
                'taxon_name': self.taxon_name,
                'taxon_common_name': self.taxon_common_name,
                'ncbi_taxid': self.ncbi_taxid,
                'genome_assembly_accession': self.genome_assembly_accession,
                'genome_annotation': self.genome_annotation,
                'gene_expression': {
                    # TODO: This part will need to be changed when we've decided
                    # on how to idenitify Document IDs.
                    # Another level needs to be added
                    'cell_names': self.cell_names,
                    'expression_scores': self.expression_scores,
                    'source_file_name': self.source_file_name,
                    'source_file_type': self.source_file_type,
                }
            }
        }}

    def get_expression_scores(self):
        """Get expression scores

        Args:
            None

        Returns:
            expression scores as list
        """
        return self.gene.get('Gene', {}).get(self.name, {}).get(
            'gene_expression')['expression_scores']
