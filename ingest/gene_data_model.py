"""
Creates a gene data model.

DESCRIPTION
This module currently takes in parameters from a file and creates a data model
for genes.

PREREQUISITES
Must have python 3.6 or higher.
"""

import sys
from itertools import islice
from typing import *

DOCUMENT_LIMIT_BYTES = 1048487


class Gene:
    def __init__(self, name: str, source_file_name: str, source_file_type: str, *,
                 gene_id: str = '', study_accession: str = '', taxon_name: str = '',
                 taxon_common_name: str = '', ncbi_taxid: str = '', genome_assembly_accession: str = '',
                 genome_annotation: str = '', cell_names: List[str] = [], expression_scores: List[float] = []) -> None:

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
        self.cell_names, self.expression_scores = self.set_expression_scores_and_cell_names(
            expression_scores, cell_names)

        self.gene_expression_subdocument = {'cell_names': self.cell_names,
                                            'expression_scores':  self.expression_scores,
                                            'source_file_name': source_file_name,
                                            'source_file_type': source_file_type,
                                            }
        self.gene = {
            'name': self.name,
            'gene_id': self.gene_id,
            'study_accession': self.study_accession,
            'taxon_name': self.taxon_name,
            'taxon_common_name': self.taxon_common_name,
            'ncbi_taxid': self.ncbi_taxid,
            'genome_assembly_accession': self.genome_assembly_accession,
            'genome_annotation': self.genome_annotation,
        }
        self.subcollection_name = 'gene_expression'

    def get_collection_name(self):
        print(f'Gene is {self.name}')
        if self.expression_scores is not None:
            print(f'expression score length is {len(self.expression_scores)}')
        return 'Gene'

    def get_subcollection_name(self):
        return 'gene_expression'

    def get_document(self):
        return self.gene

    def has_subcollection_data(self):
        return self.cell_names != None

    def get_subcollection(self):
        return self.gene_expression_subdocument

    def set_expression_scores_and_cell_names(self, expression_scores, cell_names):
        non_zero_cell_names = []
        non_zero_expression_scores = []
        for idx, value in enumerate(expression_scores):
            if(value > 0):
                non_zero_cell_names.append(cell_names[idx])
                non_zero_expression_scores.append(value)
        if len(non_zero_expression_scores) == 0:
            return None, None
        else:
            return non_zero_cell_names, non_zero_expression_scores

    def get_expression_scores(self):
        """Get expression scores

        Args:
            None

        Returns:
            Expression scores as list
        """
        return self.expression_scores

    def create_gene_subdocuments(self):
        if self.is_larger_than_doc_size():
            self.determine_amount_of_chunks()
            for gene_docs in self.chunk_gene_expression_documents():
                print(f'Gene is : {self.name}.')

    def chunk_gene_expression_documents(self):
        # sum starts at 59 because key values take up 59 bytes of disk space
        sum = 59
        start_index = 0
        float_storage = 8

        for index, cell_name in enumerate(self.cell_names):
            sum = sum + len[cell_name] + float_storage
            if (sum - 1) > DOCUMENT_LIMIT_BYTES:
                end_index = index - 1
                yield {'cell_names': self.cell_names[start_index:end_index],
                       'expression_scores':  self.expression_scores[start_index:end_index],
                       'source_file_name': source_file_name,
                       'source_file_type': source_file_type,
                       }
                sum = 59
                start_index = index

        self.number_of_chunks = int(
            self.size_of_gene_expression_document / DOCUMENT_LIMIT_BYTES)
        print(f'Number of chunks are {self.number_of_chunks}')
        return self.number_of_chunks
