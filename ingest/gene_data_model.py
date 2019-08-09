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

# Actual document limit is 1 MiB (1,048,576 bytes) per
# https://cloud.google.com/firestore/quotas#writes_and_transactions
# However, due to an unclear calculation issue, we use a limit that is roughly
# 90% smaller to avoid exceeding this Firestore constraint.
#
# TODO: Reconcile calculations and use documented size limit (1_048_576)
DOCUMENT_LIMIT_BYTES = 1_048_576


class Gene:
    COLLECTION_NAME = 'Gene'
    SUBCOLLECTION_NAME = 'gene_expression'

    def __init__(self, name: str, source_file_type: str, *,
                 gene_id: str = '', study_accession: str = '', taxon_name: str = '',
                 taxon_common_name: str = '', ncbi_taxid: str = '', genome_assembly_accession: str = '',
                 genome_annotation: str = '', cell_names: List[str] = [], expression_scores: List = [],
                 check_for_zero_values: bool = True, field_id: str = '', id: str = '') -> None:

        self.name = name.replace('"', '')
        self.gene_id = gene_id
        self.source_file_type = source_file_type,

        if check_for_zero_values:
            self.cell_names, self.expression_scores = self.set_expression_scores_and_cell_names(
                expression_scores, cell_names)
        else:
            self.cell_names = cell_names
            self.expression_scores = expression_scores
        # Subdocument that contains all cell names and expression scores for a
        # given gene
        self.subdocument = {'cell_names': self.cell_names,
                            'expression_scores':  self.expression_scores,
                            'source_file_type': source_file_type,
                            }
       # This is the top level document for the gene data model
        self.top_level_doc = {
            'field_id': field_id,
            'id': id,
            'name': self.name,
            'gene_id': self.gene_id,
            'study_accession': study_accession,
            'taxon_name': taxon_name,
            'taxon_common_name': taxon_common_name,
            'ncbi_taxid': ncbi_taxid,
            'genome_assembly_accession': genome_assembly_accession,
            'genome_annotation': genome_annotation,
        }
        print(cell_names)

    def has_subcollection_data(self):
        return self.cell_names != None

    def set_expression_scores_and_cell_names(self, expression_scores, cell_names):
        """Sets expression scores and cell names. Only keeps cell names that have
            non-zero expression values. Turns expression values into floats.

        Args:
            cell_names: List[str]
                All cell names for a given gene.
            expression_scores: List
                All expression scores for a given gene.

        Returns:
            non_zero_cell_names: List[str]
                Cell names that have non-zero expression values
            non_zero_expression_scores: List[float]
                Expression scores that are greater than 0.
        """
        non_zero_cell_names = []
        non_zero_expression_scores = []
        for idx, value in enumerate(expression_scores):
            expression_score = round(float(value), 3)
            if expression_score > 0:

                non_zero_cell_names.append(cell_names[idx])
                non_zero_expression_scores.append(expression_score)
        if len(non_zero_expression_scores) == 0:
            return None, None
        else:
            return non_zero_cell_names, non_zero_expression_scores

    def chunk_gene_expression_documents(self, document_name):
        """Partitions gene expression documents in storage sizes that are
            less than 1,048,576 bytes. Storage size calculation figures are derived from:
            # https://cloud.google.com/firestore/docs/storage-size

        Args:
            None

        Returns:
            Dictionary that consist of expression scores and cell names,
            file name and type.
        """

        # starting sum is equal to 34 (because key names take up 59 bytes) plus the
        # storage size of the source file name and file type
        size_of_document_name = len(document_name) + 1
        size_of_cell_names_field = 10 + 1
        size_of_expression_scores_field = 17 + 1
        starting_sum = 17 + \
            len(self.source_file_type) + size_of_document_name
        start_index = 0
        float_storage = 8
        sum = starting_sum

        for index, cell_name in enumerate(self.cell_names):

            cell_name_storage = len(cell_name) + 1 + size_of_cell_names_field
            expression_scores_storage = size_of_expression_scores_field + float_storage
            sum = sum + expression_scores_storage + cell_name_storage
            # Subtract 32 based off of firestore storage guidelines for strings
            # and documents
            # This and other storage size calculation figures are derived from:
            # https://cloud.google.com/firestore/docs/storage-size
            if (sum + 32) > DOCUMENT_LIMIT_BYTES or index == (len(self.cell_names) - 1):
                if index == (len(self.cell_names) - 1):
                    end_index = index
                else:
                    end_index = index - 1

                yield {'cell_names': self.cell_names[start_index:end_index],
                       'expression_scores':  self.expression_scores[start_index:end_index],
                       'source_file_type': self.source_file_type[0],
                       }
                print(sum)
               # Reset sum and add storage size at current index
                sum = starting_sum + cell_name_storage + expression_scores_storage
                start_index = index
