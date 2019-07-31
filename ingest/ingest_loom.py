"""Command-line interface for ingesting Loom files into Firestore

DESCRIPTION
This CLI passes extract and transforms functions into Ingest Pipeline thereby
allowing the Loom file to be ingested into Firestore.

PREREQUISITES
Must have python 3.6 or higher.

EXAMPLES
# Takes Loom file and stores it into Firestore
$ python ingest_loom.py <path to Loom file>
$ python ingest_loom.py ../tests/data/L5_All.agg.loom
"""

import argparse
import os
import sys
from typing import *

import loompy
import numpy as np
from gene_data_model import Gene


class Loom:
    def __init__(self, file_path):

        self.loom = loompy.connect(file_path)
        self.file_name, self.filetype = os.path.splitext(file_path)

    def extract(self, size=500) -> Iterable[Tuple[int, np.ndarray, loompy.LoomView]]:
        """Reads Loom file and extracts out batches of rows

        Args:
            size: int
                                 the chunks returned at every element of the iterator. Default
                     is 500 rows at a time

        Returns:
        ------
            Taken from https://github.com/linnarsson-lab/loompy/blob/master/loompy/loompy.py
            Iterable that yields triplets of (ix, indexes, view) where
                ix: int
                        first position / how many rows/cols have been yielded alredy
                indexes: np.ndarray[int]
                        the indexes with the same numbering of the input args cells / genes (i.e. ``np.arange(len(ds.shape[axis]))``)
                        this is ``ix + selection``
                view: LoomView
                        a view corresponding to the current chunk
        """
        for (ix, selection, view) in self.loom.scan(axis=0, batch_size=size):
            yield ix, selection, view

    def transform(self, ix, selection, view) -> List[Gene]:
        """Transforms Loom file into Firestore data model

        Args:
            Taken from https://github.com/linnarsson-lab/loompy/blob/master/loompy/loompy.py
            ix: int
                        first position / how many rows/cols have been yielded alredy
                indexes: np.ndarray[int]
                        the indexes with the same numbering of the input args cells / genes
                view: LoomView
                        a view corresponding to the current chunk

        Returns:
        ------
                    transformed_data : List[Gene]
                A list of Gene objects
        """
        transformed_data = []
        for index in selection:
            position = index - ix
            expression_scores = [float(i) for i in view[position, 0:20]]
            gene_model = Gene(view.ra.Gene[index - ix], self.file_name, self.filetype,
                              gene_id=view.ra.Accession[position],
                              expression_scores=expression_scores)
            transformed_data.append(gene_model.gene)
        return transformed_data

    # def ingest(self) -> None:
    #     """ Ingests Loom file via IngestPipeline
    #
    #     Args:
    #         Nothing
    #     Returns:
    #     ------
    #         Nothing
    #     """
    #     ingest_pipeline = ingest.connect(self.extract, self.transform)
    #     ingest_pipeline.ingest()


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         prog ='ingest_loom.py',
#         description= __doc__,
#         formatter_class=argparse.RawDescriptionHelpFormatter
#     )
#
#     #Positional argument
#     parser.add_argument(
#         "loom_file",
#         help='Absolute or relative path to loom expression file'
#     )
#     args = parser.parse_args()
#
#     loom_object = Loom(args.loom_file)
#     loom_object.ingest()
#     loom_object.loom.close()
