import copy
import json
import pprint
import random

import numpy as np
import pandas
from clusters import Clusters
from ingest_files import IngestFiles


class SubSample(IngestFiles):
    ALLOWED_FILE_TYPES = ['text/csv',
                          'text/plain', 'text/tab-separated-values']
    MAX_THRESHOLD = 100_000
    SUBSAMPLE_THRESHOLDS = [MAX_THRESHOLD, 20000, 600, 200]

    def __init__(self, file_path, annot_scope):
        IngestFiles.__init__(self,
                             file_path, self.ALLOWED_FILE_TYPES, open_as='pandas')
        self.header = list(self.file.head())
        self.annot_scope = annot_scope
        self.has_z = 'z' in (annot[0].lower()
                             for annot in self.header[0])
        self.coordinates_and_cell_names,  self.columns = self.dermine_coordinates_and_cell_names()
        columns = self.file.xs("numeric", axis=1, level=1,
                               drop_level=False).columns.tolist()
        self.file[columns] = self.file[columns].round(2)
        self.round_numeric_annot()
        # print(self.file)
        # conver_group_to_string()

    def round_numeric_annot(self):
        numeric_columns = self.file.xs("numeric", axis=1, level=1,
                                       drop_level=False).columns.tolist()
        self.file[numeric_columns] = self.file[numeric_columns].round(2)

    def dermine_coordinates_and_cell_names(self):
        if self.annot_scope is 'cluster':
            if self.has_z:
                return [header_name[0] for header_name in self.header[:3]], self.file.iloc[:, 4:]
            else:
                return [header_name[0] for header_name in self.header[:3]], self.file.iloc[:, 3:]
        else:
            return header_name[0]

    def bin(self, annotation):
        bin = {}
        if 'group' in annotation:
            # get unique values in column
            unique_values = self.file[annotation].unique()

            for col_val in unique_values:
                # get subset of data where
                subset = self.file[self.file[annotation] == col_val]
                bin[col_val] = subset[self.coordinates_and_cell_names]
        else:
            columns = copy.copy(self.coordinates_and_cell_names)
            columns.append(annotation[0])
            subset = self.file[columns]
            subset.sort_values(
                by=[annotation], inplace=True)
            for index, df in enumerate(np.array_split(subset, 20)):
                bin[str(index)] = df[self.coordinates_and_cell_names]
        return bin, annotation

    def subsample(self):

        sample_sizes = [
            sample_size for sample_size in self.SUBSAMPLE_THRESHOLDS if sample_size < len(self.file.index)]
        for bins in map(self.bin, self.columns):
            # (name of current collumn)
            annotation_name = bins[1]
            # Looks like {"Unique value #1" : dataframe, "Unique value #2": dataframe,...}
            anotation_dict = bins[0]
            for sample_size in sample_sizes:
                group_size = len(anotation_dict.keys())
                # values for the x, y, and z coordinates
                points = {k: [] for k in self.coordinates_and_cell_names}
                num_per_group = int(sample_size / group_size)

                # bin = ("unique value in column" : dataframe)
                for bin in self.return_sorted_bin(anotation_dict, annotation_name):
                    cells_left = sample_size
                    amount_of_rows = bin[1].shape[0]
                    # If the amounf of sampled values we need is larger
                    # than the whole array, take the whole array
                    if num_per_group > amount_of_rows:
                        amount_picked_rows = amount_of_rows
                    else:
                        amount_picked_rows = num_per_group
                    suffled_df = bin[1].reindex(
                        np.random.permutation(bin[1].index)).sample(n=amount_picked_rows)
                    for column in suffled_df:
                        points[column[0]].extend(
                            suffled_df[column].values.tolist())
                    # Subtract number of cells 'subsampled' from the number of cells left
                    cells_left -= amount_picked_rows
                    if bin[0] == list(anotation_dict.keys())[-2]:
                        num_per_group = cells_left
                    else:
                        group_size -= 1
                        num_per_group = int(cells_left / (group_size))
                # returns tuple = (subsampled values: dict, annotation name, sample size )
                yield (points, annotation_name, sample_size)

    def return_sorted_bin(self, bin, annot_name):

        if('group' in annot_name):
            return sorted(bin.items(), key=lambda x: len(x[1]))
        else:
            return bin.items()
