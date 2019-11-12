import copy
from typing import List, Tuple  # noqa: F401

import numpy as np
from ingest_files import IngestFiles


class SubSample(IngestFiles):
    ALLOWED_FILE_TYPES = ['text/csv', 'text/plain', 'text/tab-separated-values']
    MAX_THRESHOLD = 100_000
    SUBSAMPLE_THRESHOLDS = [MAX_THRESHOLD, 20_000, 10_000, 1_000]

    def __init__(self, cluster_file, cell_metadata_file=None):
        IngestFiles.__init__(
            self, cluster_file, self.ALLOWED_FILE_TYPES, open_as='dataframe'
        )
        self.determine_coordinates_and_cell_names()
        self.cell_metadata_df = (
            IngestFiles(
                cell_metadata_file,
                self.ALLOWED_FILE_TYPES,
                open_as='dataframe',
                header=[0, 1],
            )
            if cell_metadata_file is not None
            else None
        )
        self.preproccess()

    def prepare_cell_metadata(self):
        """ Does an inner join on cell and cluster file """
        if self.cell_metadata_df is not None:
            self.cell_metadata_df.preproccess()
            self.merge_df(
                self.file[self.coordinates_and_cell_names], self.cell_metadata_df.file
            )
            self.determine_coordinates_and_cell_names()

    def determine_coordinates_and_cell_names(self):
        """Finds column names for coordinates, annotations, and cell names"""
        self.coordinates_and_cell_names = [
            annot[0]
            for annot in self.file.columns
            if annot[0].lower() in ('z', 'y', 'x', 'name')
        ]
        # annotation column names
        self.columns = [
            annot
            for annot in self.file.columns
            if annot[0].lower() not in ('z', 'y', 'x', 'name')
        ]

    def bin(self, annotation: Tuple[str, str]):
        """Creates bins for a given group

        Args:
            annotation: Tuple[str, str]
                This is the annotation for a single column. For example annotation
                would look like ('annotation_name', 'numeric') or ('annotation_name', 'group')

        Returns:
            bin: Tuple[Dict[str: dataframe]], Tuple[str, str]]
                The first tuple contains all the bins for a given column/
                annotation. It would look like {'unique_value1': filtered dataframe where rows=unique_value1}
                for group values and there can be up to 20 bins for numeric columns.
                The second value in the tuple is structured exactly like the input value.
            """
        bin = {}
        if 'group' in annotation:
            # get unique values in column
            unique_values = self.file[annotation].unique()

            for col_val in unique_values:
                # get subset of data where row is equal to the unique value
                subset = self.file[self.file[annotation] == col_val]
                bin[col_val] = subset[self.coordinates_and_cell_names]
        else:
            columns = copy.copy(self.coordinates_and_cell_names)
            # coordinates, cell names and annotation name
            columns.append(annotation[0])
            # Subset of df where header is [cell_names, x, y, z, <annot_name>]
            subset = self.file[columns].copy()
            subset.sort_values(by=[annotation], inplace=True)
            # Generates 20 bins
            for index, df in enumerate(np.array_split(subset, 20)):
                bin[str(index)] = df[self.coordinates_and_cell_names]
        return bin, annotation

    def subsample(self):
        """Subsamples groups across a given file"""
        sample_sizes = [
            sample_size
            for sample_size in self.SUBSAMPLE_THRESHOLDS
            if sample_size < len(self.file.index)
        ]
        for bins in map(self.bin, self.columns):
            # (name of current column)
            annotation_name = bins[1]
            # Holds bins for annotation
            # Looks like {"Unique value #1" : dataframe, "Unique value #2": dataframe,...}
            anotation_dict = bins[0]
            for sample_size in sample_sizes:
                group_size = len(anotation_dict.keys())
                # Dict of values for the x, y, and z coordinates
                points = {k: [] for k in self.coordinates_and_cell_names}
                num_per_group = int(sample_size / group_size)
                cells_left = sample_size
                # bin = ("unique value in column" : dataframe)
                for bin in self.return_sorted_bin(anotation_dict, annotation_name):

                    amount_of_rows = len(bin[1].index)
                    # If the amount of sampled values is larger
                    # than the whole array, take the whole array
                    if num_per_group > amount_of_rows:
                        amount_picked_rows = amount_of_rows
                    else:
                        amount_picked_rows = num_per_group
                    shuffled_df = (
                        bin[1]
                        .reindex(np.random.permutation(bin[1].index))
                        .sample(n=amount_picked_rows)
                    )
                    for column in shuffled_df:
                        points[column[0]].extend(shuffled_df[column].values.tolist())
                    # Subtract number of cells 'subsampled' from the number of cells left
                    cells_left -= amount_picked_rows
                    # For last bin sample the number of cells left over
                    if bin[0] == list(anotation_dict.keys())[-2]:
                        num_per_group = cells_left
                    else:
                        group_size -= 1
                        num_per_group = int(cells_left / (group_size))
                # returns tuple = (subsampled values as dictionary, annotation name, sample size )
            yield (points, annotation_name, sample_size)

    def return_sorted_bin(self, bin, annot_name):
        """Sorts binned groups in order of size from smallest to largest for group annotations """

        if 'group' in annot_name:
            return sorted(bin.items(), key=lambda x: len(x[1]))
        else:
            return bin.items()
