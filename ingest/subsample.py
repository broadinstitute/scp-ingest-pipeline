import copy
from typing import List, Tuple  # noqa: F401

import numpy as np

try:
    from annotations import Annotations
    from clusters import Clusters
    from cell_metadata import CellMetadata
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .annotations import Annotations
    from .clusters import Clusters
    from .cell_metadata import CellMetadata


class SubSample(Annotations):
    ALLOWED_FILE_TYPES = ['text/csv', 'text/plain', 'text/tab-separated-values']
    MAX_THRESHOLD = 100_000
    SUBSAMPLE_THRESHOLDS = [MAX_THRESHOLD, 20_000, 10_000, 1_000]

    def __init__(self, cluster_file, cell_metadata_file=None):
        Annotations.__init__(self, cluster_file, self.ALLOWED_FILE_TYPES)
        self.preprocess()
        self.determine_coordinates_and_cell_names()
        if cell_metadata_file is not None:
            self.cell_metadata = Annotations(
                cell_metadata_file, CellMetadata.ALLOWED_FILE_TYPES
            )

    def prepare_cell_metadata(self):
        """ Does an inner join on cell and cluster file """
        if self.cell_metadata is not None:
            self.cell_metadata.preprocess()
            self.merge_df(
                self.file[self.coordinates_and_cell_headers], self.cell_metadata.file
            )
            self.determine_coordinates_and_cell_names()

    def bin(self, annotation: Tuple[str, str], scope: str):
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
        # sample the annotation along with coordinates and cell names
        columns_to_sample = copy.copy(self.coordinates_and_cell_headers)
        if scope == 'cluster':
            columns_to_sample.append(annotation[0])
        if 'group' in annotation:
            # get unique values in column
            unique_values = self.file[annotation].unique()

            for col_val in unique_values:
                # get subset of data where row is equal to the unique value
                subset = self.file[self.file[annotation] == col_val]
                bin[col_val] = subset[columns_to_sample]
        else:
            columns = copy.copy(self.coordinates_and_cell_headers)
            # coordinates, cell names and annotation name
            columns.append(annotation[0])
            # Subset of df where header is [cell_names, x, y, z, <annot_name>]
            subset = self.file[columns].copy()
            subset.sort_values(by=[annotation], inplace=True)
            # Generates 20 bins
            for index, df in enumerate(np.array_split(subset, 20)):
                bin[str(index)] = df[columns_to_sample]
        return bin, annotation

    def subsample(self, scope):
        """Subsamples groups across a given file"""
        sample_sizes = [
            sample_size
            for sample_size in self.SUBSAMPLE_THRESHOLDS
            if sample_size < len(self.file.index)
        ]
        for bins in [self.bin(col, scope) for col in self.annot_column_headers]:

            amount_of_bins = len(bins[0].keys())
            # (name of current column)
            annotation_name = bins[1]
            # Holds bins for annotation
            # Looks like {"Unique value #1" : dataframe, "Unique value #2": dataframe,...}
            annotation_dict = bins[0]
            for sample_size in sample_sizes:
                group_size = len(annotation_dict.keys())
                # Dict of values for the x, y, and z coordinates
                points = {k: [] for k in self.coordinates_and_cell_headers}
                if scope == 'cluster':
                    points[annotation_name[0]] = []
                num_per_group = int(sample_size / group_size)
                cells_left = sample_size
                # bin = ("unique value in column" : dataframe)
                for idx, bin in enumerate(
                    self.return_sorted_bin(annotation_dict, annotation_name)
                ):
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
                    # add the current observed annotation to the points dict the amount
                    # of times it has been sampled
                    # points[annotation_name] = [bin[0] for i in range(amount_picked_rows)]
                    # Subtract number of cells 'subsampled' from the number of cells left
                    cells_left -= amount_picked_rows
                    # For last bin sample the number of cells left over
                    # Subtract 2 because 0 based
                    if idx == (amount_of_bins - 2):
                        num_per_group = cells_left
                    else:
                        group_size -= 1
                        if group_size > 1:
                            num_per_group = int(cells_left / group_size)
                # returns tuple = (subsampled values as dictionary, annotation name, sample size )
                yield (points, annotation_name, sample_size)

    def return_sorted_bin(self, bin, annot_name):
        """Sorts binned groups in order of size from smallest to largest for group annotations """

        if 'group' in annot_name:
            return sorted(bin.items(), key=lambda x: len(x[1]))
        else:
            return bin.items()

    def set_data_array(self, args, kwargs):
        return Clusters.set_data_array(*args, **kwargs)
