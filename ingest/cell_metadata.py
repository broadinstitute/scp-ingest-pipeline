"""Module for ingesting cell metadata files

DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.

PREREQUISITES
Must have python 3.6 or higher.
"""
import collections
import ntpath
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

from bson.objectid import ObjectId
from mypy_extensions import TypedDict

try:
    # Used when importing internally and in tests
    from annotations import Annotations
    from ingest_files import DataArray
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .annotations import Annotations
    from .ingest_files import DataArray


class CellMetadata(Annotations):
    ALLOWED_FILE_TYPES = ['text/csv', 'text/plain', 'text/tab-separated-values']
    COLLECTION_NAME = 'cell_metadata'

    def __init__(
        self,
        file_path: str,
        study_id: ObjectId,
        study_file_id: ObjectId,
        *args,
        **kwargs,
    ):

        self.study_accession = kwargs.pop("study_accession")
        Annotations.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, study_id, study_file_id
        )
        self.cell_names = []
        # lambda below initializes new key with nested dictionary as value and avoids KeyError
        self.issues = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.ontology = defaultdict(lambda: defaultdict(list))
        self.ontology_label = dict()
        self.cells = []
        self.extra_log_params = {'study_id': self.study_id, 'duration': None}
        self.preprocess()

    # This model pertains to columns from cell metadata files
    @dataclass
    class Model(TypedDict):
        # value from column header
        name: str
        annotation_type: str
        study_file_id: ObjectId
        study_id: ObjectId
        # unique values from "group" type annotations
        values: List

    # Will evolve to do cross file validation
    def validate(self):
        """ Runs all validation checks
        """
        return all([self.is_valid_format()])

    def is_valid_format(self):
        """Validates format by calling all format validation methods"""
        return all(
            [self.validate_header_for_coordinate_values(), self.validate_format()]
        )

    def validate_header_for_coordinate_values(self):
        """Cell metadata files should not have coordinates in header
        :return: boolean True if coordinates are not in header, otherwise False
        """
        lower_cased_headers = [header.lower() for header in self.headers]
        valid = not any(
            [coordinate in ('x', 'y', 'z') for coordinate in lower_cased_headers]
        )
        if valid:
            return True
        else:
            msg = 'Header names can not be coordinate values x, y, or z (case insensitive)'
            self.store_validation_issue('error', 'format', msg)
            return False

    def transform(self):
        """ Builds cell metadata model"""
        AnnotationModel = collections.namedtuple(
            'AnnotationModel', ['annot_header', 'model']
        )
        for annot_header in self.file.columns[0:]:
            annot_name = annot_header[0]
            annot_type = annot_header[1]
            yield AnnotationModel(
                annot_header,
                self.Model(
                    {
                        'name': annot_name,
                        'annotation_type': annot_type,
                        # unique values from "group" type annotations else []
                        'values': list(self.file[annot_header].unique())
                        if annot_type == 'group'
                        else [],
                        'study_file_id': self.study_file_id,
                        'study_id': self.study_id,
                    }
                ),
            )

    def set_data_array(self, linear_data_id: str, annot_header: str):
        data_array_attrs = locals()
        del data_array_attrs['annot_header']
        del data_array_attrs['self']
        annot_name = annot_header[0]
        head, tail = ntpath.split(self.file_path)
        base_data_array_model = {
            'cluster_name': tail or ntpath.basename(head),
            'values': list(self.file[annot_header]),
            'study_file_id': self.study_file_id,
            'study_id': self.study_id,
        }
        # This is an array (or group of arrays) of every cell
        if annot_name.lower() == 'name':
            base_data_array_model.update(
                {
                    'name': 'All Cells',
                    'array_type': 'cells',
                    'linear_data_type': 'Study',
                }
            )
        # data from cell metadata file that correspond to a column of data
        else:
            base_data_array_model.update(
                {
                    'name': annot_name,
                    'array_type': 'annotations',
                    'linear_data_type': 'CellMetadatum',
                }
            )
        return DataArray(**base_data_array_model, **data_array_attrs).get_data_array()

    def yield_by_row(self) -> None:
        """ Yield row from cell metadata file"""
        for row in self.file.itertuples(index=False):
            dict_row = row._asdict()
            yield dict_row.values()
