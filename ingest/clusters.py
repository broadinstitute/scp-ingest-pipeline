import logging
from dataclasses import dataclass
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

from bson.objectid import ObjectId
from mypy_extensions import TypedDict

try:
    from ingest_files import DataArray
    from annotations import Annotations
    from monitor import setup_logger, log
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import DataArray
    from .annotations import Annotations
    from .monitor import setup_logger, log


@dataclass
class DomainRanges(TypedDict):
    x: List
    y: List
    z: List = None


class Clusters(Annotations):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    LINEAR_DATA_TYPE = 'ClusterGroup'
    COLLECTION_NAME = 'cluster_groups'
    errors_logger = setup_logger(
        __name__ + '_errors', 'errors.txt', level=logging.ERROR
    )
    # General logger for class
    info_logger = setup_logger(__name__, 'info.txt')

    my_debug_logger = log(errors_logger)

    @dataclass
    class Model(TypedDict):
        name: str
        # 3d or 2d cluster
        cluster_type: str
        # List of dictionaries that describe all extra "annotation" columns
        cell_annotations: List
        study_file_id: ObjectId
        study_id: ObjectId
        # Hash containing min/max arrays for each axis in the cluster plot
        domain_ranges: DomainRanges = None

    def __init__(
        self,
        file_path: str,
        study_id: ObjectId,
        study_file_id: ObjectId,
        name: str,
        *,
        domain_ranges: Dict = None,
        **kwargs,
    ):
        Annotations.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, study_id, study_file_id
        )
        # Lowercase coordinate headers, expected for df merge
        for i, header in enumerate(self.headers):
            if header in ['X', 'Y', 'Z']:
                self.headers[i] = self.headers[i].lower()
        self.preprocess()
        self.determine_coordinates_and_cell_names()
        self.source_file_type = "cluster"
        self.cluster_type = (
            '3d'
            if (
                "z" in self.coordinates_and_cell_headers
                or "Z" in self.coordinates_and_cell_headers
            )
            else '2d'
        )
        self.name = name
        # Check if domain_ranges is an empty dictionary
        self.domain_ranges = domain_ranges if not (not domain_ranges) else None
        self.extra_log_params = {'study_id': self.study_id, 'duration': None}

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
        """Cluster files must have coordinates 'x' and 'y' in header
        :return: boolean True if coordinates are in header, otherwise False
        """
        lower_cased_headers = [header.lower() for header in self.headers]
        for coordinate in ('x', 'y'):
            if coordinate not in lower_cased_headers:
                msg = (
                    "Header must have coordinate values 'x' and 'y' (case insensitive)"
                )
                self.store_validation_issue('error', 'format', msg)
                return False
        return True

    def transform(self):
        """ Builds cluster data model"""
        # Array of Hash objects that describe all extra "annotation" columns
        cell_annotations = []
        # Iterate through all extra "annotation" column headers
        # (I.E. column headers that are not coordinates or where TYPE !=NAME)
        for annot_headers in self.annot_column_headers:
            self.info_logger.info(
                f'Transforming header {annot_headers[0]}',
                extra={'study_id': self.study_id, 'duration': None},
            )
            annot_name = annot_headers[0]
            annot_type = annot_headers[1]
            # column_type = annot[1]
            cell_annotations.append(
                {
                    'name': annot_name,
                    'type': annot_type,
                    'values': self.file[annot_headers].unique().tolist()
                    if annot_type == 'group'
                    else [],
                }
            )
        self.info_logger.info(
            f'Created model for {self.study_id}', extra=self.extra_log_params
        )
        return self.Model(
            name=self.name,
            cluster_type=self.cluster_type,
            cell_annotations=cell_annotations,
            study_file_id=self.study_file_id,
            study_id=self.study_id,
            domain_ranges=DomainRanges(**self.domain_ranges)
            if self.domain_ranges is not None
            else None,
        )

    def get_data_array_annot(self, linear_data_id):
        for annot_header in self.file.columns:
            self.info_logger.info(
                f'Creating data array for header {annot_header[0]}',
                extra={'study_id': self.study_id, 'duration': None},
            )
            yield from Clusters.set_data_array(
                annot_header[0],
                self.name,
                self.file[annot_header].tolist(),
                self.study_file_id,
                self.study_id,
                linear_data_id,
            )
            self.info_logger.info(
                f'Attempting to load cluster header : {annot_header[0]}',
                extra={'study_id': self.study_id, 'duration': None},
            )

    def can_subsample(self):
        # TODO: Add more subsample validations
        if self.has_z:
            return len(self.header) > 4
        else:
            return len(self.header) > 3

    @staticmethod
    def set_data_array(
        name: str,
        cluster_name: str,
        values: List,
        study_file_id: str,
        study_id: str,
        linear_data_id: str,
        array_index: int = 0,
        subsample_annotation: str = None,
        subsample_threshold: int = None,
    ):
        # Add any other arguments below data_array_attr
        data_array_attr = locals()
        BASE_DICT = {'linear_data_type': Clusters.LINEAR_DATA_TYPE}

        def get_cluster_attr(annot_name):
            cluster_group_types = {
                'name': {'name': "text", 'array_type': "cells"},
                'coordinates': {
                    'name': annot_name.lower(),
                    'array_type': "coordinates",
                },
                'annot': {'name': annot_name, 'array_type': "annotations"},
            }

            if annot_name.lower() == "name":
                result = cluster_group_types.get('name')
                return result
            elif annot_name.lower() in ("x", "y", "z"):
                result = cluster_group_types.get('coordinates')
                return result
            else:
                result = cluster_group_types.get('annot')
                return result

        cluster_attr = get_cluster_attr(name)
        # Remove 'name' from function aruguments
        del data_array_attr['name']
        # Merge BASE_DICT, cluster_attr & data_array_attr and return DataArray model
        return DataArray(
            **data_array_attr, **cluster_attr, **BASE_DICT
        ).get_data_array()
