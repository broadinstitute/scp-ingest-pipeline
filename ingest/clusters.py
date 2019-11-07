from typing import Dict, Generator, List, Tuple, Union  # noqa: F401
from dataclasses import dataclass
from mypy_extensions import TypedDict

from ingest_files import DataArray
from annotations import Annotations


@dataclass
class DomainRanges(TypedDict):
    x: List
    y: List
    z: List = None


class Clusters(Annotations):
    # Need to delete after refactor
    COLLECTION_NAME = "cluster_groups"
    SUBCOLLECTION_NAME = "data"

    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    LINEAR_DATA_TYPE = 'ClusterGroup'

    @dataclass
    class Document(TypedDict):
        name: str
        # 3d or 2d cluster
        cluster_type: str
        # List of dictionaries that describe all extra "annotation" columns
        cell_annotations: List
        file_id: str
        study_id: str
        # Hash containing min/max arrays for each axis in the cluster plot
        domain_ranges: DomainRanges = None

    def __init__(
        self,
        file_path: str,
        study_id: str,
        study_file_id: str,
        name: str,
        *,
        domain_ranges: Dict = None,
    ):
        Annotations.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
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
        # self.cell_annotations = self.create_cell_annotations_field()
        self.cluster_subdocs = {}
        self.name = name
        self.domain_ranges = domain_ranges
        self.study_id = study_id
        self.study_file_id = study_file_id
        # TODO: Populate the cell_annotations array when pandas is implemented
        self.top_level_doc = {
            "name": name,
            "cell_annotations": [],
            "study_accession": study_id,
            "domain_ranges": domain_ranges,
            "points": self.amount_of_lines,
            "study_file_id": study_file_id,
        }
        self.model = self.create_cluster_model()
        # for annot_name in self.header:
        #     model = Clusters.return_cluster_subdocs(annot_name)
        #     self.cluster_subdocs.update(model)

    def update_points(self):
        self.top_level_doc["points"] = self.amount_of_lines

    def update_cell_annotations_field(self):
        if self.cell_annotations.values():
            self.top_level_doc["cell_annotations"] = list(
                self.cell_annotations.values()
            )

    def transform(self):
        """ Add data from cluster files into annotation subdocs in cluster data model"""
        for annot_header in self.file.columns:
            annot_name = annot_header[0]
            yield Clusters.set_data_array(
                annot_name,
                self.name,
                list(self.file[annot_header]),
                self.study_file_id,
                self.study_id,
                None,
            )

    # set_data_array() replaces this method
    @staticmethod
    def return_cluster_subdocs(
        annot_name, *, values=[], subsample_annotation=None, subsample_threshold=None
    ):
        """Creates cluster_subdocs"""
        cluster_subdocs = {}
        value = annot_name.lower()
        if value == "name":
            cluster_subdocs[annot_name] = Clusters.create_cluster_subdoc(
                "text",
                "cells",
                values=values,
                subsample_annotation=subsample_annotation,
                subsample_threshold=subsample_threshold,
            )
        elif value in ("x", "y", "z"):
            cluster_subdocs[annot_name] = Clusters.create_cluster_subdoc(
                value,
                "coordinates",
                values=values,
                subsample_annotation=subsample_annotation,
                subsample_threshold=subsample_threshold,
            )
        else:
            cluster_subdocs[annot_name] = Clusters.create_cluster_subdoc(
                annot_name,
                "annotations",
                values=values,
                subsample_annotation=subsample_annotation,
                subsample_threshold=subsample_threshold,
            )

        if subsample_threshold is not None:
            return cluster_subdocs[annot_name]

        else:
            return cluster_subdocs

    # set_data_array() replaces this method
    @staticmethod
    def create_cluster_subdoc(
        annot_name,
        annot_type,
        values=[],
        subsample_annotation=None,
        subsample_threshold=None,
    ):
        """Returns cluster subdoc"""

        return {
            "name": annot_name,
            "array_index": 0,
            "values": values,
            "array_type": annot_type,
            "subsample_annotation": subsample_annotation,
            "subsample_threshold": subsample_threshold,
        }

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
        # Remove 'name' from data pass into function
        del data_array_attr['name']
        # print(data_array_attr)

        # merged_array = {**data_array_attr, **cluster_attr, **BASE_DICT}
        # Merge BASE_DICT, cluster_attr & data_array_attr and return DataArray model
        return DataArray({**data_array_attr, **cluster_attr, **BASE_DICT})

    def create_cluster_model(self):
        # Array of Hash objects that describe all extra "annotation" columns
        cell_annotations = []
        # Grab all columns from df where annot_type = group
        for annot in self.file.xs('group', axis=1, level=1, drop_level=False).columns:
            col_name = annot[0]
            # column_type = annot[1]
            cell_annotations.append(
                {
                    'name': col_name,
                    'type': 'group',
                    'values': list(self.file[annot].unique()),
                }
            )
        return self.Document(
            name=self.name,
            cluster_type=self.cluster_type,
            cell_annotations=cell_annotations,
            study_file_id=self.study_file_id,
            study_id=self.study_id,
            domain_ranges=DomainRanges(**self.domain_ranges),
        )
