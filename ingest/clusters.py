import copy
from typing import Dict

from ingest_files import IngestFiles


class Clusters(IngestFiles):

    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    COLLECTION_NAME = "clusters"
    SUBCOLLECTION_NAME = "data"

    def __init__(
        self,
        file_path: str,
        file_id: str,
        study_accession: str,
        *,
        name: str,
        domain_ranges: Dict = None,
    ):
        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
        self.header = self.get_next_line(increase_line_count=False)
        # Second line in cluster is metadata_type
        self.metadata_types = self.get_next_line(increase_line_count=False)
        self.source_file_type = "cluster"
        self.has_z = "z" or "Z" in self.header
        self.cell_annotations = self.create_cell_annotations_field()
        # TODO: Populate the cell_annotations array when pandas is implemented
        self.top_level_doc = {
            "name": name,
            "cluster_type": "3d" if self.has_z else "2d",
            "cell_annotations": {},
            "study_accession": study_accession,
            "domain_ranges": domain_ranges,
            "points": self.amount_of_lines,
            "file_id": file_id,
        }
        self.cluster_subdocs = self.return_cluster_subdocs(self.header)

    def update_points(self):
        self.top_level_doc["points"] = self.amount_of_lines

    def update_cell_annotations_field(self):
        if self.cell_annotations.values():
            self.top_level_doc["cell_annotations"] = list(
                self.cell_annotations.values()
            )

    def create_cell_annotations_field(self):
        cell_annotations = {}

        def return_cell_annotations_dict(name):
            return copy.copy({"name": name, "type": "group", "values": []})

        for idx, header in enumerate(self.header[1:]):
            if self.metadata_types[
                idx + 1
            ].lower() == 'group' and header.lower() not in ('x', 'y', 'z'):
                cell_annotations[header] = {}
                dictionary = return_cell_annotations_dict(header)
                cell_annotations[header] = dictionary
        return cell_annotations

    def transform(self, row):
        """ Add data from cluster files into annotation subdocs in cluster data model"""

        for idx, column in enumerate(row):
            annotation = self.header[idx]

            # first index is cell name don't need to check annot type
            if idx != 0:
                if self.metadata_types[idx].lower() == "numeric":
                    column = round(float(column), 3)
                elif self.metadata_types[idx].lower() == "group":
                    if column not in self.cell_annotations[annotation]["values"]:
                        # Update cell_annotations field
                        self.cell_annotations[annotation]["values"].append(column)
            # perform a shallow copy

            annotation_value = copy.copy(self.cluster_subdocs[annotation]["values"])
            annotation_value.append(column)
            self.cluster_subdocs[annotation]["values"] = annotation_value

    def return_cluster_subdocs(self, headers):
        """Creates cluster_subdocs"""
        cluster_subdocs = {}
        for annot_name in self.header:
            value = annot_name.lower()
            if value == "name":
                cluster_subdocs[annot_name] = self.create_cluster_subdoc(
                    "text", "cells"
                )
            elif value in ("x", "y", "z"):
                cluster_subdocs[annot_name] = self.create_cluster_subdoc(
                    value, "coordinates"
                )
            else:
                cluster_subdocs[annot_name] = self.create_cluster_subdoc(
                    annot_name, "annotations"
                )
        return cluster_subdocs

    @staticmethod
    def create_cluster_subdoc(
        annot_name,
        header_value_type,
        *,
        values=[],
        subsample_annotation=None,
        subsample_threshold=None,
    ):
        """Returns cluster subdoc"""

        return {
            "name": annot_name,
            "array_index": 0,
            "values": values,
            "array_type": header_value_type,
            "subsample_annotation": subsample_annotation,
            "subsample_threshold": subsample_threshold,
        }

    def can_subsample(self):
        # TODO: Add more subsample validations
        if self.has_z:
            return len(self.header) > 4
        else:
            return len(self.header) > 3
