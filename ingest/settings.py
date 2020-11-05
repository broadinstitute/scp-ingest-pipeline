# File is responsible for defining globals and initializing them
from mongo_connection import MongoConnection
from bson.objectid import ObjectId
from typing import Dict

MONGO_CONNECTION = MongoConnection()


def init(study_id, study_file_id):
    global child_events
    child_events = []
    global study
    study = Study(study_id)
    global study_file
    study_file = StudyFile(study_file_id)
    global metrics_model
    metrics_model = {
        "studyAccession": study.STUDY_ACCESSION,
        "fileName": study_file.FILE_NAME,
        "fileType": study_file.FILE_TYPE,
        "fileSize": study_file.FILE_SIZE,
        "appId": "single-cell-portal",
    }


def add_child_event(event):
    child_events.append(event)


def get_study():
    global study
    return study


def get_study_file():
    global study_file
    return study_file


def update_model(props: Dict):
    global metrics_model
    metrics_model = {**metrics_model, **props}


def get_metrics_model():
    global metrics_model
    return metrics_model


class Study:
    STUDY_ACCESSION: str = None

    def __init__(self, study_id):
        Study.update_study_accession(study_id)

    @classmethod
    def update_file_type(cls, file_type):
        cls.FILE_TYPE = file_type

    @classmethod
    def update_study_accession(cls, study_id):
        study_accession = list(
            MONGO_CONNECTION._client["study_accessions"].find(
                {"study_id": ObjectId(study_id)}, {"accession": 1, "_id": 0}
            )
        ).pop()
        cls.STUDY_ACCESSION = study_accession["accession"]


class StudyFile:
    FILE_TYPE: str = None
    FILE_SIZE: int = None
    STUDY_FILE_ID = None
    FILE_NAME = None

    def __init__(self, study_file_id):
        StudyFile.update(study_file_id)

    @classmethod
    def update(cls, study_file_id):
        cls.STUDY_FILE_ID = study_file_id
        query = MONGO_CONNECTION._client["study_files"].find(
            {"_id": ObjectId("5e3dd9bcfa6429263aab24ba")}
        )
        query_results = list(query).pop()
        cls.FILE_TYPE = query_results["file_type"]
        cls.FILE_SIZE = query_results["upload_file_size"]
        cls.FILE_NAME = query_results["name"]
