# File is responsible for defining globals and initializing them
from mongo_connection import MongoConnection

MONGO_CONNECTION = MongoConnection()


class Study:
    STUDY_ACCESSION: str = None

    def __init__(self, study_id):
        Study.update_study_accession(study_id)

    @classmethod
    def update_file_type(cls, file_type):
        cls.FILE_TYPE = file_type

    @classmethod
    def update_study_accession(cls, study_id):
        study_accession = MONGO_CONNECTION["study_accessions"].find(
            {"study_id": study_id}, {"study_id": 1, "_id": 0}
        )
        cls.STUDY_ACCESSION = study_accession


class StudyFile:
    FILE_TYPE: str = None
    FILE_SIZE: int = None
    STUDY_FILE_ID = None

    def __init__(self, study_file_id):
        StudyFile.update(study_file_id)

    @classmethod
    def update(cls, study_file_id):
        cls.STUDY_FILE_ID = study_file_id
        MONGO_CONNECTION["study_accessions"].find(
            {"_id": study_file_id}, {"file_type": 1, "upload_file_size ": 1, "_id": 0}
        )


def init(study_id, study_file_id):
    global study
    study = Study(study_id)
    global study_file
    study_file = StudyFile(study_file_id)
