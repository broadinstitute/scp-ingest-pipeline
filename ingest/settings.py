# File is responsible for defining globals and initializing them
from .mongo_connection import MongoConnection


class Study:
    MONGO_CONNECTION = MongoConnection()
    STUDY_ACCESSION = None

    def __init__(self, study_id):
        Study.upadate_study_accession(study_id)

    @classmethod
    def update_study_accession(cls, study_id):
        study_accession = Study.MONGO_CONNECTION["study_accessions"].find(
            {"study_id": study_id}, {"study_id": 1, "_id": 0}
        )
        cls.STUDY_ACCESSION = study_accession


def init(study_id):
    global study
    study = Study(study_id)
