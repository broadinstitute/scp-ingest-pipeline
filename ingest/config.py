import os

# File is responsible for defining globals and initializing them
try:
    from mongo_connection import MongoConnection
except ImportError:
    from .mongo_connection import MongoConnection
from bson.objectid import ObjectId

MONGO_CONNECTION = MongoConnection()


def init(study_id, study_file_id, user_metric_uuid=None):
    global __metric_properties

    study = Study(study_id)
    study_file = StudyFile(study_file_id)
    __metric_properties = MetricProperties(study, study_file, user_metric_uuid)


def set_parent_event_name(event_name):
    global __event_name
    __event_name = event_name


def get_parent_event_name():
    global __event_name
    return __event_name


def get_metric_properties():
    global __metric_properties
    return __metric_properties


class MetricProperties:
    """Sets default properties for Mixpanel events"""

    # This is a generic write-only log token, not a secret
    USER_ID = "2f30ec50-a04d-4d43-8fd1-b136a2045079"

    def __init__(self, study, study_file, user_uuid=None):
        distinct_id = user_uuid if user_uuid else MetricProperties.USER_ID
        self.__properties = {
            "distinct_id": distinct_id,
            "studyAccession": study.accession,
            "fileName": study_file.file_name,
            "fileType": study_file.file_type,
            "fileSize": study_file.file_size,
            "trigger": study_file.trigger,
            "logger": "ingest-pipeline",
            "appId": "single-cell-portal",
        }

    def get_properties(self):
        return self.__properties

    def update(self, props):
        if props:
            self.__properties = {**self.__properties, **props}
        self.set_mixpanel_nums()

    def set_mixpanel_nums(self):
        """Derive count for each type of Mixpanel property"""
        for prop in ["errorTypes", "errors", "warningTypes", "warnings"]:
            num_prop = "num" + prop[0].upper() + prop[1:]
            if self.__properties.get(prop):
                self.__properties[num_prop] = len(self.__properties[prop])

    def append_issue(self, props):
        """Add error/warning properties to MetricsProperties
        without clobbering
        """

        if props:
            updated_props = {}
            for prop in ["errorTypes", "errors", "warningTypes", "warnings"]:
                if prop in self.__properties:
                    updated_props.setdefault(prop, []).extend(self.__properties[prop])
                if prop in props:
                    updated_props.setdefault(prop, []).extend(props[prop])
            for key in updated_props.keys():
                if key in ["errorTypes", "warningTypes"]:
                    self.__properties[key] = list(set(updated_props[key]))
                else:
                    self.__properties[key] = updated_props[key]
            self.set_mixpanel_nums()


def bypass_mongo_writes():
    """Check if developer has set environment variable to bypass writing data to MongoDB
    BYPASS_MONGO_WRITES='yes'
    """
    if os.environ.get("BYPASS_MONGO_WRITES") is not None:
        skip = os.environ["BYPASS_MONGO_WRITES"]
        if skip == "yes":
            return True
        else:
            return False
    else:
        return False


class Study:
    """Provides attributes for a given study"""

    def __init__(self, study_id):
        self.study = study_id

    @property
    def study(self):
        return self.__study

    @study.setter
    def study(self, study_id: str):
        # Annotation Object expects a proper BSON ID
        # even when the input validation is bypassed here
        try:
            study_id = ObjectId(study_id)
        except Exception:
            raise ValueError("Must pass in valid object ID for study ID.")
        # set dummy accession if running in developer mode
        if bypass_mongo_writes():
            self.accession = "SCPdev"
        else:
            study = list(
                MONGO_CONNECTION._client["study_accessions"].find(
                    {"study_id": study_id}, {"_id": 0}
                )
            )
            if not study:
                raise ValueError(
                    "Study ID is not registered with a study. Please provide a valid study ID."
                )
            else:
                self.__study = study.pop()
                self.accession = self.__study["accession"]


class StudyFile:
    "Provides attributes for a given study file"

    def __init__(self, study_file_id):
        self.study_file = study_file_id

    @property
    def study_file(self):
        return self.__study_file

    @study_file.setter
    def study_file(self, study_file_id):
        try:
            study_file_id = ObjectId(study_file_id)
        except Exception:
            raise ValueError("Must pass in valid object ID for study file ID.")
        if bypass_mongo_writes():
            # set dummy values if running in developer mode
            self.file_type = "input_validation_bypassed"
            self.file_size = 1
            self.file_name = str(study_file_id)
            self.trigger = 'dev-mode'
        else:
            query = MONGO_CONNECTION._client["study_files"].find({"_id": study_file_id})
            query_results = list(query)
            if not query_results:
                raise ValueError(
                    "Study file ID is not registered with a study. Please provide a valid study file ID."
                )
            else:
                self.__study_file = query_results.pop()
                self.file_type = self.study_file["file_type"]
                self.file_size = self.study_file["upload_file_size"]
                self.file_name = self.study_file["name"]
                # when set, remote_location is the name of the file in the bucket
                if self.study_file.get("remote_location") is not None:
                    if self.study_file["remote_location"] == "":
                        self.trigger = 'upload'
                    else:
                        self.trigger = 'sync'
                # indicate trigger state for tests/mocks
                else:
                    self.trigger = 'not set'
