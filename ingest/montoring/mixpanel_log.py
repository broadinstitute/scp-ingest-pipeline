from timeit import default_timer as timer
from typing import List, Dict  # noqa: F401
import functools
from contextlib import ContextDecorator

from .metrics_service import MetricsService


class CustomMetricTimedNode(ContextDecorator):
    def __init__(
        self,
        event_name,
        study_accession,
        functionName,
        file_type,
        file_size,
        props=None,
    ):
        self.event_name = event_name
        self.metrics_model = {
            "studyAccession": study_accession,
            "functionName": functionName,
            "fileType": file_type,
            "fileSize": file_size,
            "appId": "single-cell-portal",
            **props,
        }

    def __enter__(self):
        self.start_time = timer()
        return self

    def __exit__(self, exc, value, tb):
        if not self.start_time:
            return
        end_time = timer()
        duration = end_time - self.start_time
        self.metrics_model["perfTime"] = duration
        MetricsService.log(self.event_name, self.metrics_model)


def custom_metric(event_name, get_study_fun, get_study_file_fn, props: Dict):
    def _decorator(f):

        func_name = f.__name__

        @functools.wraps(f)
        def _wrapper(*args, **kwargs):
            study = get_study_fun()
            study_file = get_study_file_fn()
            with CustomMetricTimedNode(
                event_name,
                study.STUDY_ACCESSION,
                func_name,
                study_file.FILE_TYPE,
                study_file.FILE_SIZE,
                props=props,
            ):
                return f(*args, **kwargs)

        return _wrapper

    return _decorator
