from mixpanel import Mixpanel
from timeit import default_timer as timer
import os
import abc
import inspect
from dataclasses import dataclass
from typing import List, Dict, Generator  # noqa: F401
from mongo_connection import MongoConnection
from mypy_extensions import TypedDict
import functools
from contextlib import ContextDecorator


class MixpanelLog:
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        # Project token in the project settings dialog
        # of the Mixpanel web application
        # token = os.environ.get("MIX_PANEL_TOKEN")
        self.mp = Mixpanel("5393bdc03aa94a3480fd2f8fc583e6f4")

    def post_event(self, study_accession, event_name, props):
        self.mp.track("SCP123", event_name, props)


@dataclass
class CustomMetricTimedNode(ContextDecorator, MixpanelLog):
    MONGO_CONNECTION = MongoConnection()

    @dataclass
    class Model(TypedDict):
        def __init__(
            self,
            name: str,
            study_accession: str,
            function_name: str,
            file_path: int,
            file_type: str,
            file_size: int,
            perf_time: int = None,
        ):
            self.name = name
            self.study_accession = study_accession
            self.function_name = function_name
            self.file_path = file_path
            self.file_type = file_type
            self.file_size = file_size
            self.perf_time = perf_time

    def __init__(
        self,
        name,
        func_name,
        file_type,
        func_values=None,
        cls=None,
        props=None,
        *args,
        **kwargs,
    ):
        MixpanelLog.__init__(self)
        file_path = func_values[0]
        file_size = os.path.getsize(file_path)
        study_accession = CustomMetricTimedNode.get_study_accession(cls.study_id)
        if file_size == 0:
            raise ValueError(f"{file_size} is empty: " + str(file_size))
        self.model = CustomMetricTimedNode.Model(
            {
                "name": name,
                "study_accession": study_accession,
                "function_name": func_name,
                "file_path": file_path,
                "file_type": file_type,
                "file_size": file_size,
                "perf_time": None,
            }
        )
        self.model.update(props)

    def get_study_accession(study_id):
        return CustomMetricTimedNode.MONGO_CONNECTION["study_accessions"].find(
            {"study_id": study_id}, {"study_id": 1, "_id": 0}
        )

    def __enter__(self):
        self.start_time = timer()
        return self

    def __exit__(self, exc, value, tb):
        if not self.start_time:
            return
        end_time = timer()
        duration = end_time - self.start_time
        self.model["perf_time"] = duration
        self.post_event("", "is_sorted", self.model)


def custom_metric(name, file_type, props: Dict):
    def _decorator(f):
        func_name = f.__name__

        @functools.wraps(f)
        def _wrapper(*args, **kwargs):
            class_name = list(args).pop(0)
            if not inspect.isclass(class_name):
                raise ValueError("Decorator can only be used with class methods")
            func_values = args[1:]
            with CustomMetricTimedNode(
                name,
                func_name,
                file_type,
                func_values=func_values,
                cls=class_name,
                props=props,
            ):
                return f(*args, **kwargs)

        return _wrapper

    return _decorator
