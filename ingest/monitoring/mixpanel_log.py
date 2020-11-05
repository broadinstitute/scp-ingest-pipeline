from timeit import default_timer as timer
from typing import List, Dict  # noqa: F401
import functools
from contextlib import ContextDecorator

from settings import update_model, add_child_event


class MetricTimedNode(ContextDecorator):
    """Context manager that logs properties and performance times of Mixpanel events."""

    def __init__(
        self, function_name, perf_time_name, props=None, child_event_name=None
    ):
        props["functionName"] = function_name
        self.child_event_name = child_event_name
        self.props = props
        self.perf_time_name = perf_time_name if perf_time_name else "perfTime"

    def __enter__(self):
        self.start_time = timer()
        return self

    def __exit__(self, exc, value, tb):
        end_time = timer()
        # Calculate performance time
        duration = end_time - self.start_time
        self.props[self.perf_time_name] = round(duration, 3)
        # Add current functions properties to base model
        update_model(self.props)
        # Add and child events to list of child events captured
        if self.child_event_name:
            add_child_event(self.child_event_name)


def custom_metric(child_event_name=None, perf_time_name=None, props: Dict = {}):
    def _decorator(f):
        func_name = f.__name__

        @functools.wraps(f)
        def _wrapper(*args, **kwargs):
            with MetricTimedNode(
                func_name,
                perf_time_name,
                props=props,
                child_event_name=child_event_name,
            ):
                return f(*args, **kwargs)

        return _wrapper

    return _decorator
