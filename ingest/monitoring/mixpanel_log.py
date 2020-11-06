from timeit import default_timer as timer
from typing import List, Dict  # noqa: F401
import functools
from contextlib import ContextDecorator

from settings import add_child_event


class MetricTimedNode(ContextDecorator):
    """Context manager that logs properties and performance times of Mixpanel events."""

    def __init__(
        self,
        function_name,
        perf_time_name,
        metric_properties: Dict,
        props: Dict = None,
        child_event_name: str = None,
    ):
        props["functionName"] = function_name
        self.child_event_name = child_event_name
        self.props = props
        self.perf_time_name = perf_time_name if perf_time_name else "perfTime"
        self.metric_properties = metric_properties

    def __enter__(self):
        self.start_time = timer()
        return self

    def __exit__(self, exc, value, tb):
        end_time = timer()
        # Calculate performance time
        duration = end_time - self.start_time
        self.props[self.perf_time_name] = round(duration, 3)
        # Add current functions properties to base model
        self.metric_properties.update(self.props)
        # Add and child events to list of child events captured
        if self.child_event_name:
            add_child_event(self.child_event_name)


def custom_metric(
    get_metric_properties_fn,
    child_event_name=None,
    perf_time_name=None,
    props: Dict = {},
):
    if child_event_name or perf_time_name or props:
        # Confirms child events have a perfTime name and associated properties
        if not all([child_event_name, perf_time_name, props]):
            raise ValueError(
                "A child event name and or properties must have a perfTime name"
            )

    def _decorator(f):
        func_name = f.__name__

        @functools.wraps(f)
        def _wrapper(*args, **kwargs):
            metric_properties = get_metric_properties_fn()
            with MetricTimedNode(
                func_name,
                perf_time_name,
                metric_properties=metric_properties,
                props=props,
                child_event_name=child_event_name,
            ):
                return f(*args, **kwargs)

        return _wrapper

    return _decorator
