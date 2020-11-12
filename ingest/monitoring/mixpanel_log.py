from timeit import default_timer as timer
from typing import List, Dict  # noqa: F401
import functools
from contextlib import ContextDecorator

try:
    from monitor import testing_guard
except ImportError:
    from ..monitor import testing_guard


class MetricTimedNode(ContextDecorator):
    """Context manager that adds properties and performance times of Mixpanel events."""

    def __init__(
        self, function_name, perf_time_name, metric_properties: Dict, props: Dict = None
    ):
        props["functionName"] = function_name
        self.props = props
        self.perf_time_name = perf_time_name
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


# @testing_guard is applied here so Mixpanel events are not logged during tests.
@testing_guard
def custom_metric(
    get_metric_properties_fn, perf_time_name="perfTime", props: Dict = {}
):
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
            ):
                return f(*args, **kwargs)

        return _wrapper

    return _decorator
