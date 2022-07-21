import logging
import os
from contextlib import nullcontext

# General docs for Sentry:
# https://docs.sentry.io/platforms/python/logging/

# API docs:
# https://getsentry.github.io/sentry-python/api.html
# https://getsentry.github.io/sentry-python/integrations.html

# See logs in Sentry:
# https://sentry.io/organizations/broad-institute/issues/?project=1424198
import sentry_sdk
from sentry_sdk.integrations.logging import LoggingIntegration


def support_configs(log_file):
    """Returns handler for loggers used for SCP support"""
    formatter = logging.Formatter(
        "%(asctime)s %(name)s %(levelname)s:%(message)s", datefmt="%Y-%m-%dT%H:%M:%S%z"
    )
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    return handler


def default_configs(log_file):
    handler = logging.FileHandler(log_file)
    return handler


def setup_logger(logger_name, log_file, level=logging.DEBUG, format="default"):
    """"Sets up logging configurations and formatting"""
    logger = logging.getLogger(logger_name)
    if format == "support_configs":
        handler = support_configs(log_file)

    else:
        handler = default_configs(log_file)
        logger.propagate = True
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger


def log_exception(dev_logger, user_logger, exception):
    user_logger.critical(str(exception))
    dev_logger.exception(exception)


# Modified from https://jdkandersson.com/2019/05/19/testing-decorated-python-functions/
def testing_guard(decorator_func):
    """
    Decorator that only applies another decorator if the TESTING environment
    variable is not set.

    Args:
        decorator_func: The decorator function.

    Returns:
        Function that calls a function after applying the decorator if TESTING
        environment variable is not set and calls the plain function if it is set.
    """

    def decorator_wrapper(*decorator_args):
        def replacement(original_func):
            """Function that is called instead of original function."""

            def apply_guard(*args, **kwargs):
                """Decides whether to use decorator on function call."""
                if os.getenv("TESTING") is not None:
                    return original_func(*args, **kwargs)
                return decorator_func(*decorator_args)(original_func)(*args, **kwargs)

            return apply_guard

        return replacement

    return decorator_wrapper


def trace(fn):
    """Function decorator that enables tracing via stackdriver for performance
    metrics."""

    def trace_fn(*args, **kwargs):
        span = args[0].tracer
        if "GOOGLE_CLOUD_PROJECT" in os.environ:
            span_cm = span.span(name=f"{args[0].__class__.__name__} {fn.__name__}")
        # In the event where the environment variable is not set, use nullcontext
        # manager which does nothing
        else:
            span_cm = nullcontext()
        with span_cm:
            return fn(*args, **kwargs)

    return trace_fn


def before_send_to_sentry(event, hint):
    """Edit event properties before logging to Sentry

    Docs: https://docs.sentry.io/error-reporting/configuration/filtering/?platform=python#before-send
    """
    # Report a general logger name (`ingest_pipeline`) to Sentry,
    # rather than the function-specific logger name (e.g. `__main___errors`)
    # used internally within Ingest Pipeline.
    event["logger"] = "ingest_pipeline"
    return event


def integrate_sentry():
    """Log Ingest Pipeline errors to Sentry, by integrating with Python logger

    See also: links to Sentry resources atop this module
    """

    # Ultimately stored in Vault, passed in as environmen variable to PAPI
    sentry_DSN = os.environ.get("SENTRY_DSN")

    if sentry_DSN is None:
        # Don't log to Sentry unless its DSN is set.
        # This disables Sentry logging in development and test (i.e.,
        # environments without a SENTRY_DSN in their scp_config vault secret).
        return

    sentry_logging = LoggingIntegration(
        level=logging.ERROR,  # Capture error and above as breadcrumbs
        event_level=logging.ERROR,  # Send errors as events
    )
    sentry_sdk.init(
        dsn=sentry_DSN,
        integrations=[sentry_logging],
        attach_stacktrace=True,
        before_send=before_send_to_sentry,
    )


integrate_sentry()
