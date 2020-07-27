import logging
import os
import time
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


def setup_logger(logger_name, log_file, level=logging.DEBUG, format='default'):
    logger = logging.getLogger(logger_name)
    if format == 'default':
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s [%(study_id)s.%(name)s.%(funcName)s:%(lineno)d:%(duration)s] %(message)s'
        )
    else:
        formatter = logging.Formatter(format)
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def log(error_logger):
    error_logger = error_logger

    def debug(enter_message=None, exit_message=None):
        def wrapper(fn):
            def wrap(*args, **kwargs):
                study_id = args[0].study_id
                info_logger = setup_logger(args[0].__class__.__name__, 'info.txt')
                if enter_message is not None:
                    msg = enter_message
                else:
                    msg = f'Starting {fn.__name__}'
                try:
                    info_logger.info(
                        msg, extra={'duration': None, 'study_id': study_id}
                    )
                    start_time = time.time()
                    return_statements = fn(*args, **kwargs)  # running function
                    end_time = time.time()
                    exit_params = {
                        'duration': end_time - start_time,
                        'study_id': study_id,
                    }
                    if exit_message:
                        info_logger.info(exit_message, extra=exit_params)
                    else:
                        info_logger.info(f'Finished {fn.__name__}', extra=exit_params)
                    return return_statements
                except Exception as e:
                    error_logger.exception(e)

            return wrap

        return wrapper

    return debug


def trace(fn):
    """Function decorator that enables tracing via stackdriver for performance
    metrics."""

    def trace_fn(*args, **kwargs):
        span = args[0].tracer
        if 'GOOGLE_CLOUD_PROJECT' in os.environ:
            span_cm = span.span(name=f'{args[0].__class__.__name__} {fn.__name__}')
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
    event['logger'] = 'ingest_pipeline'
    return event


def integrate_sentry():
    '''Log Ingest Pipeline errors to Sentry, by integrating with Python logger

    See also: links to Sentry resources atop this module
    '''

    # Ultimately stored in Vault, passed in as environmen variable to PAPI
    sentry_DSN = os.environ.get('SENTRY_DSN')

    if sentry_DSN is None:
        # Don't log to Sentry unless its DSN is set.
        # This disables Sentry logging in development and test (i.e.,
        # environments without a SENTRY_DSN in their scp_config vault secret).
        return

    sentry_logging = LoggingIntegration(
        level=logging.ERROR,       # Capture error and above as breadcrumbs
        event_level=logging.ERROR  # Send errors as events
    )
    sentry_sdk.init(
        dsn=sentry_DSN,
        integrations=[sentry_logging],
        attach_stacktrace=True,
        before_send=before_send_to_sentry
    )


integrate_sentry()
