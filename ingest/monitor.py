import logging
import os
import time
from functools import wraps
from contextlib import nullcontext

import memory_profiler


def setup_logger(logger_name, log_file, level=logging.DEBUG):
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s [%(study_id)s.%(name)s.%(funcName)s:%(lineno)d:%(duration)s] %(message)s'
    )
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger.setLevel(level)
    logger.addHandler(handler)
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

fp = open('profile.log', 'w+')


def profile(fn):
    """Profiles memory usage in functions that apply this as a decorator
    """
    memory_profiler.profile(func=fn, stream=fp)
    return fn
    # @wraps(fn)
    # def profile_fn(*args, **kwargs):
    #     inner_self = args[0].__dict__
    #     print('inner_self')
    #     print(inner_self)
    #     if not inner_self['profile_memory']:
    #         return fn(*args, **kwargs)
    #     else:
    #         val = memory_profiler.profile(stream=fp)
    #         print('val')
    #         print(val)
    #         print('val.__repr__')
    #         print(val.__repr__)
    #         return fn(*args, **kwargs)
        
    # return profile_fn
