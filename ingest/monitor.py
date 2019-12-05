import logging
import time


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