# Service to record events to Mixpanel via Bard
#
# Mixpanel is a business analytics service for tracking user interactions.
# Bard is a DSP service that mediates writes to Mixpanel.
#

import json

# import os
import requests

try:
    from monitor import setup_logger
except ImportError:
    from ..monitor import setup_logger


class MetricsService:
    # Logger provides more details
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")
    BARD_HOST_URL = "https://terra-bard-dev.appspot.com/api/event"
    user_id = "2f30ec50-a04d-4d43-8fd1-b136a2045079"

    @classmethod
    def log(cls, event_name, props={}):
        props["distinct_id"] = MetricsService.user_id
        properties = {"event": event_name, "properties": props}

        post_body = json.dumps(properties)
        MetricsService.post_event(post_body)

    # Log metrics to Mixpanel via Bard web service
    #
    # Bard docs:
    # https://terra-bard-prod.appspot.com/docs/
    @staticmethod
    def post_event(props):
        try:
            r = requests.post(
                MetricsService.BARD_HOST_URL,
                headers={"content-type": "application/json"},
                data=props,
            )
            r.raise_for_status()
        # Don't want to stop parsing for logging errors. Errors will be logged and not raised.
        except requests.exceptions.HTTPError as e:
            MetricsService.dev_logger.exception(e)
            #  401 Unauthorized
        except requests.exceptions.RequestException as e:
            # Catastrophic error
            MetricsService.dev_logger.critcal(e)
        except Exception as e:
            MetricsService.dev_logger.exception(e)
