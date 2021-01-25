# Service to record events to Mixpanel via Bard
#
# Mixpanel is a business analytics service for tracking user interactions.
# Bard is a DSP service that mediates writes to Mixpanel.
#

import json

import os
import requests


try:
    from monitor import setup_logger
except ImportError:
    from ..monitor import setup_logger


if os.environ.get("BARD_HOST_URL"):
    bard_host_url = os.environ["BARD_HOST_URL"]
else:
    bard_host_url = None


class MetricsService:
    # Logger provides more details
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")
    BARD_HOST_URL = bard_host_url

    # Log metrics to Mixpanel
    @classmethod
    def log(cls, event_name, props={}):
        properties = {"event": event_name, "properties": props.get_properties()}

        post_body = json.dumps(properties)
        MetricsService.post_event(post_body)

    # Post metrics to Mixpanel via Bard web service
    # Bard docs:
    # https://terra-bard-prod.appspot.com/docs/
    @staticmethod
    def post_event(props):
        try:
            r = requests.post(
                f"{MetricsService.BARD_HOST_URL}/api/event",
                headers={"content-type": "application/json"},
                data=props,
            )

            r.raise_for_status()
        # Don't want to stop parsing for logging errors. Errors will be logged and not raised.
        except requests.exceptions.HTTPError as e:
            #  401 Unauthorized
            MetricsService.dev_logger.exception(e)
        except requests.exceptions.RequestException as e:
            # Catastrophic error
            MetricsService.dev_logger.exception(e)
        except Exception as e:
            MetricsService.dev_logger.exception(e)
