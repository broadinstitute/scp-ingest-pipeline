# Service to record events to Mixpanel via Bard
#
# Mixpanel is a business analytics service for tracking user interactions.
# Bard is a DSP service that mediates writes to Mixpanel.
#

import uuid
import json
import requests


class MetricsService:
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
    @classmethod
    def post_event(cls, props):
        requests.post(MetricsService.BARD_HOST_URL, json=props)
