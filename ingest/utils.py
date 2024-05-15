import os

def get_scp_api_origin(env):
    """Get domain etc. for SCP REST API URLs
    """
    db_name = os.environ['DATABASE_NAME']

    if env == None:
        env = db_name.split('_')[-1]

    origins_by_environment = {
        'development': 'https://localhost:3000',
        'staging': 'https://singlecell-staging.broadinstitute.org',
        'production': 'https://singlecell.broadinstitute.org'
    }

    return origins_by_environment[env]

def get_scp_api_base(env):
    """Get base URL for SCP API

    :param env override default, e.g. set to "production" in local ("development") host
    """
    origin = get_scp_api_origin(env)
    api_base = f"{origin}/single_cell/api/v1"
    return api_base


