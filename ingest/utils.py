import json
import os
import urllib

def get_scp_api_origin(env):
    """Get domain etc. for SCP REST API URLs
    """
    db_name = os.environ["DATABASE_NAME"]

    if env == None:
        env = db_name.split("_")[-1]

    origins_by_environment = {
        "development": "https://localhost:3000",
        "staging": "https://singlecell-staging.broadinstitute.org",
        "production": "https://singlecell.broadinstitute.org"
    }

    return origins_by_environment[env]

def get_scp_api_base(env):
    """Get base URL for SCP API

    :param env override default, e.g. set to "production" in local ("development") host
    """
    origin = get_scp_api_origin(env)
    api_base = f"{origin}/single_cell/api/v1"
    return api_base

def fetch_pmid_pmcid(doi, require_pmcid=False):
    """Convert Digital Object Identifier (DOI) into PMID and PMCID

    PubMed Central IDs (PMCIDs) map to open access, full-text articles
    PubMed IDs (PMIDs) map to open access abstracts

    PMCIDs are more useful than PMIDs, but there ~4x more PMIDs (37M vs 9.9M).

    Many PMC articles are completely public-access *and* machine-readable.
    """
    # Example:
    # https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?format=json&ids=10.1038/s41467-022-28372-y
    idconv_base = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
    params = "?" + "&".join(
        [
            "tool=scp-fetch-pmid-pmcid",
            "email=scp-dev@broadinstitute.org",
            "format=json",
            f"ids={doi}",
        ]
    )
    idconv_url = idconv_base + params
    with urllib.request.urlopen(idconv_url) as response:
        data = response.read().decode("utf-8")
    idconv_json = json.loads(data)
    record = idconv_json["records"][0]

    pmid = record["pmid"]
    pmcid = None

    if "pmcid" in record:
        pmcid = record["pmcid"]
    elif require_pmcid:
        msg = f"PubMed Central ID (PMCID) not found for DOI \"{doi}\""
        raise ValueError(msg)

    pmcid = record["pmcid"]
    return [pmid, pmcid]


