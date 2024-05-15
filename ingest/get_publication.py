"""Get scientific article publication for SCP study, via SCP API or inference

Many studies have publications, but not all of them are explicitly specified
using SCP's "Publication" fields.  This module fetches publications through
a series of techniques:

1.  From the study itself, i.e. SCP API (/site/studies/{accession} endpoint)
  A.  `publications` array -- ideal place for article cross-references
      Introduced in Q1 2022, https://github.com/broadinstitute/single_cell_portal_core/pull/1379

  B.  `external_resources` array -- not tailored for publications
      Introduced in Q1 2019, https://github.com/broadinstitute/single_cell_portal_core/pull/221

  C.  `description` string -- oldest, rawest place to note a publication
      Introduced in Q3 2016 at SCP launch

2.  PubMed Central (PMC) API
  A.  Search PubMed for references to SCP domain, see if study title or
      authors roughly match publication title or authors.

3.  BiorXiv API
  A.  Apply same technique as PMC API, but for BiorXiv API.
"""

import requests

from utils import get_scp_api_base

biorxiv_url = "https://www.biorxiv.org"

env = None # Switch to "production" or "staging" to easily work with those locally

def get_study(accession):
    """Get JSON object for study
    """
    api_base = get_scp_api_base(env)
    study_api_url = f"{api_base}/site/studies/{accession}"
    response = requests.get(study_api_url)
    study_json = response.json()
    return study_json

def get_publication_from_study(accession):
    study_json = get_study(accession)
    publications = study_json["publications"]
    external_resources = study_json["external_resources"]
    description = study_json["description"]


    if len(publications) > 0:

    elif len(external_resources) > 0:

    else:


accession = "SCP2560"



