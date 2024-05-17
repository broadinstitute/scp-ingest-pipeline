"""Parse or infer scientific article publications for SCP studies

Many studies have publications but not all of them use canonical "Publication"
fields.  This module detects publications by text mining non-canonical sources.

A publication can reasonably be defined as "a work including substantial narrative scientific
prose that specifically contextualizes the study's data".
"""
import urllib

import ftfy
import requests

from utils import get_scp_api_base, fetch_pmid_pmcid

# Techniques:
# 1.  From the study itself, i.e. SCP API (/site/studies/{accession} endpoint)
#   A.  `publications` array -- ideal place for article cross-references
#       Introduced in Q1 2022, https://github.com/broadinstitute/single_cell_portal_core/pull/1379
#
#   B.  `external_resources` array -- not tailored for publications
#       Introduced in Q1 2019, https://github.com/broadinstitute/single_cell_portal_core/pull/221
#
#   C.  `description` string -- oldest, rawest place to note a publication
#       Introduced in Q3 2016 at SCP launch
#
# 2.  PubMed Central (PMC) API
#   A.  Search PubMed for references to SCP domain, see if study title or
#       authors roughly match publication title or authors.
#
# 3.  BiorXiv API
#   A.  Apply same technique as PMC API, but for BiorXiv API.
#
# Out of first 60 (visualizable) of 1633 studies in SCP as of 2024-05-16,
# 5-7 have directly detectable publications (i.e. via data in study itself):
# SCP2510* -> EX https://zenodo.org/records/10667499 -> https://www.nature.com/articles/s41590-024-01792-2
# SCP2484 -> DE https://doi.org/10.1093/bfgp/elac044 -> https://academic.oup.com/bfg/article/22/3/263/6874511
# SCP2454 -> DE https://doi.org/10.1101/2023.09.18.558077 -> https://www.biorxiv.org/content/10.1101/2023.09.18.558077v1
#     ^ Presumably same for 2-5 associated studies, linked from there (perhaps sometimes transitively)
# SCP2450 -> DE (private, abstract)
# SCP2393 -> PU, https://doi.org/doi:10.1126/sciadv.add9668, https://www.ncbi.nlm.nih.gov/pmc/articles/37756410
# SCP2384* -> ER https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244451 -> https://pubmed.ncbi.nlm.nih.gov/37873456 -> https://www.biorxiv.org/content/10.1101/2023.10.09.561581v1
# SCP2369 -> ER https://www.science.org/doi/10.1126/sciadv.adh9570
#
# * Zenodo and GEO links are not HTTP redirects to publications, like DOI links.
#
# This initial data suggests:
#   1.  ~10% of studies (5-7/60) have publications that are directly detectable.
#   2.  ~15-20% of study publications (1/5-7) use canonical "Publication" fields
#
# Benefits to refining non-canonical publication data:
#   1.  Ease findability of publication content -- it'd help to link in a consistent place and format
#   2.  Improve visual quality of study summary -- unconventional formatting can distract from exploration
#   3.  Increase findability of studies by (publication) author -- informal delimiters can confuse text search
#   4.  Enable text mining publication content, for e.g. relevant genes, canonicalizing custom abbreviations
#
# Benefits to inferring publications where study lacks even non-canonical publication data
#   1.  Enable finding publication -- linking from study where none was available would greatly lower barriers
#   2.  Plus all benefits of refining non-canonical publication data

biorxiv_url = "https://www.biorxiv.org"

env = None # Switch to "production" or "staging" to easily work with those locally

publication_bases = [
    "doi.org", "biorxiv.org", "nature.com", "science.org"
]

doi_stems_by_domain = {
    "nature.com": "10.1038"
}

publication_parseable_bases = [
    'zenodo.org', 'ncbi.nlm.nih.gov/geo'
]

def get_study(accession):
    """Get JSON object for study

    Docs: https://singlecell.broadinstitute.org/single_cell/api/v1#/Site/site_study_view_path
    """
    api_base = get_scp_api_base(env)
    study_api_url = f"{api_base}/site/studies/{accession}"
    response = requests.get(study_api_url)
    study_json = response.json()
    return study_json

def get_doi(url):
    """Convert a URL to a DOI (digital object identifier)

    Example DOI: 10.1093/bfgp/elac044
    """
    doi = None

    pub_base = next((pb for pb in publication_bases if pb in url), None)
    if pub_base is None:
        return None

    if pub_base == "doi.org":
        # E.g. https://doi.org/10.1093/bfgp/elac044 -> 10.1093/bfgp/elac044
        doi = url.split("https://doi.org/")[1]
    elif pub_base == "science.org":
        # E.g. https://www.science.org/doi/10.1126/sciadv.adh9570 -> 10.1126/sciadv.adh9570
        doi = url.split('/doi/')[1]
    elif pub_base in doi_stems_by_domain:
        doi_stem = doi_stems_by_domain[pub_base]
        last_url_segment = url.split("/")[-1]
        doi = f"{doi_stem}/{last_url_segment}"

    return doi

def fetch_citation(doi, bibliographic_data):
    """Get a ready-to-display citation string

    Example: (TODO)
    """
    # E.g. 10.1126/sciadv.adh9570 -> 10.1126%2Fsciadv.adh9570
    safe_doi = urllib.parse.quote_plus(doi)

    # Nature citation style is familiar to the single cell community,
    # and "no-et-al" provides the full list of authors, which
    # could be useful for SCP global search.
    style = "nature-no-et-al"

    # Example:
    # https://citation.crosscite.org/format?doi=10.1126%2Fsciadv.adh9570&style=nature-no-et-al&lang=en-US
    # UI: https://citation.crosscite.org
    crosscite_url = (
        "https://citation.crosscite.org/format" +
        f"?doi={safe_doi}&style={style}&lang=en-US"
    )
    crosscite_response = requests.get(crosscite_url)
    raw_cite = crosscite_response.text()

    # Fix occasional mangling of non-ASCII strings in Crosscite responses
    base_cite = ftfy.fix_encoding(raw_cite)

    # Get archive IDs; often shown for convenient bibliographic reference
    archives = []
    pmcid = bibliographic_data["pmcid"] # PubMed Central ID (full text)
    pmid = bibliographic_data["pmid"] # PubMed ID (abstracts)
    if pmcid:
        archives.append(f"PMCID: {pmcid}")
    if pmid:
        archives.append(f"PMID: {pmcid}")
    archive_ids.append(f"DOI: {doi}")
    archive_ids = archives.join("; ")

    citation = f"{base_cite} {archive_ids}"

    return citation

def fetch_bibliographic_data(doi):
    """Return publication title, journal, PMCID, authors, year, etc.

    (Abstract, references, and much else is also gathered, but not returned;
    these might be useful for future development.)
    """

    # Example:
    # https://api.crossref.org/works/10.1126/sciadv.adh9570/transform/application/vnd.citationstyles.csl+json
    crossref_url = (
        "https://api.crossref.org/works/" +
        f"{doi}/transform/application/vnd.citationstyles.csl+json"
    )
    crossref_response = requests.get(crossref_url)
    raw_biblio = crossref_response.json()
    date_ymd = raw_biblio["published"]["date-parts"][0] # e.g. "[2023, 8, 25]"

    # Get PubMed ID, and PubMed Central ID
    [pmid, pmcid] = fetch_pmid_pmcid(doi)

    bibliographic_data = {
        # Canonical publication fields
        "title": raw_biblio["title"], # required
        "journal": raw_biblio["container-title"], # required
        "pmcid": pmcid, # optional

        # Other useful fields
        "pmid": pmid,
        "journal_short": raw_biblio["container-title-short"],
        "year": date_ymd[0],

        # Example `authors`:
        # [{
        #    "ORCID": "http://orcid.org/0000-0002-5958-4616",
        #    "authenticated-orcid": true,
        #    "given": "Qijun",
        #    "family": "Tang",
        #    "sequence": "first",
        #    "affiliation": [{
        #       "name": "Department of Biology, University of Virginia, Charlottesville, VA 22904, USA."
        #       }]
        #     }, ...]
        "authors": raw_biblio["authors"]
    }
    return bibliographic_data

def parse_publications_from_resources(external_resources):
    """Return rich, canonical publication objects from resources, if possible
    """
    publications = []

    for external_resource in external_resources:
        url = external_resource["url"]
        # title = external_resource["title"]
        # description = external_resource["description"]
        doi = get_doi(url)

        if doi:
            # Lacking DOI means publication is pre-preprint or absent, so skip
            continue

        bibliographic_data = fetch_bibliographic_data(doi)
        citation = fetch_citation(doi, bibliographic_data)

        publication = {
            "title": bibliographic_data["title"], # required
            "journal": bibliographic_data["journal"], # required
            "url": url, # required

            "citation": citation, # optional, helps search
            "pmcid": bibliographic_data["pmcid"] # optional, helps text mining, search
        }
        publications.append(publication)

    return publications

def get_publications_for_study(accession):
    """Parse or infer publications corresponding to a study
    """
    study_json = get_study(accession)
    canonical_publications = study_json["publications"]
    external_resources = study_json["external_resources"]
    description = study_json["description"]

    publications = []

    # Canonical SCP publication format:
    #   Required: title, journal, URL
    #   Optional: pmcid, citation, preprint (Boolean, default: true)

    if len(canonical_publications) > 0:
        publications += canonical_publications
    if len(external_resources) > 0:
        publications += parse_publications_from_resources(external_resources)

accession = "SCP2560"



