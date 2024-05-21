"""Parse or infer scientific article publications for SCP studies

Many studies have publications but not all of them use canonical "Publication"
fields.  This module detects publications by text mining non-canonical sources.

A publication can reasonably be defined as "a work including substantial narrative scientific
prose that specifically contextualizes the study's data".
"""
import json
from pathlib import Path
import re
import time
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
# TODO:
# 2.  PubMed Central (PMC) API
#   A.  Search PubMed for references to SCP domain, see if study title or
#       authors roughly match publication title or authors.
#
# 3.  BiorXiv API
#   A.  Apply same technique as PMC API, but for BiorXiv API.
#
# Out of first 60 (visualizable) of 1633 studies in SCP as of 2024-05-16,
# 5-7 have directly detectable publications (i.e. via data in study itself):
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
    print(f"Get study object from from SCP API for {accession}")
    api_base = get_scp_api_base(env)
    study_api_url = f"{api_base}/site/studies/{accession}"
    response = requests.get(study_api_url, verify=False)
    study_json = response.json()
    return study_json

def get_doi(url):
    """Convert a URL to a DOI (digital object identifier)

    Example DOI: 10.1093/bfgp/elac044
    """
    doi = None
    print(f"Get DOI for URL: {url}")

    pub_base = next((pb for pb in publication_bases if pb in url), None)
    if pub_base is None:
        print("Cannot determine DOI from publication URL")
        return None

    if pub_base == "doi.org":
        # E.g. https://doi.org/10.1093/bfgp/elac044 -> 10.1093/bfgp/elac044
        doi = url.split("doi.org/")[1]
    elif pub_base == "science.org":
        # E.g. https://www.science.org/doi/10.1126/sciadv.adh9570 -> 10.1126/sciadv.adh9570
        doi = url.split('/doi/')[1]
    elif pub_base in doi_stems_by_domain:
        doi_stem = doi_stems_by_domain[pub_base]
        if "nature.com" in url and "/figures/" in url:
            # Don't attempt to get DOIs for images and their captions
            # https://www.nature.com/articles/s42003-022-04196-w/figures/1
            return None
        last_url_segment = url.split("/")[-1]
        last_url_segment = last_url_segment.replace(".html", "")
        doi = f"{doi_stem}/{last_url_segment}"

    print("Fetched DOI", doi)
    return doi

def fetch_citation(doi, bibliographic_data):
    """Get a ready-to-display citation string

    Example (for SCP2369, DOI 10.1126/sciadv.adh9570):
    Tang, Q., Godschall, E., Brennan, C. D., Zhang, Q., Abraham-Fan, R.-J., Williams, S. P., Güngül, T. B., Onoharigho, R., Buyukaksakal, A., Salinas, R., Sajonia, I. R., Olivieri, J. J., Calhan, O. Y., Deppmann, C. D., Campbell, J. N., Podyma, B. & Güler, A. D. Leptin receptor neurons in the dorsomedial hypothalamus input to the circadian feeding network. Sci. Adv. 9, (2023). PMCID: PMC10456850; PMID: 37624889; DOI: 10.1126/sciadv.adh9570

    TODO: Remove HTML from citation, e.g. "<i>in vivo</i>" for SCP2454
    """
    print("Fetch citation for DOI", doi)
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
    status_code = crosscite_response.status_code
    if status_code == 500:
        # Occurs with e.g. bioRxiv publications
        return ""

    raw_cite = crosscite_response.text.strip()

    # Fix occasional mangling of non-ASCII strings in Crosscite responses
    base_cite = ftfy.fix_encoding(raw_cite)

    if base_cite.startswith("1. "):
        # Fixes e.g. "1. Sundell, T., Grimstad, K." for SCP2484
        base_cite = base_cite.replace("1. ", "", 1)

    if " . . . " in base_cite:
        # Fix "Shannon, E., . . . Murphy, G. J." for SCP2393
        base_cite = base_cite.replace(" . . . ", "", 1)

    # Get archive IDs; often shown for convenient bibliographic reference
    archives = []
    pmcid = bibliographic_data["pmcid"] # PubMed Central ID (full text)
    pmid = bibliographic_data["pmid"] # PubMed ID (abstracts)
    if pmcid:
        archives.append(f"PMCID: {pmcid}")
    if pmid:
        archives.append(f"PMID: {pmid}")
    archives.append(f"DOI: {doi}")
    archive_ids = "; ".join(archives)

    citation = f"{base_cite} {archive_ids}"

    return citation

def fetch_bibliographic_data(doi):
    """Return publication title, journal, PMCID, authors, year, etc.

    (Abstract, references, and much else is also gathered, but not returned;
    these might be useful for future development.)
    """
    print("Fetch bibliographic data for DOI", doi)
    # Example:
    # https://api.crossref.org/works/10.1126/sciadv.adh9570/transform/application/vnd.citationstyles.csl+json
    crossref_url = (
        "https://api.crossref.org/works/" +
        f"{doi}/transform/application/vnd.citationstyles.csl+json" +
        "?mailto=scp-dev@broadinstitute.org"
    )
    crossref_response = requests.get(crossref_url)
    raw_biblio = crossref_response.json()

    date_ymd = raw_biblio["published"]["date-parts"][0] # e.g. "[2023, 8, 25]"

    # Get PubMed ID, and PubMed Central ID
    [pmid, pmcid] = fetch_pmid_pmcid(doi)

    journal = raw_biblio["container-title"]
    if journal == []:
        # e.g. "bioRxiv" for DOI "10.1101/2023.09.18.558077"
        journal = raw_biblio["institution"][0]["name"]
    journal_short = raw_biblio.get("container-title-short", journal)

    bibliographic_data = {
        # Canonical publication fields
        "title": raw_biblio["title"], # required
        "journal": journal_short, # required
        "pmcid": pmcid, # optional

        # Other useful fields
        "pmid": pmid,
        "journal_long": journal,
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
        "authors": raw_biblio["author"]
    }
    return bibliographic_data

def parse_publication_from_url(url):
    """Return rich, canonical publication objects from URL, if possible
    """
    doi = get_doi(url)

    if not doi:
        # Lacking DOI means publication is pre-preprint or absent
        return None

    bibliographic_data = fetch_bibliographic_data(doi)
    citation = fetch_citation(doi, bibliographic_data)

    publication = {
        "title": bibliographic_data["title"], # required
        "journal": bibliographic_data["journal"], # required
        "url": url, # required

        "citation": citation, # optional, helps search
        "pmcid": bibliographic_data["pmcid"] # optional, helps text mining, search
    }
    return publication

def parse_publications_from_resources(external_resources):
    """Return publication objects for resources' URLs, if possible
    """
    print("Parse publications from resources")
    publications = []

    for external_resource in external_resources:
        url = external_resource["url"]
        publication = parse_publication_from_url(url)
        if publication:
            publications.append(publication)

    return publications

def parse_publications_from_description(description):
    """Return publication objects for free-text description, if possible
    """
    print("Parse publications from description")
    publications = []

    # From https://stackoverflow.com/a/15518253
    url_regex = '<a\s*href=[\'|"](.*?)[\'"].*?>'
    urls = re.findall(url_regex, description)

    for url in urls:
        publication = parse_publication_from_url(url)
        if publication:
            publications.append(publication)

    return publications

def get_publications_for_study(accession, study_json=None):
    """Parse or infer publications corresponding to a study
    """
    if not study_json:
        study_json = get_study(accession)
    canonical_publications = study_json["publications"]
    external_resources = study_json["external_resources"]
    description = study_json["full_description"]

    publications = []

    if len(canonical_publications) > 0:
        publications += canonical_publications
    if len(external_resources) > 0:
        publications += parse_publications_from_resources(external_resources)
    if len(publications) == 0:
        publications += parse_publications_from_description(description)

    return publications

def fetch_and_merge_study(accession, studies_json, studies_json_path):
    time.sleep(2)
    study_json = get_study(accession)
    studies_json.append(study_json)

    with open(studies_json_path, "w") as f:
        f.write(json.dumps(studies_json))
    print(f"Wrote study {accession} to {studies_json_path}")

    return studies_json

def get_public_study_objects(reuse_studies_json=False):
    """Get JSON objects for all public studies

    Docs: https://singlecell.broadinstitute.org/single_cell/api/v1#/Site/site_study_view_path
    """
    print(f"Get all public study objects from from SCP API")

    api_base = get_scp_api_base(env)
    study_api_url = f"{api_base}/site/studies"
    response = requests.get(study_api_url, verify=False)
    studies_brief_json = response.json()

    if not reuse_studies_json:
        studies_json = []

        for study_brief in studies_brief_json:
            accession = study_brief["accession"]
            studies_json = fetch_and_merge_study(
                accession, studies_json, studies_json_path
            )
    else:
        print(f"Reusing cached public studies JSON from: {studies_json_path}")

        with open(studies_json_path) as f:
            studies_json = json.loads(f.read())

        cached_accessions = [s["accession"] for s in studies_json]

        for study_brief in studies_brief_json:
            accession = study_brief["accession"]
            if accession in cached_accessions:
                print(f"Using cached data for study {accession}")
                continue
            studies_json = fetch_and_merge_study(
                accession, studies_json, studies_json_path
            )


    print(f"Finished writing public studies to {studies_json_path}")

    return studies_json

def fetch_and_merge_publications(
    accession, study_json, publications, publications_json_path
):
    study_publications = get_publications_for_study(accession, study_json)

    with open(publications_json_path, "w") as f:
        publications.append({
            "accession": accession,
            "publications": study_publications,
        })
        f.write(json.dumps(publications))
    print(f"Wrote study publications for {accession} to {publications_json_path}")

    return publications

def get_all_study_publications(studies_json, reuse_publications_json=False):
    publications = [] # publications for all studies
    publications_json = []
    publications_json_path = "publications.json"
    cache_file = Path(publications_json_path)
    if not cache_file.exists():
        reuse_publications_json = False
    else:
        with open(publications_json_path) as f:
            publications_json = json.loads(f.read())

    print('publications_json', publications_json)
    cached_accessions = [p["accession"] for p in publications_json]

    for study_json in studies_json:
        # Get publications for this study
        accession = study_json["accession"]
        if not reuse_publications_json or accession not in cached_accessions:
            publications = fetch_and_merge_publications(
                accession, study_json, publications, publications_json_path
            )
        else:
            print(f"Using cached publications for {accession}")
            # Use cached study publications
            study_publications = next(
                p for p in publications_json if p["accession"] == accession
            )
            publications.append(study_publications)

    print(f"Fetched all publications for all studies; wrote to {publications_json_path}")
    return publications

# SCP2510* -> EXT https://zenodo.org/records/10667499 -> https://www.nature.com/articles/s41590-024-01792-2
# SCP2484 -> DES https://doi.org/10.1093/bfgp/elac044 -> https://academic.oup.com/bfg/article/22/3/263/6874511
# SCP2454 -> DES https://doi.org/10.1101/2023.09.18.558077 -> https://www.biorxiv.org/content/10.1101/2023.09.18.558077v1
#     ^ Presumably same for 2-5 associated studies, linked from there (perhaps sometimes transitively)
# SCP2450** -> DES (private, abstract)
# SCP2393 -> PUB https://doi.org/doi:10.1126/sciadv.add9668, https://www.ncbi.nlm.nih.gov/pmc/articles/37756410
# SCP2384* -> EXT https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244451 -> https://pubmed.ncbi.nlm.nih.gov/37873456 -> https://www.biorxiv.org/content/10.1101/2023.10.09.561581v1
# SCP2369 -> EXT https://www.science.org/doi/10.1126/sciadv.adh9570
#
# * Zenodo and GEO links are not HTTP redirects to publications, like DOI links.
# ** Private studies can be supported via gcloud auth

accession = "SCP2454"
#env = None # Switch to "production" or "staging" to easily work with those locally
env = "production"

publications = []

studies_json_path = "studies.json"
reuse_studies_json = True # Use previously-cached studies
reuse_publications_json = True # Use previously-cached publications

def run():
    studies_json = get_public_study_objects(reuse_studies_json)
    publications = get_all_study_publications(studies_json)

    print('publications')
    print(publications)


