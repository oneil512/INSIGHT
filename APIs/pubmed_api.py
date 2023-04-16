pubmed_api = """

API Examples

Here is an example showing how to find a list of ids of pubmed articles about breast cancer, and then get the abstracts of those studies.

The retmax parameter controls how many abstracts you will get back.
The retstart parameter controls where to start looking for results in the result set. This is useful if you have searched this query before and would like to get further results
The sort parameter can take [relevance|pub_date]

===

import json
from Bio import Entrez

search_term = "breast cancer"
retmax=6
retstart=0
sort='relevance'

search_handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax, retstart=retstart, sort=sort)
search_results = Entrez.read(search_handle)
search_handle.close()

pubmed_ids = search_results["IdList"]

fetch_handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="abstract", retmode="xml")
abstracts = fetch_handle.read()
fetch_handle.close()

ret = abstracts

===

""".strip()