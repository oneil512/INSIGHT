pubmed_api = """

API Examples

Here is an example showing how to find a list of ids of pubmed articles about breast cancer, and then get the abstracts of those studies.

===

from Bio import Entrez
import json

search_term = "breast cancer"

search_handle = Entrez.esearch(db="pubmed", term=search_term, retmax=5)
search_results = Entrez.read(search_handle)
search_handle.close()

pubmed_ids = search_results["IdList"]

fetch_handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="abstract", retmode="text")
abstracts = fetch_handle.read()
fetch_handle.close()

ret = abstracts

===

""".strip()