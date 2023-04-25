pubmed_api = """

API RULES
    - Every code block that calls to the mygene api must start with this line:
        from Bio import Entrez
    - Every code block must end by assigning the output to a variable called 'ret'

API PARAMETERS
    - query_term: What is being searched
    - retmax: Max results to return
    - retstart: Where to start searching from
    - sort: How to sort results

API EXAMPLES

Here is an example showing how to find a list of ids of pubmed articles about breast cancer, and then get the abstracts of those studies.

===

from Bio import Entrez

query_term = "breast cancer"
retmax=6
retstart=0
sort='relevance'

search_handle = Entrez.esearch(db="pubmed", term=query_term, retmax=retmax, retstart=retstart, sort=sort)
search_results = Entrez.read(search_handle)
search_handle.close()

pubmed_ids = search_results["IdList"]

fetch_handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="abstract", retmode="xml")
abstracts = fetch_handle.read()
fetch_handle.close()

ret = abstracts

===

""".strip()
