pubmed_api = """

RULES
    - The only thing you should change is the query_term, retmax, and retstart

API PARAMETERS
    - query_term: What is being searched
    - retmax: Max results to return
    - retstart: Where to start searching from

API EXAMPLES

Here is an example showing how to find a list of 10 ids of pubmed articles about breast cancer.

===

query_term = 'breast cancer'
retmax = 10
retstart = 0


""".strip()
