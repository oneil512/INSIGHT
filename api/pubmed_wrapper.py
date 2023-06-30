def pubmed_wrapper(query_term, retmax, retstart):
    from Bio import Entrez

    search_handle = Entrez.esearch(db="pubmed", term=query_term, retmax=retmax, retstart=retstart)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    pubmed_ids = search_results["IdList"]

    fetch_handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="abstract")
    abstracts = fetch_handle.read()
    fetch_handle.close()

    return abstracts