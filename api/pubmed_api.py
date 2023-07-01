pubmed_api = """

RULES
    - The only thing you should change is the query_term, retmax, and retstart
    - Prioritize querying relavent clinical trials.
    - query_term must be under 4 words. Strive to only include the most salient words that capture the task. For example, for a task: Find studies on the safety and efficacy of revumenib (SNDX-5613) in patients with relapsed or refractory acute leukemia, a good query_term might be: "revumenib leukemia". This is because it is the shortest query that captures the main point of the task.
    - Long query_term usually produce no results.
    - Avoid any abbreviations.

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
