myvariant_api = """

RULES
    - Decide what query_term to use. It could be an rsID, a ClinVar ID, or a variant in the format chr1:g.35367G>A

API PARAMETERS
    - query_term: What is being searched

API EXAMPLES

Here is an example showing how to find information on the variant chr7:g.140453134T>C

Note that chr7:g.140453134T>C" describes a variant on chromosome 7 at position 140453134, where the usual thymine (T) is replaced by a cytosine (C).

You can also search with an rsID or ClinVar

===

query_term = 'chr1:g.35367G>A'


""".strip()
