myvariant_api = """

API RULES
    - Every code block that calls to the mygene api must start with these two lines:
        import myvariant
        mv = myvariant.MyVariantInfo()
    - Every code block must end by assigning the output to a variable called 'ret'
    - Do not change the variable names.
    - The only thing you should change is the query_term

API PARAMETERS
    - query_term: What is being searched

API EXAMPLES

Here is an example showing how to find information on the variant chr7:g.140453134T>C

Note that chr7:g.140453134T>C" describes a variant on chromosome 7 at position 140453134, where the usual thymine (T) is replaced by a cytosine (C).

You can also search with an rsID or ClinVar

===

import myvariant
mv = myvariant.MyVariantInfo()

query_term = 'chr1:g.35367G>A'

ret = mv.getvariant(query_term)

===

""".strip()
