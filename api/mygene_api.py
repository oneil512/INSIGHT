mygene_api = """

RULES
    - You must decide the parameters to use
    - Forward slashes in the query must by escaped with a backwards slash:
        Instead of query_term="abc/123" do query_term="abc\/123"
    - The only thing you should change is the query_term, size, and from_

API PARAMETERS
    - size: How many genes you will get back
    - from_: Where in the results to start selecting from
    - query_term: What is being searched

API EXAMPLES

Here is an example showing the parameters used to query mygene for genes related to diabetes:

===

query_term="diabetes"
size=10
from_=0


""".strip()
