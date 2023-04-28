mygene_api = """

API RULES
    - Every code block that calls to the mygene api must start with these two lines:
        import mygene
        mg = mygene.MyGeneInfo()
    - Every code block must end by assigning the output to a variable called 'ret'
    - Forward slashes in the query must by escaped with a backwards slash:
        Instead of query_term="abc/123" do query_term="abc\/123"
    - Only use the fields you see in the examples. Do not create any of your own fields
    - Do not change the variable names.
    - The only thing you should change is the query_term, size, and from_

API PARAMETERS
    - size: How many genes you will get back
    - from_: Where in the results to start selecting from
    - query_term: What is being searched

API EXAMPLES

Here is an example showing how to find genes that are associated with diabetes and then find more information on those genes.

===

import mygene
mg = mygene.MyGeneInfo()

query_term="diabetes"
fields="symbol,name,entrezgene,ensemblgene"
size=10
from_=0

gene_results = mg.query(query_term, fields=fields, species="human", size=size, from_=from_)
hits = gene_results['hits']

gene_info_list = []

for gene in hits:
    gene_info = mg.getgene(gene['_id'], fields=['name','symbol','type_of_gene','genomic_pos_hg19','refseq','taxid','summary','pathway'])
    gene_info_list.append(gene_info)

ret = gene_info_list

===

""".strip()
