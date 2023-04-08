gene_api = """

API Rules
    - Every code block that calls to the mygene api must start with these two lines:
        import mygene
        mg = mygene.MyGeneInfo()
    - Every code block must end by assigning the output to a variable called 'ret'

API Examples

Here is an example showing how to find genes that are associated with diabetes and then find more information on those genes.

===

import mygene
mg = mygene.MyGeneInfo()

query_term = "diabetes"
gene_results = mg.query(query_term, fields="symbol,name,entrezgene,ensemblgene", species="human", size=10)
hits = gene_results['hits']

gene_info_list = []

for gene in hits:
    gene_info = mg.getgene(gene['_id'], fields=['name','symbol','type_of_gene','generif','genomic_pos_hg19','refseq','taxid', 'pathway'])
    gene_info_list.append(gene_info)

ret = gene_info_list

===

""".strip()