
def mygene_wrapper(query_term, size, from_):
    import mygene
    mg = mygene.MyGeneInfo()
    fields="symbol,name,entrezgene,ensemblgene"


    gene_results = mg.query(query_term, fields=fields, species="human", size=size, from_=from_)
    hits = gene_results['hits']

    gene_info_list = []

    for gene in hits:
        gene_info = mg.getgene(gene['_id'], fields=['name','symbol','type_of_gene','genomic_pos_hg19','refseq','taxid','generif','summary','pathway'])
        gene_info_list.append(gene_info)

    return gene_info_list