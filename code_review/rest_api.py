import json
import urllib.request

import pprint

# favourite gene of Rebecca
queryfield = "symbol"
query = "zeb1"
returnfields = "entrezgene,ensembl.gene"

# # favourite gene of jos
# queryfield = "ensembl.gene"
# query = "ENSG00000007372"
# returnfields = "HGNC,symbol,entrezgene"

request = f"http://mygene.info/v3/query?q={queryfield}:{query}&fields={returnfields}"
print(request)


with urllib.request.urlopen(request) as url:
    data = json.loads(url.read().decode())
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(data)
