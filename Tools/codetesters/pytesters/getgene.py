#!/usr/bin/env python

import ensembl

serverRegistry = ensembl.get_registry(host='ensembldb.ensembl.org', user='anonymous')

coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')

geneAdaptor = coreDBAdaptor.get_adaptor('gene')

geneAdaptor.fetch_by_stable_id("
