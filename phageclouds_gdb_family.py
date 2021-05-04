#!/usr/bin/env python

from neo4j import GraphDatabase
import pandas as pd
from collections import defaultdict
from pyvis import network as net
import sys, argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mplcol
from ete3 import NCBITaxa

class Neo4jConnection:
    def __init__(self, uri="bolt://127.0.0.1:7687", user='neo4j', pwd='phagedb'):
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            print("Failed to create the driver:", e)

    def close(self):
        if self.__driver is not None:
            self.__driver.close()

    def query(self, query):
        assert self.__driver is not None, "Driver not initialized!"
        session = None
        response = None
        print("CYPHER - running query:", query)
        try:
            session = self.__driver.session()
            response = list(session.run(query))
        except Exception as e:
            print("Query failed:", e)
        finally:
            if session is not None:
                session.close()
        return response

    def query2df(self, query):
        res = self.query(query)
        df = pd.DataFrame([dict(_) for _ in res])
        return df

def color_phage_tax(row, family, taxon):
    ncbi = NCBITaxa()
    taxid = ncbi.get_name_translator([family]).get(family)[0]
    descendant_taxa = ncbi.get_descendant_taxa(taxid)
    lineages_dict = ncbi.get_lineage_translator(descendant_taxa)
    phage_taxa_set = set()
    for key,value in lineages_dict.items():
        phage_taxa_set = phage_taxa_set.union({ncbi.get_taxid_translator([item]).get(item) for item in value if ncbi.get_rank([item]).get(item) == taxon})
    taxa_cmap = plt.cm.get_cmap('tab20', len(phage_taxa_set))
    color_dict = dict()
    for i,j in enumerate(phage_taxa_set):
        color_dict[j] = mplcol.to_hex(taxa_cmap(i))
    if row['Source'] == 'NCBI':
        for key,value in color_dict.items():
            if key in row['Lineage']:
                return value
        else:
            return '#000000'
    else:
        return '#FFFFFF'

def extract_phage_tax(row, taxon, taxdict):
    ncbi = NCBITaxa()
    taxid = taxdict.get(row['Phage'], None)
    if not taxid:
        return None
    else:
        taxid_lineage = ncbi.get_lineage(taxid)
        target_taxon = {ncbi.get_rank([item]).get(item):ncbi.get_taxid_translator([item]).get(item) for item in taxid_lineage}.get(taxon, None)
        return target_taxon



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract and draw phage clouds for user-defined phage family and colour nodes based on genus or subfamily')
    parser.add_argument('-f', '--fam', dest = 'family', help = 'Phage family to search', required = True)
    parser.add_argument('-t', '--tax', dest = 'taxon', help = 'Color nodes by "subfamily" or "genus" membership', required = True)
    parser.add_argument('-d', '--dist', dest = 'dist', type = float, help = 'User-defined distance threshold (default: 0.15)', default=0.15)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        phage_family = args.family
        phage_tax = args.taxon
        dist_thres = args.dist
        acc_taxid_file = 'phages.accessions.txt.taxid'
        conn = Neo4jConnection()
        query = """MATCH (p:PhageGenome {{source:"NCBI"}}) WHERE p.taxonomy CONTAINS "{}" OPTIONAL MATCH (p)-[r:sharesDNA]->(q:PhageGenome) WHERE r.distance <= {} WITH collect(p.accession) AS target_phages, collect(q.accession) AS connected_phages RETURN target_phages + [x IN connected_phages WHERE NOT x IN target_phages] AS phage_nodes;""".format(phage_family, dist_thres)
        target_phages = set(conn.query(query)[0].get('phage_nodes'))
        query = """MATCH (p:PhageGenome) WHERE p.accession IN {} RETURN p.accession AS Phage, p.source AS Source, p.genome_size AS Genome_size, p.taxonomy AS       Lineage;""".format(list(target_phages))
        nodes_df = conn.query2df(query)
        query = """MATCH (p:PhageGenome)-[r:sharesDNA]->(q:PhageGenome) WHERE p.accession IN {} AND q.accession IN {} AND r.distance <= {} RETURN p.accession AS Source, q.accession AS Target, r.distance AS Distance;""".format(list(target_phages), list(target_phages), dist_thres)
        edges_df = conn.query2df(query)
        nodes_df['Color'] = nodes_df.apply(color_phage_tax, args = (phage_family, phage_tax), axis = 1)
        with open(acc_taxid_file) as taxidfile:
            acc_taxid_dict = {line.rstrip().split(',')[0]:line.rstrip().split(',')[1] for line in taxidfile}
        nodes_df['Target_taxon'] = nodes_df.apply(extract_phage_tax, args = (phage_tax, acc_taxid_dict), axis = 1)
        size_scale_factor = 3000
        node_metadata = defaultdict(dict)
        for row in nodes_df.itertuples(index = False):
            node_size = int(row.Genome_size / size_scale_factor)
            colorsdict = {'border':'#000000', 'background':row.Color}
            node_metadata[row.Phage].update({'size':node_size, 'color':colorsdict, 'title':f'Source: {row.Source}<br>Genome size: {row.Genome_size:_}<br>{phage_tax}: {row.Target_taxon}'})
        pyvis_graph = net.Network(notebook = False, height='1500px', width='1500px')
        pyvis_graph.force_atlas_2based()
        for node, node_attrs in node_metadata.items():
            pyvis_graph.add_node(node, **node_attrs)
        edges_metadata = []
        for row in edges_df.itertuples(index = False):
            edges_metadata.append((row.Source, row.Target, {'weight':row.Distance, 'title':row.Distance}))
        for source,target,edge_attrs in edges_metadata:
            if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
                edge_attrs['value'] = dist_thres - edge_attrs['weight'] + 0.1
                edge_attrs['color'] = 'lightgray'
                pyvis_graph.add_edge(source, target, **edge_attrs)
        pyvis_graph.show_buttons()
        pyvis_graph.save_graph(f'{phage_family}_{"".join(str(dist_thres).split("."))}_{phage_tax}_clouds.html')
