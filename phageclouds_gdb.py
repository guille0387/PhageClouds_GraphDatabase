#!/usr/bin/env python

from neo4j import GraphDatabase
import pandas as pd
from collections import defaultdict
from pyvis import network as net
import argparse, sys

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract and draw phage clouds for user defined phage taxon')
    parser.add_argument('-t', '--tax', dest = 'taxon', help = 'Phage taxon to search', required = True)
    parser.add_argument('-d', '--dist', dest = 'dist', type = float, help = 'Distance threshold (default : 0.25)', default = 0.25)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        phage_tax = args.taxon
        dist_thres = args.dist
        conn = Neo4jConnection()

        query = """MATCH (a:PhageGenome {{source:'NCBI'}})-[r:sharesDNA]->(b:PhageGenome)
                WHERE a.taxonomy CONTAINS '{}' AND r.distance <= {}
                RETURN a.accession AS {}_phage, b.accession AS target_phage;""".format(phage_tax, dist_thres, phage_tax)

        result_df = conn.query2df(query)

        node_set = set(result_df[f'{phage_tax}_phage'].values).union(set(result_df['target_phage'].values))

        query = """MATCH (a:PhageGenome)-[r:sharesDNA]->(b:PhageGenome)
                WHERE a.accession in {} AND b.accession in {} AND r.distance <= {}
                RETURN a.accession AS Source, b.accession AS Target,
                r.distance as Distance;""".format(list(node_set), list(node_set), dist_thres)

        edges_df = conn.query2df(query)

        query = """MATCH (a:PhageGenome) WHERE a.accession in {}
                RETURN a.accession as Phage, a.source as Source, a.genome_size as Genome_size,
                a.taxonomy CONTAINS '{}' as Phage_is_{};""".format(list(node_set), phage_tax, phage_tax)

        nodes_df = conn.query2df(query)

        node_metadata = defaultdict(dict)

        def node_color(row_data):
            if row_data['Source'] == 'NCBI' and row_data[f'Phage_is_{phage_tax}']:
                return 'green'
            elif row_data['Source'] == 'NCBI' and not row_data[f'Phage_is_{phage_tax}']:
                return 'red'
            elif row_data['Source'] == 'Tara':
                return 'cyan'
            elif row_data['Source'] == 'GPD_Isolate':
                return 'pink'
            elif row_data['Source'] == 'GPD_Metagenome':
                return 'purple'
            else:
                return 'yellow'

        nodes_df['color'] = nodes_df.apply(node_color, axis = 1)

        size_scale_factor = 3000

        for row in nodes_df.itertuples(index = False):
            node_size = int(row.Genome_size / size_scale_factor)
            colorsdict = {'border':'#000000', 'background':row.color}
            node_metadata[row.Phage].update({'size':node_size, 'color':colorsdict})

        pyvis_graph = net.Network(notebook = False, height='1500px', width='1500px')

        for node, node_attrs in node_metadata.items():
            pyvis_graph.add_node(node, **node_attrs)

        edges_metadata = []

        for row in edges_df.itertuples(index = False):
            edges_metadata.append((row.Source, row.Target, {'weight':row.Distance}))

        for source,target,edge_attrs in edges_metadata:
            if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
                edge_attrs['value'] = dist_thres - edge_attrs['weight'] + 0.1
                edge_attrs['color'] = 'lightgray'
                pyvis_graph.add_edge(source, target, **edge_attrs)

        pyvis_graph.show_buttons()

        pyvis_graph.save_graph(f'{phage_tax}_{"".join(str(dist_thres).split("."))}_clouds.html')
