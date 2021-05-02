#!/usr/bin/env python

from neo4j import GraphDatabase
import pandas as pd
from collections import defaultdict
from pyvis import network as net
import sys, argparse

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
    parser = argparse.ArgumentParser(description='Extract and draw phage clouds for user defined bacterial host')
    parser.add_argument('-g', '--genus', dest = 'host', help = 'host genus used for searching phage clouds', required = True)
    parser.add_argument('-t', '--thres', dest = 'dist', type = float, help = 'threshold applied to phage intergenomic distances (default: 0.25)', default = 0.25)
    parser.add_argument('--harsh', dest = 'harsh', action = 'store_true', help = 'exclude GTDB_predicted_prophages')
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        target_host = args.host
        dist_thres = args.dist
        conn = Neo4jConnection()
        if args.harsh:
            query = """MATCH (a:PhageGenome)-[r:sharesDNA]->(b:PhageGenome) WHERE (a)-[:infects]->(:Host {{genus:"{}"}}) AND r.distance <= {} AND a.source <> "GTDB_predicted_prophages" WITH collect(a.accession) as a_list, collect(b.accession) as b_list RETURN a_list + [x IN b_list WHERE NOT x IN a_list] AS node_list;""".format(target_host, dist_thres)
        else:
            query = """MATCH (a:PhageGenome)-[r:sharesDNA]->(b:PhageGenome) WHERE (a)-[:infects]->(:Host {{genus:"{}"}}) AND r.distance <= {} WITH collect(a.accession) as a_list, collect(b.accession) as b_list RETURN a_list + [x IN b_list WHERE NOT x IN a_list] AS node_list;""".format(target_host, dist_thres)
        node_set = set(conn.query(query)[0]['node_list'])
        query = """MATCH (a:PhageGenome) WHERE a.accession IN {} OPTIONAL MATCH (a)-[:infects]->(h:Host)
                RETURN a.accession AS Phage, a.source AS Source, a.genome_size AS Genome_size, a.genus AS Phage_genus, h.genus AS Host;""".format(list(node_set))
        nodes_df = conn.query2df(query)
        query = """MATCH (a:PhageGenome)-[r:sharesDNA]->(b:PhageGenome) WHERE a.accession IN {} AND b.accession IN {} AND r.distance <= {} RETURN a.accession as Start, b.accession as End, r.distance as Distance;""".format(list(node_set), list(node_set), dist_thres)
        edges_df = conn.query2df(query)
        node_color_dict = {'NCBI':'#8acb4a', 'Tara':'#39dede', 'GTDB_predicted_prophages':'#f1e653', 'GPD_Isolate':'#9b4aed', 'GPD_Metagenome':'#c734df'}
        node_metadata = defaultdict(dict)
        nodes_df['Color'] = nodes_df.apply(lambda x: node_color_dict[x['Source']], axis = 1)
        size_scale_factor = 3000
        for row in nodes_df.itertuples(index = False):
            node_size = int(row.Genome_size / size_scale_factor)
            colorsdict = {'border':'#000000', 'background':row.Color}
            phage_genus = row.Phage_genus
            host_genus = row.Host
            node_metadata[row.Phage].update({'size':node_size, 'color':colorsdict, 'title':f'Target host genus: {host_genus}<br>Phage genus: {phage_genus}<br>Genome size: {row.Genome_size:_} bp'})
        pyvis_graph = net.Network(notebook = False, height='1500px', width='1500px')
        pyvis_graph.force_atlas_2based()
        for node, node_attrs in node_metadata.items():
            pyvis_graph.add_node(node, **node_attrs)
        edges_metadata = []
        for row in edges_df.itertuples(index = False):
            edges_metadata.append((row.Start, row.End, {'weight':row.Distance, 'title':row.Distance}))
        for source,target,edge_attrs in edges_metadata:
            if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
                edge_attrs['value'] = dist_thres - edge_attrs['weight'] + 0.1
                edge_attrs['color'] = 'lightgray'
                pyvis_graph.add_edge(source, target, **edge_attrs)
        pyvis_graph.show_buttons()
        if args.harsh:
            pyvis_graph.save_graph(f'{target_host}_{"".join(str(dist_thres).split("."))}_clouds_v2_harsh.html')
        else:
            pyvis_graph.save_graph(f'{target_host}_{"".join(str(dist_thres).split("."))}_clouds_v2.html')

# CODE GRAVEYARD
        # query = """MATCH (a:PhageGenome) WHERE (a)-[:infects]->(:Host {{genus:"{}"}})
        #         RETURN a.accession AS Phage, a.source AS Source,
        #         a.genome_size AS Genome_size;""".format(target_host)
        # nodes_df = conn.query2df(query)
        # query = """MATCH (a:PhageGenome)-[r:sharesDNA]->(b:PhageGenome)
        #         WHERE (a)-[:infects]->(:Host {{genus:"{}"}}) AND
        #         (b)-[:infects]->(:Host {{genus:"{}"}}) AND r.distance <= {}
        #         RETURN a.accession as Start, b.accession as End,
        #         r.distance as Distance;""".format(target_host, target_host, dist_thres)
        # edges_df = conn.query2df(query)
        # node_metadata = defaultdict(dict)
        # def node_color(row_data):
        #     if row_data['Source'] == 'NCBI':
        #         return 'green'
        #     elif row_data['Source'] == 'GTDB_predicted_prophages':
        #         return 'yellow'
        #     else:
        #         return 'purple'
        # nodes_df['Color'] = nodes_df.apply(node_color, axis = 1)
        # size_scale_factor = 3000
        # for row in nodes_df.itertuples(index = False):
        #     node_size = int(row.Genome_size / size_scale_factor)
        #     colorsdict = {'border':'#000000', 'background':row.Color}
        #     node_metadata[row.Phage].update({'size':node_size, 'color':colorsdict})
        # pyvis_graph = net.Network(notebook = False, height='1500px', width='1500px')
        # pyvis_graph.force_atlas_2based()
        # for node, node_attrs in node_metadata.items():
        #     pyvis_graph.add_node(node, **node_attrs)
        # edges_metadata = []
        # for row in edges_df.itertuples(index = False):
        #     edges_metadata.append((row.Start, row.End, {'weight':row.Distance}))
        # for source,target,edge_attrs in edges_metadata:
        #     if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
        #         edge_attrs['value'] = dist_thres - edge_attrs['weight'] + 0.1
        #         edge_attrs['color'] = 'lightgray'
        #         pyvis_graph.add_edge(source, target, **edge_attrs)
        # pyvis_graph.show_buttons()
        # pyvis_graph.save_graph(f'{target_host}_{"".join(str(dist_thres).split("."))}_clouds.html')
