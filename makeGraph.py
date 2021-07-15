import pickle
import networkx as nx
import math
import pandas as pd

'''
Get the parent number (move UP the tree one level/generation).
'''
def up(n):
    sep = '.'
    n = n.split(sep)
    if len(n) > 1:
        n.pop()
        n = sep.join(n)
    else:
        return n[0][0]
    return n


'''
Get indications headings.
'''
def get_headings(drug):
    if 'CHEMBL' in drug:
        indications = chembl[chembl['chembl_id'] == drug.upper()]
    else:
        indications = chembl[chembl['pref_name'] == drug.upper()]
    headings = sorted(list(set(indications.mesh_heading)))
    return headings


'''
Make a graph where the nodes are MeSH headings and the directed edges form the heirarchy. Each node has as an attribute the list of drugs that include that node.
'''
def make_graph(drugs):

    node_drug_dict = {}

    # initialize DAG
    G = nx.DiGraph()

    for drug in drugs:
        drug_node_dict[drug] = set()
        for heading in get_headings(drug):
            # add heading node
            if heading not in list(G.nodes):
                G.add_node(heading)
                node_drug_dict[heading] = set()
                
            drug_node_dict[drug].add(heading)
            node_drug_dict[heading].add(drug)

            # add all parents
            for n in mesh_headings[heading]:
                c_heading = heading
                for i in range(n.count('.') + 1):
                    p = up(n) # get parent number
                    p_heading = mesh_numbers[p] # get parent heading

                    # add parent heading node
                    if p_heading not in list(G.nodes):
                        G.add_node(p_heading)
                        node_drug_dict[p_heading] = set()
                    
                    drug_node_dict[drug].add(p_heading)
                    node_drug_dict[p_heading].add(drug)
                    
                    if not G.has_edge(p_heading, c_heading):
                        G.add_edge(p_heading, c_heading) # add directed edge from parent to child

                    n = p # move up path
                    c_heading = p_heading
    
    nx.set_node_attributes(G, node_drug_dict, 'drugs')

    return G


'''
Compute information accretion for each node in a graph.
'''
def compute_ia(G):
    # annotate nodes with information accretion value from probability
    node_ia_dict = {}

    for node in G.nodes:
        n_drugs = len(G.nodes[node]['drugs']) # get the number of drugs with this heading 

        ## get the number of drugs with all parent headings

        parents = G.predecessors(node) # get parent nodes

        # get drugs for each parent
        drug_sets = []
        for parent in parents:
            drug_sets.append(G.nodes[parent]['drugs'])

        n_p_drugs = len(set().union(*drug_sets)) # count the number of drugs that have all parent headings

        # compute probability
        if n_p_drugs == 0:
            prob = 1
        else:
            prob = n_drugs / n_p_drugs

        # compute information accretion
        ia = -math.log2(prob) # does this work for single values or does it have to be an array?

        node_ia_dict[node] = ia

    nx.set_node_attributes(G, node_ia_dict, 'ia')

    return G


'''
Main.
'''
if __name__ == "__main__":

    # load in MeSH pickle files
    f = open("data/mesh_headings.pkl", "rb")
    mesh_headings = pickle.load(f)
    f = open("data/mesh_numbers.pkl", "rb")
    mesh_numbers= pickle.load(f)

    # load in ChEMBL
    chembl = pd.read_csv('data/chembl_indications.tsv', sep='\t')

    drugs = list(set(chembl['pref_name'])) # get drug list

    drug_node_dict = {}

    G = make_graph(drugs) # build graph
    G = compute_ia(G) # compute and add ia values

    # save graph to pickle file
    nx.write_gpickle(G, "data/chembl.gpickle")

    # save drug_node_dict to pickle
    dict_file = open ('data/drug_node_dict.pickle', 'wb')
    pickle.dump(drug_node_dict, dict_file)