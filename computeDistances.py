import numpy as np
import argparse
import pickle
import networkx as nx
import math
from itertools import combinations, permutations
import pandas as pd
from multiprocessing import Pool, Manager


'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Compute the semantic distance between drugs in ChEMBL using the MeSH headings of their indications.')
    parser.add_argument('-n', '--num-cpus', help="The number of cpus to use. Default is 4.", default=4, type=int, required=False)
    parser.add_argument('-r', '--rows', help="The number of rows of ChEMBL indications to use. Default is all.", type=int, required=False)
    args = parser.parse_args()
    return args


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
Calculate mis-information between two drug graphs.
'''
def compute_mi(drug1, drug2):
    nodes = np.setdiff1d(list(drug_node_dict[drug1]), list(drug_node_dict[drug2]))

    mi = 0
    for n in nodes:
        mi += G.nodes[n]['ia']

    return mi


'''
Compute remaining uncertainty between two drug graphs.
'''
def compute_ru(drug1, drug2):
    nodes = np.setdiff1d(list(drug_node_dict[drug2]), list(drug_node_dict[drug1]))

    ru = 0
    for n in nodes:
        ru += G.nodes[n]['ia']

    return ru


'''
Calculate te semantic distance (sd) between two drugs by summing the mis-information and remaining uncertainty values.
'''
def semantic_distance(drug1, drug2):
    sd = compute_mi(drug1, drug2) + compute_ru(drug1, drug2)
    return sd


'''
Calculate the semantic distance between every pairwise combination of drugs, no repeats.
'''
def run_comparisons(comparisons):

    # make new array and add new column for mi value
    zeros = np.zeros((np.shape(comparisons)[0],1))
    comparisons = np.column_stack((comparisons, zeros))


    for i,c in enumerate(comparisons):
        distance = semantic_distance(c[0], c[1])
        comparisons[i][2] = distance

        # write to table
        # res.loc[c[0], c[1]] = distance
        # res.loc[c[1], c[0]] = distance

    return comparisons


'''
Initializer for multiprocessing to generate global variables to use in each proces.
'''
def initializer():
    global chembl
    global mesh_headings
    global mesh_numbers
    # global res

    # load in ChEMBL
    chembl = pd.read_csv('data/chembl_indications.tsv', sep='\t')

    # load in MeSH pickle files
    f = open("data/mesh_headings.pkl", "rb")
    mesh_headings = pickle.load(f)
    f = open("data/mesh_numbers.pkl", "rb")
    mesh_numbers= pickle.load(f)

    # res = pd.DataFrame(columns=drug_list)
    # for drug in drug_list:
    #     res = res.append(pd.Series(name=drug))


'''
Main function for each process. Builds a tree for the drugs included in the comparisons, computes the information accretion (ia)
value for each node, and then computes all comparisons.
'''
def main(comparisons):
    global drug_node_dict
    global G

    drug_node_dict = {} # initialize dict key: drug value: nodes in its network
    drugs = set([j for i in comparisons for j in i]) # extract individual drugs from comparisons to build graph
    
    # main process
    G = make_graph(drugs) # build graph
    G = compute_ia(G) # compute and add ia values
    distances = run_comparisons(comparisons) # compute semantic distances between all drugs

    return distances


'''
Write results to table. this needs to be faster!
'''
def write_results(distances):
    
    res = pd.DataFrame(columns=all_drugs)
    for drug in all_drugs:
        res = res.append(pd.Series(name=drug))

    for d in distances:
        # get values
        drug1 = d[0]
        drug2 = d[1]
        distance = d[2]

        # write to table
        res.loc[drug1, drug2] = distance
        res.loc[drug2, drug1] = distance

    res.to_csv('drug_distances.csv', sep=',', index=True, na_rep=0, index_label='Drug')


'''
Main.
'''
if __name__ == "__main__":

    print('setting things up...')

    # get arguments
    args = parseArgs() # parse arguments
    n = args.num_cpus # number of processes/cpus to use
    rows = args.rows # number of rows of indications to use in ChEMBL (for testing on less data)

    print('loading data...')

    # load in ChEMBL
    chembl = pd.read_csv('data/chembl_indications.tsv', sep='\t') # Global.
    all_drugs = list(set(chembl['pref_name'][:rows])) # chembl[chembl['therapeutic_flag'] == 1]['pref_name']

    # generate all drug comparisons
    all_comparisons = list(combinations(all_drugs,2))

    # split drugs into nearly equal parts for every process
    l = len(all_comparisons) // n
    if len(all_comparisons) % n != 0:
        l += 1
    lists = [all_comparisons[i:i + l] for i in range(0, len(all_comparisons), l)]


    print('running processes...')

    # compute distances across n threads
    with Pool(n, initializer, ()) as p:
        distances = p.map(main, lists)

    # write results to csv
    # res.to_csv('drug_distances.csv', sep=',', index=True, na_rep=0, index_label='Drug')

    distances = [j for i in distances for j in i]

    print('writing results...')
    write_results(distances)