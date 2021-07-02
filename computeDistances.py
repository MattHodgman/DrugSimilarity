import numpy as np
import pickle
import networkx as nx
import math
from itertools import combinations, permutations
import pandas as pd
import sys
from multiprocessing import Pool
import operator as op
from functools import reduce

# n choose 2
def ncr(n):
    r = 2
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

# get the parent number (move UP the tree one level/generation)
def up(n):
    sep = '.'
    n = n.split(sep)
    if len(n) > 1:
        n.pop()
        n = sep.join(n)
    else:
        return n[0][0]
    return n

# get indications headings
def get_headings(drug):
    if 'CHEMBL' in drug:
        indications = chembl[chembl['chembl_id'] == drug.upper()]
    else:
        indications = chembl[chembl['pref_name'] == drug.upper()]
    headings = sorted(list(set(indications.mesh_heading)))
    return headings

# make a graph where the nodes are MeSH headings and the directed edges form the heirarchy.
# each node has as an attribute the list of drugs that include that node
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


# compute information accretion for each node in a graph
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



# calculate mis-information between two drug graphs
def compute_mi(drug1, drug2):
    nodes = np.setdiff1d(list(drug_node_dict[drug1]), list(drug_node_dict[drug2]))

    mi = 0
    for n in nodes:
        mi += G.nodes[n]['ia']

    return mi


# compute remaining uncertainty between two drug graphs
def compute_ru(drug1, drug2):
    nodes = np.setdiff1d(list(drug_node_dict[drug2]), list(drug_node_dict[drug1]))

    ru = 0
    for n in nodes:
        ru += G.nodes[n]['ia']

    return ru


# calculate te semantic distance between two drugs by summing the mis-information and remaining uncertainty values
def semantic_distance(drug1, drug2):
    sd = compute_mi(drug1, drug2) + compute_ru(drug1, drug2)
    return sd


# calculate the semantic distance between every pairwise combination of drugs, no repeats
def run_comparisons(comparisons):

    # make all comparisons
    comparisons = np.array(comparisons)

    # make new array and add new column for mi value
    # zeros = np.zeros((np.shape(comparisons)[0],1))
    # comparisons = np.column_stack((comparisons, zeros))

    for i,c in enumerate(comparisons):
        distance = semantic_distance(c[0], c[1])

        # write to table
        res.loc[c[0], c[1]] = distance
        res.loc[c[1], c[0]] = distance


# write results to table
# this needs to be faster!
def write_results(distances):
    
    res = pd.DataFrame(columns=drug_list)
    for drug in drug_list:
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

def initializer():
    global chembl
    global mesh_headings
    global mesh_numbers
    global res

    # load in ChEMBL
    chembl = pd.read_csv('data/chembl_indications.tsv', sep='\t')
    drug_list = list(set(chembl['pref_name'][:int(sys.argv[2])]))

    # load in MeSH pickle files
    f = open("data/mesh_headings.pkl", "rb")
    mesh_headings = pickle.load(f)
    f = open("data/mesh_numbers.pkl", "rb")
    mesh_numbers= pickle.load(f)

    res = pd.DataFrame(columns=drug_list)
    for drug in drug_list:
        res = res.append(pd.Series(name=drug))


def main(comparisons):
    global drug_node_dict
    global drugs_array
    global G

    drugs = set([j for i in comparisons for j in i])

    drug_node_dict = {}
    drugs_array = np.array(list(drugs))

    G = make_graph(drugs)
    G = compute_ia(G)
    run_comparisons(comparisons)

# main
if __name__ == "__main__":

    # num processes
    n = int(sys.argv[1])

    # num rows
    r = int(sys.argv[2])

    # load in ChEMBL
    chembl = pd.read_csv('data/chembl_indications.tsv', sep='\t')
    drug_list = list(set(chembl['pref_name'][:r])) # chembl[chembl['therapeutic_flag'] == 1]['pref_name']

    all_comparisons = list(combinations(drug_list,2))


    # split drugs into nearly equal parts for every process
    l = len(all_comparisons) // n
    if len(all_comparisons) % n != 0:
        l += 1
    lists = [all_comparisons[i:i + l] for i in range(0, len(all_comparisons), l)]

    # compute distances across n threads
    with Pool(n, initializer, ()) as p:
        p.map(main, lists)


    res.to_csv('drug_distances.csv', sep=',', index=True, na_rep=0, index_label='Drug')

    # distances = [j for i in distances for j in i]

    # print('writing results')
    # write_results(distances)