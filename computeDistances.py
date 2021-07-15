from os import sep
import numpy as np
import csv
import argparse
import pickle
import networkx as nx
import math
import pandas as pd
from multiprocessing import Pool


'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Compute the semantic distance between drugs in ChEMBL using the MeSH headings of their indications.')
    parser.add_argument('-n', '--num-cpus', help="The number of cpus to use. Default is 4.", default=4, type=int, required=True)
    parser.add_argument('-g', '--graph', help='A pickle file of a networkx graph of ChEMBL.', type=str, required=True)
    parser.add_argument('-d', '--drugs', help='A file containing a list of drugs for computation', type=str, required=True)
    parser.add_argument('-a', '--all-drugs', help='A file containing a list of all drugs in ChEMBL.', type=str, required=True)
    parser.add_argument('-p', '--drug-node-dict', help='A pickle file containing drug_node_dict.', type=str, required=True)
    parser.add_argument('-i', '--id', help='A unique id to use for this scripts output file.', type=str, required=True)
    parser.add_argument('-o', '--output', help='Output directory.', type=str, required=True)
    args = parser.parse_args()
    return args


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

    # run comparisons
    for i,c in enumerate(comparisons):
        distance = semantic_distance(c[0], c[1]) # compute distance
        comparisons[i][2] = distance # add to data

    return comparisons


'''
NEW
'''
def run_comparisons_new(drugs, all_drugs):

    results = []
    for drug1 in drugs:
        distances = [drug1]
        for drug2 in all_drugs:
            distance = semantic_distance(drug1, drug2) # compute distance
        distances.append(distance)
        results.append(distances)

    return results
    

'''
Initializer for multiprocessing to generate global variables to use in each proces.
'''
def initializer():
    global all_drugs
    global drug_node_dict
    global G

    # load dict -> key: drug, value: nodes in its network
    f = open(args.drug_node_dict, 'rb')
    drug_node_dict = pickle.load(f)

    # load in all drugs
    f = open(args.all_drugs, 'r')
    all_drugs = [line.rstrip() for line in f]

    # read in graph
    G = nx.read_gpickle(args.graph)


'''
Main function for each process. Computes all comparisons.
'''
def main(drugs):
    
    # generate comparisons (list of tuples?)
    '''
    comparisons = []
    for drug1 in drugs:
        for drug2 in all_drugs:
            comparisons.append((drug1, drug2))

    distances = run_comparisons(comparisons) # compute semantic distances between all drugs
    '''

    results = run_comparisons_new(drugs, all_drugs)

    return results


'''
Write results to table.
'''
def write_results(distances):
    
    # initialize df
    res = pd.DataFrame(columns=all_drugs)
    for drug in drugs:
        res = res.append(pd.Series(name=drug))

    for d in distances:
        # get values
        drug1 = d[0]
        drug2 = d[1]
        distance = d[2]

        # write to table
        res.loc[drug1, drug2] = distance
        res.loc[drug2, drug1] = distance

    # write to tsv. exclude header if not first child job so they can be easily merged later
    if args.id != 'aa':
        res.to_csv(f'{args.output}/drug_distances_{args.id}.tsv', sep='\t', index=True, na_rep=0, index_label='Drug', header=False)
    else:
        res.to_csv(f'{args.output}/drug_distances_{args.id}.tsv', sep='\t', index=True, na_rep=0, index_label='Drug')
    

'''
NEW
'''
def write_results_new(distances):
    
    with open(f'{args.output}/drug_distances_{args.id}.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(distances)

    
'''
Main.
'''
if __name__ == "__main__":

    global args # hopefully this allows the subprocesses to access them

    # get arguments
    args = parseArgs() # parse arguments
    n = args.num_cpus # number of processes/cpus to use

    print('loading data...')

    # load in small drug list
    f = open(args.drugs, 'r')
    drugs = [line.rstrip() for line in f]

    # load in all drugs
    f = open(args.all_drugs, 'r')
    all_drugs = [line.rstrip() for line in f]

    # divide drugs for processes
    num_drugs = math.ceil(len(drugs) / n)
    lists = [drugs[i:i + num_drugs] for i in range(0, len(drugs), num_drugs)]
  
    print('running processes...')

    # compute distances across n threads
    with Pool(n, initializer, ()) as p:
        distances = p.map(main, lists)

    # merge results
    distances = [j for i in distances for j in i]

    print('writing results...')
    write_results_new(distances)