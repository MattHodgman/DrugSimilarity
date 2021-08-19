'''
for all pairwise comparisons
    get the list of indications of both drugs
    compute similarity (intersection / smaller of two sets)

its that simple lol
'''

from itertools import combinations
import numpy as np
import pandas as pd
import argparse


'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Compute the semantic distance between drugs in ChEMBL using the MeSH headings of their indications.')
    parser.add_argument('-n', '--num-cpus', help="The number of cpus to use. Default is 4.", default=4, type=int, required=False)
    parser.add_argument('-c', '--chembl', help='A TSV file containing ChEMBL drug indication information. Must contain the drug name under the column \'pref_name\' and the indication name under \'mesh_heading\'.', type=str, required=True)
    parser.add_argument('-o', '--output', help='Output directory.', default='.', type=str, required=False)
    args = parser.parse_args()
    return args


'''
Get a drug's indications
'''
def getIndications(drug):
    indications = chembl[chembl['pref_name'] == drug]['mesh_heading'].tolist()
    return indications


'''
Compute the similarity of two drugs based on their indications.
Similarity = intersection / smaller set
'''
def similarity(drug1, drug2):
    # get the indications for each drug
    indications1 = getIndications(drug1)
    indications2 = getIndications(drug2)

    intersection = len(list(set(indications1) & set(indications2))) # get the size of the intersection of both sets

    # get the size of the smaller set
    if len(indications1) < len(indications2):
        smaller = len(indications1)
    else:
        smaller = len(indications2)

    # compute similarity components
    similarity = intersection / smaller

    return similarity


'''
The big O squared way of doing this
'''
def bruteForce():

    drugs = sorted(chembl['pref_name'].unique().tolist())

    # compute and write
    with open(f'{out}/similarities.tsv', 'w') as f:
        
        # write header
        for d in drugs:
            f.write(d + '\t')
        
        f.write('\n') # newline

        # compute and write similarities
        for drug1 in drugs:
            f.write(drug1) # write drug name
            for drug2 in drugs:
                s = str(similarity(drug1, drug2)) # compute similarity
                f.write('\t' + s) # write similarity
            f.write('\n') # newline


'''
Vectorization, fast but needs >350GB mem lol
'''
def vectorIt():
    # get data
    df = pd.crosstab(chembl.pref_name,[chembl.mesh_heading]).T
    rawdata = df.to_numpy()
    ndrugs = rawdata.shape[1]
    r,c = np.tril_indices(ndrugs,-1)

    # Use those indicees to slice out respective columns 
    p1 = rawdata[:,c]
    p2 = rawdata[:,r]

    # Perform n11 and n00 vectorized computations across all indexed columns
    n11v = ((p1==1) & (p2==1)).sum(0) # the number of indications that both drugs have
    nmin = np.minimum((p1==1).sum(0), (p2==1).sum(0)) # get the size of the smaller set of indications

    # Finally, setup output array and set final division computations
    out = np.eye(ndrugs)
    out[c,r] = n11v / (nmin)

    df2 = pd.DataFrame(out, columns = df.columns, index = df.columns)
    df2.to_csv(f'{out}/similarities.tsv', sep='\t')


'''
Main.
'''
if __name__ == "__main__":

    # get arguments
    args = parseArgs() # parse arguments
    n = args.num_cpus # number of processes/cpus to use
    out = args.output # get the output dir path

    # load ChEMBL
    chembl = pd.read_csv(args.chembl, sep='\t')