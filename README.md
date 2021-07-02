# DrugSimilarity
A python script for comparing drugs based on the semantic distance between the [MeSH](https://www.nlm.nih.gov/mesh/meshhome.html) headings of their [ChEMBL](https://www.ebi.ac.uk/chembl/) indications.


### Usage
```
python3 computeDistances.py -n 10
```


### Parameter Reference
```
usage: computeDistances.py [-h] [-n NUM_CPUS] [-r ROWS]

Compute the semantic distance between drugs in ChEMBL using the MeSH headings of their indications.

optional arguments:
  -h, --help            show this help message and exit
  -n NUM_CPUS, --num-cpus NUM_CPUS
                        The number of cpus to use. Default is 4.
  -r ROWS, --rows ROWS  The number of rows of ChEMBL indications to use. Default is all.
```


### Details
computeDistances.py utilizes [NetworkX](https://networkx.org/), multiprocessing, and curated files from [ChEMBL](https://www.ebi.ac.uk/chembl/) and [MeSH](https://www.nlm.nih.gov/mesh/meshhome.html) to quickly calculate the semantic distance between all drugs in [ChEMBL](https://www.ebi.ac.uk/chembl/).


### Future Directions
Data generated from this script will be used to add a new similarity metric to the [Small Molecule Suite](https://labsyspharm.shinyapps.io/smallmoleculesuite/). We also anticipate
clustering drugs based on their semantic distance to verify accurate similarity computations.