# DrugSimilarity
A python script for comparing drugs based on the semantic distance between the [MeSH](https://www.nlm.nih.gov/mesh/meshhome.html) headings of their [ChEMBL](https://www.ebi.ac.uk/chembl/) indications.

### Details
computeDistances.py utilizes [NetworkX](https://networkx.org/), multiprocessing, and curated files from [ChEMBL](https://www.ebi.ac.uk/chembl/) and [MeSH](https://www.nlm.nih.gov/mesh/meshhome.html) to quickly calculate the semantic distance between all drugs in [ChEMBL](https://www.ebi.ac.uk/chembl/).

### Future Directions
Data generated from this script will be used to add a new similarity metric to the [Small Molecule Suite](https://labsyspharm.shinyapps.io/smallmoleculesuite/). We also anticipate
clustering drugs based on their semantic distance to verify accurate similarity computations.