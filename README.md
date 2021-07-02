# DrugSimilarity
A python script for comparing drugs based on the semantic distance between the MeSH headings of their ChEMBL indications.

### Details
computeDistances.py utilizes networkx, multiprocessing, and curated files from ChEMBL and MeSH to quickly calculate the semantic distance between all drugs in ChEMBL.

### Future Directions
Data generated from this script will be used to add a new similarity metric to the [Small Molecule Suite](https://labsyspharm.shinyapps.io/smallmoleculesuite/). We also anticipate
clustering drugs based on their semantic distance to verify accurate similarity computations.
