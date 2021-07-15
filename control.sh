num_jobs=$1 # the number of jobs to submit

# make drug list
tail -n+2 data/chembl_indications.tsv | cut -d$'\t' -f2 | sort | uniq > data/drugs

# split drug list into a file for each job
split -n l/${num_jobs} data/drugs 'data/drugs_'

# make output dir for tsvs
mkdir output

# run a job on each drug list
ls data/drugs_* | xargs -I {} sbatch child.sh {} output