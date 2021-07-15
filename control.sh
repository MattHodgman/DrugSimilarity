num_jobs=$1 # the number of jobs to submit

# make drug list
cut -d$'\t' -f2 data/chembl_indications.tsv | sort | uniq > drugs

# split drug list into a file for each job
split -n l/${num_jobs} drugs 'data/drugs_'

# run a job on each drug list
ls data/drugs_* | xargs -I {} sbatch child.sh {}