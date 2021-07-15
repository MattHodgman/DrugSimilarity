num_jobs=$1 # the number of jobs to submit
chembl=$2 # chembl_indications.tsv

# make drug list
cut -d$'\t' -f2 ${chembl} | sort | uniq > drugs

# split drug list into a file for each job
split -n l/${num_jobs} drugs 'drugs_'

# run a job on each drug list
ls drugs_* | xargs -I {} sbatch job.sh -d {} -a drugs -g chembl.gpickle -p drug_node_dict.pickle