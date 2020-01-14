#!/bin/bash

# submit the job by $bash file_name

#n_ens=2  # number of ensemble members
max_it=3 # number of iteration

# first job "truth" - no dependencies
id_truth=$(sbatch --parsable -A esm run_1m_truth.sh)

# if job truth succeeds setup the eki class
id_init_ens=$(sbatch --parsable --dependency=afterok:$id_truth -A esm run_1m_eki_init.sh)

for i in $(seq 1 $max_it); do

    if [ $i -eq 1 ]; then
        # array of jobs depending on a single job eki class setup
        id_ens=$(sbatch --parsable --dependency=afterok:$id_init_ens -A esm --array=1-2 run_1m_ensemble.sh)
    else
        id_ens=$(sbatch --parsable --dependency=afterok:$id_after_ens -A esm --array=1-2 run_1m_ensemble.sh)
    fi

    # last job depending on the success of the array of jobs
    id_after_ens=$(sbatch --parsable --dependency=afterok:$id_ens -A esm run_1m_after_ensemble.sh)

done
