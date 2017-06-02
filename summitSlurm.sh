#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by sbatch
# as arguments.  For more details about usage of these arguments see "man sbatch"

# Set a walltime for the job. The time format is HH:MM:SS

# Run for 24 hours:
#SBATCH --time=23:59:59

# Select one nodes and processors

#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 1

# Set the output file name to [jobid].out (or leave as default of slurm-[jobid].out)

#SBATCH -o hr.out
#SBATCH -e hr.err

# Allocation
#SBATCH -A ucb-summit-smr

# Select the janus-short QOS (comperable to a queue)
#SBATCH --qos=condo

# partition
#SBATCH --partition=shas

# Load any modules you need here
module load slurm
module load matlab/matlab-2013b

# Execute the program.
echo "Time is `date`"
echo "Submit dir ${SLURM_SUBMIT_DIR}"
echo "Running ${SLURM_NNODES} nodes. ${SLURM_NTASKS_PER_NODE} tasks per node. ${SLURM_CPUS_PER_TASK} processors per task"
echo "In dir `pwd`"
touch jobRunning.txt
# Run matlab program
matlab -nodesktop -nosplash \
  -r  "try, runHardRod, catch, exit(1), end, exit(0);" \
  2>&1 | tee ${SLURM_JOB_NAME}.out
echo "Finished. Matlab exit code: $?" 
mv jobRunning.txt jobFinish.txt
echo "Time is `date`"
