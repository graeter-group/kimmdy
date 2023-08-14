
#!/bin/env bash
#SBATCH --job-name=alanine_hat
#SBATCH --output=kimmdy-job.log
#SBATCH --error=kimmdy-job.log
# #SBATCH -p <your-partition>.p
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --mincpus=40
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --gpus 1
#SBATCH --mail-type=ALL
# #SBATCH --mail-user=<your-email
echo "Job Name:"
echo "$SLURM_JOB_NAME"
JOB_NAME = "$SLURM_JOB_NAME"

# Setup up your environment here
# modules.sh might load lmod modules, set environment variables, etc.
source ./_modules.sh

JOB=alanine_hat
CYCLE=20
SUBMIT="jobscript.sh $JOB"

START=$(date +"%s")

END=$(date +"%s")
LEN=$((END-START))
HOURS=$((LEN/3600))
echo "$LEN seconds ran"
echo "$HOURS full hours ran"
let "CYCLE--"
if [ $HOURS -lt $CYCLE ]; then
  echo "last cycle was just $HOURS h long and therefore finito"
  rm $targetRunning
  exit 3
else
  echo "cycle resubmitting"
  sbatch -J $SLURM_JOB_NAME $SUBMIT
  exit 2
fi
