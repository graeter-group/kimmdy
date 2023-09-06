jobscript = """
#!/bin/env bash
#SBATCH --job-name={config.out.name}
#SBATCH --output=kimmdy-job.log
#SBATCH --error=kimmdy-job.log
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --mincpus=40
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --gpus 1
#SBATCH --mail-type=ALL
# #SBATCH -p <your-partition>.p
# #SBATCH --mail-user=<your-email


# Setup up your environment here
# modules.sh might load lmod modules, set environment variables, etc.
# source ./_modules.sh

CYCLE=24

START=$(date +"%s")

kimmdy -p {config.out}/kimmdy.cpt

END=$(date +"%s")

LEN=$((END-START))
HOURS=$((LEN/3600))

echo "$LEN seconds ran"
echo "$HOURS full hours ran"

let "CYCLE--"
if [ $HOURS -lt $CYCLE ]; then
  echo "last cycle was just $HOURS h long, KIMMDY is done."
  exit 3
else
  echo "cycle resubmitting"
  # sbatch ./jobscript.sh
  ./jobscript.sh
  exit 2
fi
"""
