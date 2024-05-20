jobscript = """
#!/bin/env bash
#SBATCH --job-name={config.out.name}
#SBATCH --output={config.out.name}-job.log
#SBATCH --error={config.out.name}-job.log
#SBATCH --time={config.max_hours}:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mincpus=20
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --gpus 1
#SBATCH --mail-type=ALL
# #SBATCH -p <your-partition>.p
# #SBATCH --mail-user=<your-email>


# Setup up your environment here
# modules.sh might load lmod modules, set environment variables, etc.
source ./_modules.sh

CYCLE={config.max_hours}
CYCLE_buffered=$(echo "scale=2; $CYCLE - 0.08" | bc)


START=$(date +"%s")

timeout ${{CYCLE_buffered}}h kimmdy -i {config.input_file}

END=$(date +"%s")

LEN=$((END-START))
HOURS=$((LEN/3600))

echo "$LEN seconds ran"
echo "$HOURS full hours ran"

if [ $HOURS -lt $CYCLE ]; then
  echo "last cycle was just $HOURS h long, KIMMDY is done."
  exit 3
else
  echo "jobscript resubmitting"
  sed -i.bak "s/\(run_directory:\s*\).*/\\1'{config.out.name}'/" {config.input_file}
  kimmdy --generate-jobscript
  sbatch ./jobscript.sh
  exit 2
fi
"""
