jobscript = """
#!/bin/env bash
#SBATCH --job-name={config.out.name}
#SBATCH --output={config.out.name}-job.log
#SBATCH --error={config.out.name}-job.log
#SBATCH --time={config.max_hours}:00:00
#SBATCH --nodes={config.slurm.N}
#SBATCH --ntasks-per-node={config.slurm.ntasks_per_node}
#SBATCH --mincpus={config.slurm.mincpus}
#SBATCH --exclusive
#SBATCH --cpus-per-task={config.slurm.cpus_per_task}
#SBATCH --gpus={config.slurm.gpus}
#SBATCH --mail-type=ALL
# # uncomment these to use:
# #SBATCH --partition <your-partition>.p
# #SBATCH --mail-user=<your-email>


# Setup up your environment here
# modules.sh might load lmod modules, set environment variables, etc.
if [ -f ./_modules.sh ]; then
    source ./_modules.sh
fi

CYCLE={config.max_hours}
CYCLE_buffered=$(echo "scale=2; $CYCLE - 0.08" | bc)


START=$(date +"%s")

timeout ${{CYCLE_buffered}}h kimmdy -i {config.input_file} --restart

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
  echo "jobscript resubmitting"
  {config.slurm.runcmd} ./jobscript.sh
  exit 2
fi
"""
