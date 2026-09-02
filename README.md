# KIMMDY

This is a development branch for KIMMDY with the new and updated barrier-hat-models based on MACE.

It includes further many updates. I recommend investigating the git history if you want to assess them.

An exemplary .yaml including newer functions can be seen below.

[![](https://github.com/graeter-group/kimmdy/actions/workflows/tests.yml/badge.svg?branch=release-please--branches--main)](https://github.com/graeter-group/kimmdy/actions/workflows/tests.yml/?branch=release-please--branches--main)

Reactive MD pipeline for GROMACS using Kinetic Monte Carlo / Molecular Dynamics (KIMMDY)

Please create a fork or a personal branch if you want to do any changes, as this branch version is still being developed.

Further please keep this branch and clones private for now, as it may contain structures with a pending patent. This includes the use of coding agents that train on used conversations, to be sure.

## Installation

**Note**: KIMMDY requires [GROMACS](https://www.gromacs.org/) to be installed, or be available as a loadable plugin on your cluster.
Some reactions need a GROMACS version patched with [PLUMED](https://www.plumed.org/).
The gromacs version name should then contain `MODIFIED` or `plumed`.

Further this new branch requires the new HAT-plugin found in a different repository (https://github.com/graeter-group/KIMMDY-MACE-HAT-plugin)

### Installation on MPCDF login nodes (optimized for Otter)

Prerequisites:

- SSH access to GitHub
- Access to both private repositories

#### 1. Create The uv Environment

```bash
module purge
module load anaconda/3/2023.03

eval "$(conda shell.bash hook)"

conda create \
    --channel conda-forge \
    python=3.11 \
    uv=0.9.10 \
    --yes

conda activate uv
```

This environment only provides `uv`. The actual KIMMDY environment will be the repository’s `.venv`.

#### 2. Clone KIMMDY

Do not initialize submodules immediately because the committed `.gitmodules` still points to the public HAT repository.

```bash
mkdir -p "$HOME/software/KIMMDY"
cd "$HOME/software/KIMMDY"

git clone \
    --branch MACE_HAT \
    git@github.com:graeter-group/KIMMDY-MACE-HAT.git \
    kimmdy

cd kimmdy
```

#### 3. Configure The Private HAT Submodule

```bash
git config \
    submodule.plugins/kimmdy-hat.url \
    git@github.com:graeter-group/KIMMDY-MACE-HAT-plugin.git

git submodule update --init --recursive
```

Verify that everything was initialized:

```bash
git submodule status --recursive
```

There should be no submodule lines beginning with `-`.

#### 4. Install KIMMDY

```bash
uv sync --extra all
```

This installs:

- KIMMDY
- all reaction plugins
- the private editable HAT plugin
- the local modified MACE implementation
- Grappa
- development and testing dependencies
- no legacy TensorFlow/PaiNN stack

#### 5. Install GPU-Enabled Torch


The login node does not expose the compute-node GPU, so select the tested backend explicitly:

Activate the kimmdy environment:

```bash
source .venv/bin/activate
```

Install torch-cuda with the fixed version required by GRAPPA:

```bash
uv pip install --reinstall "torch==2.0.0" --torch-backend=cu118
```

Verify installed packages:

```bash
uv pip show \
    kimmdy \
    kimmdy-hat \
    mace-torch \
    kimmdy-grappa \
    grappa-ff \
    torch
```

#### 6. Use In Slurm Jobs

The relevant initialization block is:

```bash
module purge
module load anaconda/3/2023.03
module load gromacs-plumed/2024.6

eval "$(conda shell.bash hook)"
conda activate "$HOME/conda-envs/uv"

source "$HOME/software/KIMMDY/kimmdy/.venv/bin/activate"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
```

A KIMMDY configuration for this cluster can then use:

```yaml
gromacs_alias: gmx_mpi
mdrun_prefix: "srun --ntasks=1"
```

After running `uv sync` again later, reinstall the CUDA Torch wheel:

```bash
uv pip install --reinstall "torch==2.0.0" --torch-backend=cu118
```

#### 7. Development

If you plan to modify the software, create a separate development branch in the main KIMMDY repository and in every submodule you intend to modify, such as the `kimmdy-hat` repository.

Changes made inside a submodule must be committed in that submodule first. The resulting submodule commit must then be recorded and committed in the parent KIMMDY repository.

## Documentation

The documentation is available [here](https://graeter-group.github.io/kimmdy/).

### Getting started

Head over to the [getting started](https://graeter-group.github.io/kimmdy/guide/tutorials/getting-started.html) tutorial.

## Development

The above instructions install all repositories in editable mode, allowing changes.

Run kimmdy either using `uv run kimmdy`, or activate the virtual environment uv has created:

```bash
source .venv/bin/activate
kimmdy --help
```

Conventions:

- code style: black
- docstrings: numpy
- [Conventional commit](https://www.conventionalcommits.org/en/v1.0.0/) messages when possible for pretty release notes.

For more information, head to
[the documentation](https://graeter-group.github.io/kimmdy/guide/how-to/contribute.html)


## Exemplary yamls and submit script:

### Charmm-GUI-style .itp topolgy

<details>
<summary>Example: CHARMM lipid HAT with MACE and Eyring rates</summary>

```yaml
# Execute simulations rather than only printing the planned tasks.
dryrun: false

# Each replica has its own directory, so this name can be identical in every rep.
name: "ERY_DOPA_MACE"

# Allows enough tasks for setup plus ten reaction cycles.
max_tasks: 100

# Do not parameterize the complete membrane during setup.
# Grappa is invoked after reactions for affected bonded neighborhoods.
parameterize_at_setup: false

# MPI-enabled, CUDA-capable GROMACS executable.
gromacs_alias: "gmx_mpi"

# Launch each mdrun as one Slurm MPI task.
mdrun_prefix: "srun --ntasks=1"

# Match -ntomp to SLURM_CPUS_PER_TASK in the accompanying job script.
gmx_mdrun_flags: >-
  -maxh 23.5
  -dlb auto
  -pin on
  -ntomp 32
  -pme gpu
  -bonded gpu
  -nb gpu

# CHARMM force-field directory relative to repXX.
ff: "../../Topology/toppar.ff"

# Starting topology shared by all replicas.
top: "../../Topology/topol.top"

# Equilibrated starting coordinates for this replica.
gro: "step6.6_equilibration.gro"

# Index generated for the radicalized system.
ndx: "../index_rad.ndx"

# Oxygen atoms of the three initial hydroxyl radicals.
# These are one-based GROMACS atom IDs.
radicals: "22905 22907 22909"

topology:
  reactive:
    # Keep solvent and ions outside the merged Reactive moleculetype.
    # OHR remains reactive because it participates in HAT.
    exclude: "TIP3 POT CLA"

  # Use CHARMM-style molecule ITPs as product templates.
  product_template_sources:
    - itp

  # Discover molecule templates in the configured force-field directory.
  product_template_itp_glob: "*.itp"

  # Permit hydrogenated OHR to match the OHRH product template,
  # including the configured ORA/HRA to OT/HT atom-type transition.
  product_template_transitions:
    OHR: OHRH

mds:
  sample:
    # Sampling performed before each reaction search.
    mdp: "../../mdp_templates/md_sample.mdp"

  equilibrate:
    # Conventional equilibration after each accepted reaction.
    mdp: "../../mdp_templates/md_eq.mdp"

  slow_growth:
    # Lambda relaxation while applying the hydrogen transfer.
    mdp: "../../mdp_templates/md_growth.mdp"

changer:
  topology:
    # Reparameterize affected bonded neighborhoods after topology changes.
    parameterization: grappa
    parameterization_kwargs:
      # Explicitly use the Grappa model trained for radical chemistry.
      grappa_tag: "grappa-1.4.1-radical"

      # Interpret the input nonbonded parameters as CHARMM-style parameters.
      charge_model: "charmm"

  coordinates:
    # MD configuration used for lambda coordinate changes.
    md: slow_growth

    # Interpolate the complete set of changed topology parameters.
    slow_growth: full

reactions:
  hat_reaction:
    # Select transition-state-theory rates using the Eyring equation.
    rate_equation: eyring

    eyring_equation:
      # Must be consistent with the MD temperature.
      temperature: 310

      # Probability of successful barrier crossing.
      # A value of 1.0 corresponds to conventional transition-state theory.
      transmission_coefficient: 1.0

    # Use the MACE barrier-prediction model.
    model: mace

    # Generic capping is required for MACE input structures.
    cap: true

    # Evaluate all ten available ensemble members.
    ensemble_size: 10

    # Use the GPU visible in the Slurm allocation.
    mace_device: cuda

    # Number of structures evaluated in one MACE batch.
    mace_batch_size: 64

    # Select model files from the bundled model directory.
    mace_model_pattern: "*.model"

    # Preserve the dtype used to train the models.
    mace_default_dtype: "float64"

    # Unit produced by the MACE barrier models.
    mace_barrier_unit: eV

    # Maximum donor-H search distance in angstrom.
    h_cutoff: 3

    # Inspect every stored trajectory frame.
    polling_rate: 1

    # Keep at most 100 geometrically best frames per reaction identity.
    n_unique: 100

    # Use the sampled XTC trajectory for candidate generation.
    trajectory_format: xtc

    # Transfer hydrogen through the configured slow-growth procedure.
    change_coords: lambda

    # Use rejection-free kinetic Monte Carlo.
    kmc: rfkmc

    # Retain generated prediction-system data for subsequent analysis.
    keep_structures: true

    # Write capped structures as PDB files for visual inspection.
    write_pdbs: true

sequence:
  # Perform exactly ten sample -> reaction -> equilibration cycles.
  - mult: 10
    tasks:
      - sample
      - reactions
      - equilibrate
```

</details>

### Classic .rtp style topology

<details>
<summary>Example: RTP protein topology with MACE and Eyring rates</summary>

```yaml
dryrun: false
name: "PROTEIN_HAT_MACE"
max_tasks: 100

# Parameterize only reaction-affected regions.
parameterize_at_setup: false

# HPC-specific GROMACS execution.
gromacs_alias: "gmx_mpi"
mdrun_prefix: "srun --ntasks=1"
gmx_mdrun_flags: >-
  -maxh 23.5
  -dlb auto
  -pin on
  -ntomp 32
  -pme gpu
  -bonded gpu
  -nb gpu

# Standard GROMACS force-field directory containing aminoacids.rtp.
ff: "amber99sb-star-ildnp.ff"

top: "topol.top"
gro: "equilibrated.gro"
ndx: "index.ndx"

# One-based GROMACS IDs of the initial radical atoms.
radicals: "1234"

topology:
  reactive:
    # Adjust these names to match the solvent and ions in the topology.
    exclude: "SOL NA CL"

  # Use the conventional residue topology database.
  product_template_sources:
    - rtp

mds:
  sample:
    mdp: "md_sample.mdp"

  equilibrate:
    mdp: "md_eq.mdp"

  slow_growth:
    mdp: "md_growth.mdp"

changer:
  topology:
    parameterization: grappa
    parameterization_kwargs:
      grappa_tag: "grappa-1.4.1-radical"
      charge_model: "amber99"

  coordinates:
    md: slow_growth
    slow_growth: full

reactions:
  hat_reaction:
    rate_equation: eyring

    eyring_equation:
      temperature: 310
      transmission_coefficient: 1.0

    model: mace
    cap: true
    ensemble_size: 10

    mace_device: cuda
    mace_batch_size: 64
    mace_model_pattern: "*.model"
    mace_default_dtype: "float64"
    mace_barrier_unit: eV

    h_cutoff: 3
    polling_rate: 1
    n_unique: 100
    trajectory_format: xtc
    change_coords: lambda
    kmc: rfkmc

    keep_structures: true
    write_pdbs: true

sequence:
  - mult: 10
    tasks:
      - sample
      - reactions
      - equilibrate
```

> **RTP limitation:** RTP matching restores the transferred hydrogen identity, atom type, and relevant charges. Complete classical bonded-parameter restoration is currently implemented only for matching single-residue ITP templates; RTP products retain Grappa bonded parameters.

</details>

### Submit script (Otter)

<details>
<summary>Example: SLURM submit script for kimmdy run on Otter using 1 GPU</summary>

```bash
#!/bin/bash -l
#SBATCH -o kimmdy_ex.out.%j
#SBATCH -e kimmdy_ex.err.%j
#SBATCH -D ./
#SBATCH -J kimmdy_ex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=s.otter
#SBATCH --mem=30G
#SBATCH --gres=gpu:l40s:1
#SBATCH --time=24:00:00
#SBATCH --mail-type=none

set -euo pipefail

module purge
module load anaconda/3/2023.03
module load gromacs-plumed/2024.6

eval "$(conda shell.bash hook)"
conda activate uv
source ~/software/KIMMDY/kimmdy/.venv/bin/activate

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTHONUNBUFFERED=1

echo "Job: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "OpenMP threads: $OMP_NUM_THREADS"

echo "KIMMDY run completed successfully."
```

<details>