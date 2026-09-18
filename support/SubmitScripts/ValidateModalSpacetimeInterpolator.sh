#!/bin/bash
# Distributed under the MIT License.
# See LICENSE.txt for details.

# Example Slurm job array for validating one node-written volume file per job.
# Submit with, for example:
#
#   SPECTRE_BUILD_DIR=/u/$USER/build-Release \
#   DATA_DIR=/path/to/Segment_0001 \
#   OUTPUT_DIR=/path/to/validation-results \
#   sbatch support/SubmitScripts/ValidateModalSpacetimeInterpolator.sh
#
# The default is a Lapse-only smoke test. Set VALIDATION_COMPONENTS=all to
# validate every non-coordinate component, or provide a comma-separated list.
# Submit from the same module environment used to build SpECTRE. In particular,
# the build must use Boost 1.81 or newer. To load a project-specific environment
# in the job, set SPECTRE_ENVIRONMENT_SCRIPT to a script that defines
# spectre_load_modules.

#SBATCH --job-name=ModalSpacetimeValidation
#SBATCH --output=ModalSpacetimeValidation-%A_%a.out
#SBATCH --error=ModalSpacetimeValidation-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --partition=p.urania
#SBATCH --array=0-3
#SBATCH --no-requeue

set -euo pipefail

: "${SPECTRE_BUILD_DIR:?Set SPECTRE_BUILD_DIR to the build directory}"
: "${DATA_DIR:?Set DATA_DIR to the directory containing BbhVolume*.h5}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR for CSV, JSON, and log output}"

if [[ -n "${SPECTRE_ENVIRONMENT_SCRIPT:-}" ]]; then
  source "${SPECTRE_ENVIRONMENT_SCRIPT}"
  spectre_load_modules
fi

if [[ ! -x "${SPECTRE_BUILD_DIR}/bin/spectre" ]]; then
  echo "Cannot execute ${SPECTRE_BUILD_DIR}/bin/spectre" >&2
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p "${OUTPUT_DIR}"

shard="${SLURM_ARRAY_TASK_ID}"
input_file="${DATA_DIR}/BbhVolume${shard}.h5"
if [[ ! -r "${input_file}" ]]; then
  echo "Cannot read ${input_file}" >&2
  exit 1
fi

# Staging a shard to node-local storage makes the many small HDF5 reads much
# cheaper. Disable this with STAGE_TO_LOCAL_SCRATCH=0 if local space is tight.
stage_to_local_scratch="${STAGE_TO_LOCAL_SCRATCH:-1}"
if [[ "${stage_to_local_scratch}" == "1" ]]; then
  scratch_dir="${TMPDIR:-/tmp}/modal-spacetime-${SLURM_JOB_ID}-${shard}"
  mkdir -p "${scratch_dir}"
  cp "${input_file}" "${scratch_dir}/"
  input_file="${scratch_dir}/BbhVolume${shard}.h5"
fi

args=(
  "${SPECTRE_BUILD_DIR}/bin/spectre"
  validate-modal-spacetime-interpolator
  "${input_file}"
  --samples-per-element "${SAMPLES_PER_ELEMENT:-1}"
  --verbosity "${VALIDATION_VERBOSITY:-quiet}"
  --output "${OUTPUT_DIR}/validation-shard-${shard}.csv"
)

validation_components="${VALIDATION_COMPONENTS:-Lapse}"
if [[ "${validation_components}" != "all" ]]; then
  IFS=',' read -r -a components <<< "${validation_components}"
  for component in "${components[@]}"; do
    args+=(--var "${component}")
  done
fi
if [[ -n "${MAX_VALIDATION_OBSERVATIONS:-}" ]]; then
  args+=(--max-validation-observations "${MAX_VALIDATION_OBSERVATIONS}")
fi

echo "Running: ${args[*]}"
"${args[@]}"
