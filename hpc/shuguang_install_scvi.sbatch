#!/bin/bash
#SBATCH -J brainomics_scvi_env
#SBATCH -p cpu-partition
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/install_scvi_%j.out
#SBATCH -e /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/install_scvi_%j.err

set -euo pipefail

repo_dir="/path/to/hpc/home/repositories/BrainOmicsData"
conda_prefix="/path/to/hpc/home/miniforge3"
scvi_prefix="/path/to/hpc/home/envs/brainomics-scvi"
wheelhouse="/path/to/hpc/home/software_sources/python/scvi-tools-1.5.0.post1-cp312-linux-x86_64"
runtime_lib_dir="/path/to/hpc/home/software_sources/lib/brainomics-scvi-runtime"
scvi_version="1.5.0.post1"

if [ ! -f "$runtime_lib_dir/runtime.sha256" ]; then
  echo "Missing pinned scVI runtime libraries: $runtime_lib_dir/runtime.sha256" >&2
  exit 1
fi
(
  cd "$runtime_lib_dir"
  sha256sum -c runtime.sha256
)
export LD_LIBRARY_PATH="$runtime_lib_dir:${LD_LIBRARY_PATH:-}"

mkdir -p "$repo_dir/hpc/logs" "$(dirname "$scvi_prefix")"

if [ ! -f "$wheelhouse/wheelhouse.sha256" ]; then
  echo "Missing offline scVI wheelhouse: $wheelhouse" >&2
  exit 1
fi
(
  cd "$wheelhouse"
  sha256sum -c wheelhouse.sha256
)

if [ ! -x "$scvi_prefix/bin/python" ]; then
  "$conda_prefix/bin/python3" -m venv "$scvi_prefix"
fi

"$scvi_prefix/bin/python" -m pip install \
  --no-index \
  --find-links "$wheelhouse" \
  --upgrade \
  --requirement "$wheelhouse/requirements.txt"

imaging_library="$($scvi_prefix/bin/python - <<'PY'
import pathlib
import PIL

matches = sorted(pathlib.Path(PIL.__file__).parent.glob("_imaging.*.so"))
if len(matches) != 1:
    raise RuntimeError(f"Expected one Pillow imaging library, found {matches}")
print(matches[0])
PY
)"
if ldd "$imaging_library" | grep -q "not found"; then
  ldd "$imaging_library" >&2
  exit 1
fi

"$scvi_prefix/bin/python" - <<'PY'
import anndata
import numpy
import pandas
import scipy
import scvi
import torch

assert scvi.__version__ == "1.5.0.post1"
assert tuple(map(int, __import__("sys").version_info[:2])) == (3, 12)
print("scvi-tools", scvi.__version__)
print("torch", torch.__version__)
print("CUDA build", torch.version.cuda)
print("anndata", anndata.__version__)
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("scipy", scipy.__version__)
PY

"$scvi_prefix/bin/python" -m pip freeze \
  > "$scvi_prefix/brainomics-scvi-freeze.txt"
