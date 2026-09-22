"""Verify pinned package versions; optional real imports, never training."""
import argparse
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--imports', action='store_true', help='Import the six analysis packages; does not run analysis.')
args = parser.parse_args()
here = Path(__file__).resolve().parent
assert sys.version_info[:3] == (3, 12, 8), sys.version
normalize = lambda s: s.lower().replace('_', '-').replace('.', '-')
expected = {}
for line in (here / 'python-scvi-versions.txt').read_text().splitlines():
    name, version = line.split('==')
    expected[normalize(name)] = version
actual = {normalize(d.metadata['Name']): d.version for d in metadata.distributions()}
errors = {n: {'expected': v, 'observed': actual.get(n)} for n, v in expected.items() if actual.get(n) != v}
extras = sorted(set(actual) - set(expected) - {'pip', 'wheel'})
assert not errors and not extras, {'version_errors': errors, 'unexpected_packages': extras}
result = {'state': 'PASS_PINNED_PYTHON_VERSIONS', 'python': sys.version,
          'analysis_packages': len(expected), 'prefix': sys.prefix,
          'scope': 'Installed package versions only; no data read or computation.'}
if args.imports:
    import anndata
    import numpy
    import pandas
    import scipy
    import scvi
    import torch
    assert scvi.__version__ == '1.5.0.post1'
    assert torch.__version__ == '2.5.1+cu124' and torch.version.cuda == '12.4'
    result.update(state='PASS_PINNED_VERSIONS_AND_IMPORTS', torch=torch.__version__, cuda_build=torch.version.cuda,
                  scope='Package versions and actual imports; CUDA build metadata does not establish GPU availability or successful training. No data read or computation.')
print(json.dumps(result, indent=2))
