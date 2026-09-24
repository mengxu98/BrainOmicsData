"""Check the Python version and installed package versions."""
import argparse
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--imports', action='store_true', help='Import the main analysis packages.')
args = parser.parse_args()
here = Path(__file__).resolve().parent
if sys.version_info[:3] != (3, 12, 8):
    raise RuntimeError(f'Python 3.12.8 is required; observed {sys.version}')
normalize = lambda s: s.lower().replace('_', '-').replace('.', '-')
expected = {}
for line in (here / 'python-scvi-freeze.txt').read_text().splitlines():
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    name, version = line.split('==')
    if normalize(name) in expected:
        raise ValueError(f'Duplicate package in Python lock: {name}')
    expected[normalize(name)] = version
actual = {normalize(d.metadata['Name']): d.version for d in metadata.distributions()}
errors = {n: {'expected': v, 'observed': actual.get(n)} for n, v in expected.items() if actual.get(n) != v}
extras = sorted(set(actual) - set(expected))
if errors or extras:
    raise RuntimeError({'version_errors': errors, 'unexpected_packages': extras})
result = {'state': 'ready', 'python': sys.version,
          'packages': len(expected), 'prefix': sys.prefix}
if args.imports:
    import anndata
    import numpy
    import pandas
    import scipy
    import scvi
    import torch
    assert scvi.__version__ == '1.5.0.post1'
    assert torch.__version__ == '2.5.1+cu124' and torch.version.cuda == '12.4'
    result.update(imports='ready', torch=torch.__version__,
                  cuda_build=torch.version.cuda)
print(json.dumps(result, indent=2))
