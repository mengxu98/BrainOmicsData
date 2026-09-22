"""Compare exported H5AD values/order/metadata against their source shards."""
import sys
from pathlib import Path
import anndata as ad
import pandas as pd
from scipy import io, sparse
import numpy as np

root, output = map(Path, sys.argv[1:])
manifest = pd.read_csv(root / 'expression/shard_manifest.tsv', sep='\t')
metadata = pd.read_csv(root / 'metadata/metadata.tsv.gz', sep='\t', dtype=str).set_index('Cells')
for row in manifest.itertuples(index=False):
    shard = root / 'expression' / row.Relative_Directory
    x = ad.read_h5ad(output / f'{row.Dataset}.h5ad')
    expected = sparse.csr_matrix(io.mmread(shard / 'matrix.mtx.gz').T)
    assert x.shape == expected.shape
    assert (x.X != expected).nnz == 0
    cells = pd.read_csv(shard / 'barcodes.tsv.gz', sep='\t', header=None, dtype=str)[0].tolist()
    genes = pd.read_csv(shard / 'features.tsv.gz', sep='\t', header=None, dtype=str)[0].tolist()
    assert x.obs_names.tolist() == cells
    assert x.var_names.tolist() == genes
    for column in metadata:
        a = x.obs[column].astype('string').fillna('<NA>').tolist()
        b = metadata.loc[cells, column].astype('string').fillna('<NA>').tolist()
        assert a == b, column
    for name, key in [('integrated_pca','X_integrated_pca'), ('integrated_umap','X_integrated_umap'), ('unintegrated_umap','X_unintegrated_umap')]:
        path = root / 'embeddings' / f'{name}.tsv.gz'
        if path.exists():
            ref = pd.read_csv(path, sep='\t').set_index('cell_id').loc[cells].to_numpy()
            np.testing.assert_allclose(x.obsm[key], ref, rtol=1e-12, atol=1e-12)
print('Python reader round-trip verified: counts, cell/gene order, metadata and embeddings')
