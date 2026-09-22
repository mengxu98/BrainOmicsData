import sys,tempfile,gzip,csv
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'processing'))
import h5py,numpy as np
from scipy import sparse
from percent_mito_calculate_h5ad import compute
root=Path(tempfile.mkdtemp(prefix='mito_h5ad_fixture_'));path=root/'source.h5ad'
x=sparse.csr_matrix(np.array([[1,2,3,4],[0,1,1,8],[0,0,0,0]],dtype=float))
with h5py.File(path,'w') as h:
 g=h.create_group('raw/X');g.attrs['shape']=x.shape;g.attrs['encoding-type']='csr_matrix'
 for k,v in [('data',x.data),('indices',x.indices),('indptr',x.indptr)]:g.create_dataset(k,data=v)
 h.create_dataset('raw/var/feature_name',data=np.array(['MT-ND1','ICMT-DT','MTND1P23','GENE'],dtype='S'))
 h.create_dataset('obs/_index',data=np.array(['A','B','C'],dtype='S'));h['obs'].attrs['_index']='_index'
 h.create_dataset('X',data=np.ones((3,4))*999)
compute(path,root,'fixture')
with gzip.open(root/'fixture_source_qc.tsv.gz','rt') as f:r=list(csv.DictReader(f,delimiter='\t'))
assert [z['percent_mito'] for z in r]==['10','0','']
assert [z['Status'] for z in r]==['COMPUTED_SOURCE_COUNTS','COMPUTED_SOURCE_COUNTS','ZERO_TOTAL']
with h5py.File(path,'r+') as h:h['raw/X/data'][0]=.5
try:compute(path,root,'invalid')
except AssertionError:pass
else:raise AssertionError('fractional counts accepted')
print('PASS: H5AD raw/X not normalized X; full denominator; pseudogenes excluded; zero total NA; fractional rejection')
