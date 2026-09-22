#!/usr/bin/env python3
"""Mitochondrial fraction from the source H5AD raw/X counts (never normalized X).

Usage: percent_mito_calculate_h5ad.py RAW_ROOT OUTPUT_DIR
"""
import argparse, csv, gzip, json
from pathlib import Path
import h5py
import numpy as np
from scipy import sparse
from percent_mito_inventory_sources import strings, feature_result


def compute(path, out, dataset):
    records=[]
    with h5py.File(path) as h:
        x=h['raw/X'];var=h['raw/var'];genes=strings(var['feature_name'])
        obs=h['obs'];key=obs.attrs['_index'];key=key.decode() if isinstance(key,bytes) else key
        cells=strings(obs[key]);shape=tuple(x.attrs['shape'])
        assert shape==(len(cells),len(genes)) and len(set(cells))==len(cells)
        mt=np.char.startswith(genes.astype(str),'MT-');assert mt.any()
        encoding=x.attrs['encoding-type'];encoding=encoding.decode() if isinstance(encoding,bytes) else encoding
        assert encoding=='csr_matrix', f'Unsupported matrix encoding {encoding}: do not reinterpret axes'
        pointer=x['indptr'][:];total=np.zeros(len(cells));numerator=np.zeros(len(cells))
        for start in range(0,len(cells),4096):
            end=min(start+4096,len(cells));lo=int(pointer[start]);hi=int(pointer[end])
            values=x['data'][lo:hi];indices=x['indices'][lo:hi]
            assert np.isfinite(values).all() and (values>=0).all() and (values==np.floor(values)).all()
            block=sparse.csr_matrix((values.astype(np.float64),indices,pointer[start:end+1]-lo),shape=(end-start,len(genes)))
            total[start:end]=np.asarray(block.sum(axis=1)).ravel()
            numerator[start:end]=np.asarray(block[:,mt].sum(axis=1)).ravel()
        pct=np.full(len(cells),np.nan);np.divide(100*numerator,total,out=pct,where=total>0)
        assert ((pct[np.isfinite(pct)]>=0)&(pct[np.isfinite(pct)]<=100)).all()
        final=out/f'{dataset}_source_qc.tsv.gz';tmp=out/f'{dataset}_source_qc.partial.tsv.gz'
        with gzip.open(tmp,'wt') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['Dataset','Source_Object_Cell_ID','percent_mito','MT_Counts','Total_Counts','Source_File','Source_Layer','Status'])
            for cell,p,n,t in zip(cells,pct,numerator,total):
                w.writerow([dataset,cell,format(p,'.17g') if t>0 else '',format(n,'.17g'),format(t,'.17g'),str(path),'raw/X','COMPUTED_SOURCE_COUNTS' if t>0 else 'ZERO_TOTAL'])
        tmp.replace(final)
        row=dict(Dataset=dataset,Source_File=str(path),Kind='h5ad',Source_Bytes=path.stat().st_size,Count_Semantics='Original raw/X Gene Expression counts; exon/intron scope as supplied, not independently verified',**feature_result(genes),Source_Cells=len(cells))
        row['Status']='COMPUTED_SOURCE_COUNTS'
        with open(out/f'{dataset}_h5ad_inventory.tsv','w') as f:
            w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
        print('SOURCE_QC_COMPLETE',dataset,len(cells),int(mt.sum()),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('raw_root');p.add_argument('output');a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    for ds,name in [('Wang_2025','Wang_2025_RNA.h5ad'),('Velmeshev_2023','Velmeshev_2023_full_cellxgene.h5ad')]: compute(Path(a.raw_root)/ds/name,out,ds)
