#!/usr/bin/env python3
"""Mitochondrial fraction from the original GEO CellRanger H5 files.

The author H5AD for this source carries no MT genes, so the GEO libraries are used.
Library identity comes from the canonical source metadata; barcodes are never matched
across libraries and unavailable barcodes stay blank.

Usage: percent_mito_calculate_gse168408.py DATA_ROOT OUTPUT_DIR
"""
import csv,gzip,io,argparse
from pathlib import Path
import h5py,numpy as np
from scipy import sparse
from percent_mito_inventory_sources import strings,feature_result
from percent_mito_summarize_inventory import write
p=argparse.ArgumentParser();p.add_argument('data_root');p.add_argument('out');a=p.parse_args()
root=Path(a.data_root);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
with gzip.open(root/'processed/GSE168408/metadata_canonical.tsv.gz','rt') as f:
    rows=[{k:r[k] for k in ['Original_Cell_ID','Original_Library_ID','Original_Source_Record_ID']} for r in csv.DictReader(f,delimiter='\t')]
files=list((root/'raw/GSE168408').rglob('*_raw_feature_bc_matrix.h5.gz'))
bygeo={x.name.split('_')[0]:x for x in files};assert len(bygeo)==len(files)
geos=sorted({r['Original_Source_Record_ID'] for r in rows});assert len(geos)==27 and set(geos)<=set(bygeo)
inventory=[];all_scores=[];missing=[]
for geo in geos:
    path=bygeo[geo];target=[r for r in rows if r['Original_Source_Record_ID']==geo]
    libraries={r['Original_Library_ID'] for r in target};assert len(libraries)==1
    library=next(iter(libraries));suffix='-'+library
    assert all(r['Original_Cell_ID'].endswith(suffix) for r in target)
    wanted=[r['Original_Cell_ID'][:-len(suffix)]+'-1' for r in target]
    assert len(set(wanted))==len(wanted)
    with gzip.open(path,'rb') as f: buf=io.BytesIO(f.read())
    with h5py.File(buf) as h:
        g=h['matrix'];genes=strings(g['features/name']);bars=strings(g['barcodes'])
        ids={b:i for i,b in enumerate(bars)};assert len(ids)==len(bars)
        keep=np.ones(len(genes),dtype=bool)
        if 'feature_type' in g['features']:keep=strings(g['features/feature_type'])=='Gene Expression'
        mt=np.char.startswith(genes.astype(str),'MT-')&keep;assert mt.any()
        found=[i for i,b in enumerate(wanted) if b in ids];indices=[ids[wanted[i]] for i in found]
        x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))[:,indices].astype(np.float64)
        assert np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all()
        total=np.asarray(x[keep,:].sum(axis=0)).ravel();numerator=np.asarray(x[mt,:].sum(axis=0)).ravel()
        pct=np.full(len(found),np.nan);np.divide(100*numerator,total,out=pct,where=total>0)
        assert ((pct[np.isfinite(pct)]>=0)&(pct[np.isfinite(pct)]<=100)).all()
        lookup={j:(pct[k],numerator[k],total[k]) for k,j in enumerate(found)}
        for i,r in enumerate(target):
            present=i in lookup;pct_i,num,tot=lookup.get(i,(np.nan,np.nan,np.nan))
            state='SOURCE_CELL_ABSENT' if not present else 'ZERO_TOTAL' if tot==0 else 'COMPUTED_SOURCE_COUNTS'
            rec=dict(Dataset='GSE168408',Source_Object_Cell_ID=r['Original_Cell_ID'],percent_mito='' if not np.isfinite(pct_i) else format(pct_i,'.17g'),MT_Counts='' if not present else num,Total_Counts='' if not present else tot,Source_File=str(path),Source_Layer='CellRanger matrix; Gene Expression only',Status=state)
            all_scores.append(rec)
            if not present:missing.append(dict(Original_Cell_ID=r['Original_Cell_ID'],GEO=geo,Source_Barcode=wanted[i]))
        inv=dict(Dataset='GSE168408',Source_File=str(path),Kind='original_10x_h5',Source_Bytes=path.stat().st_size,Count_Semantics='Original CellRanger Gene Expression counts; exon/intron scope as supplied, not independently verified',**feature_result(genes[keep]),Source_Cells=len(bars),Requested_Cells=len(target),Matched_Cells=len(found))
        inv['Status']='COMPUTED_SOURCE_COUNTS';inventory.append(inv)
        write(out/'GSE168408_10x_inventory.tsv',inventory)
        print(geo,'MT',int(mt.sum()),'matched',len(found),'/',len(target),flush=True)
    del x,buf
assert len(all_scores)==len(rows) and len({r['Source_Object_Cell_ID'] for r in all_scores})==len(rows)
final=out/'GSE168408_source_qc.tsv.gz';tmp=out/'GSE168408_source_qc.partial.tsv.gz'
with gzip.open(tmp,'wt') as f:
    w=csv.DictWriter(f,fieldnames=list(all_scores[0]),delimiter='\t');w.writeheader();w.writerows(all_scores)
tmp.replace(final)
if missing:write(out/'GSE168408_cells_absent_from_10x.tsv',missing)
print('SOURCE_QC_COMPLETE GSE168408',len(rows),'missing',len(missing),flush=True)
