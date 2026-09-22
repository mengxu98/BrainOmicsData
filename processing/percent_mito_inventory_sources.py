#!/usr/bin/env python3
"""Read the source feature universes, not the exported gene panel. Counts are not rewritten.

RDS/RData sources are read by inventory_r_objects.R. Unresolved Ensembl identifiers do
not by themselves indicate missing mitochondrial genes.

Usage: percent_mito_inventory_sources.py RAW_ROOT OUTPUT_TSV
"""
import argparse, csv, gzip, json, re
from pathlib import Path
import h5py
import numpy as np


def text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def strings(node):
    if isinstance(node, h5py.Group):
        cats = strings(node['categories']); codes = node['codes'][:]
        return np.array([cats[i] if i >= 0 else '' for i in codes])
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in node[:]])


def feature_result(genes):
    symbols = [str(g).strip().strip('"').split('|')[-1] for g in genes]
    mt = [g for g in symbols if g.startswith('MT-')]
    unresolved = sum(bool(re.match(r'^ENSG\d', g)) for g in symbols)
    state = 'MT_PRESENT_COUNTS_NOT_YET_VALIDATED' if mt else ('UNRESOLVED_GENE_IDS' if unresolved else 'NO_MT_FEATURES')
    return dict(Features=len(symbols), MT_Genes=len(mt), MT_Symbols=';'.join(sorted(set(mt))), Status=state)


def inventory(root):
    root = Path(root); rows = []
    def add(ds, path, kind, genes, **extra):
        row = dict(Dataset=ds, Source_File=str(path), Kind=kind, Source_Bytes=path.stat().st_size,
                   Count_Semantics='Source-supplied expression matrix; exon/intron scope not independently verified',
                   **feature_result(genes), **extra)
        rows.append(row); print(ds, path.name, row['Features'], row['MT_Genes'], row['Status'], flush=True)
    # Sources whose features are in columns: reading the header is sufficient.
    for ds, name, sep in [('AllenM1','matrix.csv',','),('GSE294786','GSE294786_all_counts_transposed.tsv.gz','\t')]:
        p=root/ds/name
        with text(p) as f: genes=next(csv.reader(f, delimiter=sep))[1:]
        add(ds,p,'cell_rows',genes,Delimiter=sep)
    # Dense gene-row files: inspect EVERY row, not a sample of rows.
    patterns = {'EGAS00001006537':'*raw_count_gEX_matrix.txt.gz', 'GSE104276':'rawData/*UMI_count_NOERCC.xls.gz',
                'GSE81475':'counts.txt', 'GSE97942':'*_counts.txt'}
    for ds, pattern in patterns.items():
        for p in sorted((root/ds).glob(pattern)):
            with text(p) as f:
                next(f); genes=[line.split(None,1)[0].strip('"') for line in f if line.strip()]
            add(ds,p,'gene_rows',genes)
    # Independent small single-cell files: no inference from a representative file.
    for p in sorted((root/'GSE67835/GSE67835').glob('*.csv')):
        with text(p) as f: genes=[line.split('\t',1)[0].strip() for line in f if line.strip()]
        add('GSE67835',p,'single_cell',genes)
    for ds, name, group in [('GSE168408','RNA-all_full-counts-and-downsampled-CPM.h5ad','X'),
                            ('Velmeshev_2023','Velmeshev_2023_full_cellxgene.h5ad','raw/X'),
                            ('Wang_2025','Wang_2025_RNA.h5ad','raw/X')]:
        p=root/ds/name
        with h5py.File(p) as h:
            var=h['raw/var' if group.startswith('raw/') else 'var']
            key=next((x for x in ['feature_name','gene_symbols','gene_name','Gene','gene','_index'] if x in var),None)
            if key is None: key=var.attrs['_index']; key=key.decode() if isinstance(key,bytes) else key
            genes=strings(var[key]); add(ds,p,'h5ad',genes,Matrix_Group=group,Gene_Field=key,Var_Fields=';'.join(var.keys()))
    # Gene Expression only: exclude ATAC peaks from numerator AND denominator.
    for ds in ['GSE217511','GSE296073','PRJCA015229','ROSMAP']:
        paths=sorted((root/ds).rglob('*features.tsv.gz'))
        for p in paths:
            if ds=='PRJCA015229' and not p.parent.parent.name.startswith('HM'): continue
            if ds=='ROSMAP' and p.parent.name!='RNA': continue
            with text(p) as f: records=list(csv.reader(f,delimiter='\t'))
            records=[r for r in records if len(r)<3 or r[2]=='Gene Expression']
            genes=[r[1] if len(r)>1 else r[0] for r in records]
            add(ds,p,'10x_features',genes)
    for ds,name in [('GSE186538','rawData/GSE186538_Human_genes.txt.gz'),('GSE207334','GSE207334_Multiome_rna_genes.txt.gz')]:
        p=root/ds/name
        with text(p) as f: genes=[line.strip() for line in f if line.strip()]
        add(ds,p,'mtx_genes',genes)
    p=root/'GSE212606/GSM6657986_gene_annotation.csv'
    with text(p) as f: genes=[r['gene_short_name'] for r in csv.DictReader(f)]
    add('GSE212606',p,'mtx_gene_annotation',genes)
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('raw_root');p.add_argument('output');a=p.parse_args()
    rows=inventory(a.raw_root); out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    columns=list(dict.fromkeys(k for r in rows for k in r))
    tmp=out.with_suffix('.partial')
    with open(tmp,'w') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
    tmp.replace(out)
    print('FEATURE_INVENTORY_COMPLETE',len(rows),flush=True)
if __name__=='__main__': main()
