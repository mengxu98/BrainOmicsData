#!/usr/bin/env python3
"""Summarize the per-source mitochondrial inventories for the released cohort.

Reads the component-level inventories written by the calculators, the cohort table
(Dataset and retained cells) and the per-source QC files, then writes one row per
dataset plus the component evidence.

Usage: percent_mito_summarize_inventory.py QC_DIR COHORT_TSV OUT_TSV
"""
import argparse, csv
from pathlib import Path


def read(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    columns = list(dict.fromkeys(k for r in rows for k in r))
    tmp = path.with_suffix('.partial')
    with open(tmp, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def component_rows(qc):
    components = []
    for path in sorted(qc.glob('*inventory.tsv')):
        components += read(path)
    for path in sorted((qc / 'gse296073_10x').glob('*inventory.tsv')):
        components += read(path)
    return components


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('qc')
    parser.add_argument('cohort')
    parser.add_argument('out')
    args = parser.parse_args()
    qc = Path(args.qc)
    components = component_rows(qc)
    rows = []
    for entry in read(Path(args.cohort)):
        dataset = entry['Dataset']
        found = [x for x in components if x['Dataset'] == dataset]
        mt = [int(x['MT_Genes']) for x in found]
        qc_file = qc / f'{dataset}_source_qc.tsv.gz'
        if dataset == 'GSE296073':
            qc_file = qc / 'gse296073_10x' / f'{dataset}_source_qc.tsv.gz'
        if not found:
            state = 'PENDING_SOURCE_INSPECTION'
        elif any(x['Status'] == 'UNRESOLVED_GENE_IDS' for x in found):
            state = 'UNRESOLVED_GENE_IDS'
        elif max(mt) == 0:
            state = 'NO_MT_FEATURES'
        elif qc_file.exists():
            state = 'SOURCE_QC_COMPUTED'
        else:
            state = 'MT_PRESENT_COMPUTATION_PENDING'
        rows.append(dict(
            Dataset=dataset,
            Retained_Cells=entry.get('Cells', ''),
            Source_Files=';'.join(sorted({x['Source_File'] for x in found})),
            Source_Components=len(found),
            MT_Genes_Min=min(mt) if mt else '',
            MT_Genes_Max=max(mt) if mt else '',
            Status=state,
            Formula='100 * source MT counts / source total Gene Expression counts',
            MT_Definition='Source gene symbol starts with MT-; nuclear pseudogenes excluded',
            Count_Semantics='Exon/intron scope retained as supplied; no re-quantification from reads',
            NA_Policy='Blank for sources without MT features or with zero total counts; never zero-filled',
            Note='Source stored under raw/GSE207334/Ma_Sestan_mat.rds' if dataset == 'Ma_et_al_2022' else ''
        ))
    write(Path(args.out), rows)
    for row in rows:
        print(row['Dataset'], row['MT_Genes_Max'], row['Status'])


if __name__ == '__main__':
    main()
