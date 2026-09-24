#!/usr/bin/env python3
"""Bring an exported deposit package to the released structure.

Every step below corresponds to a released change; running it on an already
normalized package is a no-op apart from regenerating the checksums.

Usage: normalize_package.py PACKAGE_DIR
"""
import csv, datetime, hashlib, json, pathlib, sys

FILE_MANIFEST_COLUMNS = ['Relative_Path', 'Bytes', 'Dimensions', 'Schema_Version', 'License', 'SHA256', 'MD5']
PACKAGE_SCHEMA_VERSION = '2.0.0'

DROP_STATUS_COLUMNS = ['Public_Deposition_Status', 'License_Scope', 'Article_License_Scope', 'Redistribution_Status']
DROP_MT_COLUMNS = ['Source_Files', 'Source_Component_Count', 'MT_Gene_Count_Min', 'MT_Gene_Count_Max',
                   'MT_Genes', 'Percent_Mito_Status', 'Percent_Mito_Definition', 'Percent_Mito_NA_Policy']
REQUIRED_MANIFEST_COLUMNS = ['Dataset', 'Cells', 'Donors', 'Specimens', 'Libraries', 'Original_CellType_Source',
                             'Source_Accession', 'Source_Repository', 'Publication_DOI', 'Title', 'First_Author',
                             'Publication_Year', 'Source_URL', 'Repository_Record_URL', 'Processed_Input_Route',
                             'Processed_Input_URLs', 'Processed_Record_URL', 'Data_License',
                             'Derived_Matrix_Scope', 'Article_License_URL']
FIELD_DEFINITIONS = {
 'Cells': 'Public cell or nucleus identifier',
 'Dataset': 'Source dataset identifier',
 'Technology': 'Source sequencing technology',
 'Sequence': 'Source assay or modality string as reported; predominantly scRNA-seq or snRNA-seq, including RNA arms of multiome assays',
 'Sample_ID': 'Public specimen identifier',
 'Original_Sample_ID': 'Source sample or specimen identifier; use with Dataset',
 'Donor_ID': 'Public donor identifier; shared donors reconciled across sources',
 'Library_ID': 'Public library identifier; blank when unavailable',
 'BrainRegion': 'Standardized anatomical region',
 'Age': 'Age standardized from the source value and unit where available, otherwise the source-reported age string; unitless values and ranges retained',
 'Sex': 'Source-standardized sex; unresolved values are left blank',
 'Cluster': 'RPCA cluster, C0-C74',
 'seurat_clusters': 'Numeric part of Cluster, 0-74',
 'CellType': 'Adopted annotation, 12 major types',
 'AgeIntervalID': 'Age interval, S1-S15',
 'AgeInterval': 'Developmental or adult interval name',
 'AgeRange': 'Definition of the age interval',
 'percent_mito': '100 x mitochondrial gene counts / all Gene Expression counts in the original source matrix, joined by source library and cell identifier; blank when the source has no MT features, the source cell is absent, or the source total is zero',
 'n_counts': 'Total counts in the released 24,659-gene matrix',
 'n_genes': 'Detected genes in the released matrix',
 'Original_CellType': 'Source-study label without harmonization; blank when not reported',
}


def read_tsv(path):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return reader.fieldnames, list(reader)


def write_tsv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, '') for column in columns})


def manifest_attributes(root, rel, shards, licenses, cell_count):
    """Describe the released file axes and the applicable licence evidence."""
    dimensions = ''
    license_text = 'MIT'
    parts = pathlib.PurePosixPath(rel).parts
    if len(parts) == 4 and parts[:2] == ('expression', 'shards'):
        dataset, component = parts[2:]
        shard = shards[dataset]
        cells, features = int(shard['Cells']), int(shard['Features'])
        if component == 'matrix.mtx.gz':
            dimensions = f'{features}x{cells}'
        elif component == 'features.tsv.gz':
            dimensions = f'{features}x3'
        elif component == 'barcodes.tsv.gz':
            dimensions = f'{cells}x1'
        else:
            raise ValueError(f'unrecognized expression component: {rel}')
        license_text = licenses[dataset]
    elif rel == 'expression/shard_manifest.tsv':
        dimensions = f'{len(shards)}x5'
    elif rel == 'metadata/metadata.tsv.gz':
        dimensions = f'{cell_count}x21'
    elif rel == 'metadata/metadata_dictionary.tsv':
        dimensions = '21x2'
    elif rel == 'metadata/feature_metadata.tsv':
        dimensions = f'{next(iter(shards.values()))["Features"]}x4'
    elif rel == 'provenance/dataset_manifest.tsv':
        with open(root/rel, newline='') as handle:
            columns = next(csv.reader(handle, delimiter='\t'))
        dimensions = f'{len(licenses)}x{len(columns)}'
    elif rel.startswith('embeddings/') and rel.endswith('.tsv.gz'):
        coordinates = 2 if rel.endswith('_umap.tsv.gz') else 50
        dimensions = f'{cell_count}x{coordinates}'
    elif rel == 'validation/lisi.tsv.gz':
        dimensions = f'{cell_count}x8'
    elif rel not in ('README.md', 'scripts/readers.zip'):
        raise ValueError(f'unrecognized package file: {rel}')
    schema = '' if rel in ('README.md', 'scripts/readers.zip') else PACKAGE_SCHEMA_VERSION
    return {'Dimensions': dimensions, 'Schema_Version': schema, 'License': license_text}


def enrich_manifest(root, rows):
    _, shard_rows = read_tsv(root/'expression/shard_manifest.tsv')
    _, source_rows = read_tsv(root/'provenance/dataset_manifest.tsv')
    shards = {row['Dataset']: row for row in shard_rows}
    licenses = {row['Dataset']: row['Data_License'] for row in source_rows}
    if len(shards) != 22 or set(shards) != set(licenses):
        raise ValueError('source datasets differ between shard and provenance manifests')
    cell_count = sum(int(row['Cells']) for row in shard_rows)
    if cell_count != 2602031 or {int(row['Features']) for row in shard_rows} != {24659}:
        raise ValueError('shard dimensions differ from the released cohort')
    for row in rows:
        row.update(manifest_attributes(root, row['Relative_Path'], shards, licenses, cell_count))


def update_manifest_metadata_only(root):
    """Refresh manifest metadata without hashing unchanged large payloads."""
    manifest_path = root/'provenance/file_manifest.tsv'
    columns, rows = read_tsv(manifest_path)
    if not set(['Relative_Path', 'Bytes', 'SHA256', 'MD5']).issubset(columns):
        raise ValueError('existing file manifest lacks checksum columns')
    enrich_manifest(root, rows)
    for row in rows:
        rel = row['Relative_Path']
        path = root/rel
        if rel == 'README.md':
            data = path.read_bytes()
            row['Bytes'] = len(data)
            row['SHA256'] = hashlib.sha256(data).hexdigest()
            row['MD5'] = hashlib.md5(data).hexdigest()
        elif path.is_file() and path.stat().st_size != int(row['Bytes']):
            raise ValueError(f'payload size changed; full normalization required: {rel}')
    write_tsv(manifest_path, FILE_MANIFEST_COLUMNS, rows)
    md5_lines = [f'{row["MD5"]}  {row["Relative_Path"]}' for row in rows]
    md5_lines.append(f'{hashlib.md5(manifest_path.read_bytes()).hexdigest()}  provenance/file_manifest.tsv')
    (root/'md5sum.txt').write_text('\n'.join(md5_lines) + '\n')
    print(json.dumps({'payload_files': len(rows), 'schema_version': PACKAGE_SCHEMA_VERSION}))


def main():
    root = pathlib.Path(sys.argv[1]).resolve()
    if '--manifest-metadata-only' in sys.argv[2:]:
        update_manifest_metadata_only(root)
        return
    seal_path = pathlib.Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else root.parent/'package_seal.json'
    changed = []

    # 1. cluster annotation duplicates Cluster/CellType/counts already in metadata.
    cluster_annotation = root/'metadata/cluster_annotation.tsv'
    if cluster_annotation.is_file():
        cluster_annotation.unlink(); changed.append('removed metadata/cluster_annotation.tsv')

    # 2. cell-type label dictionary: keep only the per-dataset source column.
    dictionary = root/'provenance/celltype_raw_dictionary.tsv'
    manifest_path = root/'provenance/dataset_manifest.tsv'
    columns, rows = read_tsv(manifest_path)
    if dictionary.is_file():
        _, label_rows = read_tsv(dictionary)
        summary = {}
        for row in label_rows:
            entry = summary.setdefault(row['Dataset'], {'labels': set(), 'column': set(), 'semantics': set()})
            entry['labels'].add(row['Original_CellType']); entry['column'].add(row['Source_Column']); entry['semantics'].add(row['Source_Semantics'])
        if 'Original_CellType_Source' not in columns:
            columns.insert(columns.index('Source_Accession'), 'Original_CellType_Source')
        for row in rows:
            entry = summary.get(row['Dataset'])
            row['Original_CellType_Source'] = '' if not entry else (
                f"{'; '.join(sorted(entry['column']))} ({len(entry['labels'])} labels); {'; '.join(sorted(entry['semantics']))}")
        dictionary.unlink(); changed.append('merged provenance/celltype_raw_dictionary.tsv into dataset_manifest.tsv')

    # 3. drop columns already covered by metadata or by the manifest itself.
    drop = [c for c in columns if c in DROP_STATUS_COLUMNS or c in DROP_MT_COLUMNS]
    if drop:
        columns = [c for c in columns if c not in drop]
        changed.append('dropped columns: ' + ', '.join(drop))
    missing = [c for c in REQUIRED_MANIFEST_COLUMNS if c not in columns]
    if missing:
        raise SystemExit('dataset_manifest.tsv is missing required columns: ' + ', '.join(missing))
    write_tsv(manifest_path, columns, rows)

    # 4. metadata dictionary: generated from the released field definitions.
    dictionary_path = root/'metadata/metadata_dictionary.tsv'
    write_tsv(dictionary_path, ['Field', 'Definition'],
              [{'Field': field, 'Definition': definition} for field, definition in FIELD_DEFINITIONS.items()])
    changed.append('generated metadata/metadata_dictionary.tsv')

    # 5. checksum files live under provenance/, md5sum.txt stays at the package root.
    root_manifest = root/'file_manifest.tsv'
    if root_manifest.is_file():
        target = root/'provenance/file_manifest.tsv'
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(root_manifest.read_bytes()); root_manifest.unlink()
        changed.append('moved file_manifest.tsv into provenance/')
    manifest_path = root/'provenance/file_manifest.tsv'
    excluded = {'md5sum.txt', 'provenance/file_manifest.tsv'}
    payload = sorted(str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and str(p.relative_to(root)) not in excluded)
    manifest_rows, md5_lines = [], []
    for rel in payload:
        data = (root/rel).read_bytes()
        digest256, digest5 = hashlib.sha256(data).hexdigest(), hashlib.md5(data).hexdigest()
        manifest_rows.append({'Relative_Path': rel, 'Bytes': len(data), 'SHA256': digest256, 'MD5': digest5})
        md5_lines.append(f'{digest5}  {rel}')
    enrich_manifest(root, manifest_rows)
    write_tsv(manifest_path, FILE_MANIFEST_COLUMNS, manifest_rows)
    md5_lines.append(f"{hashlib.md5(manifest_path.read_bytes()).hexdigest()}  provenance/file_manifest.tsv")
    (root/'md5sum.txt').write_text('\n'.join(md5_lines) + '\n')
    seal = {'state': 'PASS', 'files_hashed': len(manifest_rows), 'sha256': 'in provenance/file_manifest.tsv',
            'md5': 'md5sum.txt', 'manifest': 'provenance/file_manifest.tsv',
            'sealed_at': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
            'note': 'not part of the deposit: reproduction.zip, Supplementary_Data.xlsx, objects/ (plotting object removed)'}
    seal_path.write_text(json.dumps(seal) + '\n')
    print(json.dumps({'changed': changed, 'payload_files': len(manifest_rows)}, indent=1))


if __name__ == '__main__':
    main()
