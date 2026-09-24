#!/usr/bin/env python3
"""Contract check for a released ScienceDB package directory.

Usage: test_package_contract.py [--strict] [PACKAGE_DIR]
Default: $BRAINOMICS_PACKAGE_DIR or <repo>/../../data/BrainOmicsData/ScienceDB.
Exits 0 with a SKIP message when no package directory is present.
"""
import csv, gzip, hashlib, os, pathlib, sys

REQUIRED_MANIFEST_COLUMNS = ['Dataset', 'Cells', 'Donors', 'Cross_Dataset_Donors', 'Specimens', 'Libraries',
                             'Original_CellType_Source', 'Source_Accession', 'Source_Repository', 'Publication_DOI',
                             'Title', 'First_Author', 'Publication_Year', 'Source_URL', 'Repository_Record_URL',
                             'Processed_Input_Route', 'Processed_Input_URLs', 'Processed_Record_URL',
                             'Data_License', 'Derived_Matrix_Scope', 'Article_License_URL']
FILE_MANIFEST_COLUMNS = ['Relative_Path', 'Bytes', 'Dimensions', 'Schema_Version', 'License', 'SHA256', 'MD5']
FORBIDDEN_PATHS = ['metadata/cluster_annotation.tsv', 'provenance/celltype_raw_dictionary.tsv',
                   'Supplementary_Data.xlsx', 'reproduction.zip', 'objects', 'metadata/dataset_summary.tsv']
FORBIDDEN_WORDS = ['Current_', '2,712,452', '2712452']


def main():
    repo = pathlib.Path(__file__).resolve().parent.parent
    default = repo.parents[1]/'data/BrainOmicsData/ScienceDB'
    args = [a for a in sys.argv[1:] if a != '--strict']
    strict = '--strict' in sys.argv[1:]
    root = pathlib.Path(args[0] if args else os.environ.get('BRAINOMICS_PACKAGE_DIR', default)).resolve()
    if not root.is_dir():
        print(f'SKIP package contract: no package at {root}')
        return 0
    problems = []
    for rel in FORBIDDEN_PATHS:
        if (root/rel).exists():
            problems.append(f'forbidden path present: {rel}')
    # manifest + checksum chain
    manifest = root/'provenance/file_manifest.tsv'
    if not manifest.is_file():
        problems.append('missing provenance/file_manifest.tsv')
    else:
        with open(manifest, newline='') as handle:
            file_rows = list(csv.DictReader(handle, delimiter='\t'))
        if list(file_rows[0].keys()) != FILE_MANIFEST_COLUMNS:
            problems.append('file_manifest.tsv columns differ')
        if len(file_rows) != 82 or len({row['Relative_Path'] for row in file_rows}) != len(file_rows):
            problems.append(f'file_manifest.tsv must list 82 unique payload files, found {len(file_rows)}')
        for row in file_rows:
            rel = row['Relative_Path']
            if not row['License']:
                problems.append(f'missing licence in file manifest: {rel}')
            if rel not in ('README.md', 'scripts/readers.zip'):
                if row['Schema_Version'] != '2.0.0' or not row['Dimensions']:
                    problems.append(f'missing dimensions or schema version: {rel}')
        not_mirrored = []
        for row in file_rows:
            path = root/row['Relative_Path']
            if not path.is_file():
                not_mirrored.append(row['Relative_Path'])
                if strict:
                    problems.append(f'manifest lists missing file: {row["Relative_Path"]}')
                continue
            if hashlib.sha256(path.read_bytes()).hexdigest() != row['SHA256']:
                problems.append(f'sha256 mismatch: {row["Relative_Path"]}')
        if not_mirrored:
            print(f'note: {len(not_mirrored)} large payload files are not mirrored locally')
    # dataset manifest
    dataset_manifest = root/'provenance/dataset_manifest.tsv'
    if dataset_manifest.is_file():
        with open(dataset_manifest, newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            header, rows = reader.fieldnames, list(reader)
        missing = [c for c in REQUIRED_MANIFEST_COLUMNS if c not in header]
        if missing:
            problems.append('dataset_manifest.tsv missing columns: ' + ', '.join(missing))
        if len(rows) != 22 or len({r['Dataset'] for r in rows}) != 22:
            problems.append(f'dataset_manifest.tsv must hold 22 datasets, found {len(rows)}')
        if manifest.is_file():
            licenses = {r['Dataset']: r['Data_License'] for r in rows}
            for item in file_rows:
                parts = pathlib.PurePosixPath(item['Relative_Path']).parts
                if len(parts) == 4 and parts[:2] == ('expression', 'shards'):
                    if item['License'] != licenses.get(parts[2]):
                        problems.append(f'expression shard licence differs from source manifest: {item["Relative_Path"]}')
    else:
        problems.append('missing provenance/dataset_manifest.tsv')
    shard_manifest = root/'expression/shard_manifest.tsv'
    if shard_manifest.is_file() and manifest.is_file():
        with open(shard_manifest, newline='') as handle:
            shards = {row['Dataset']: row for row in csv.DictReader(handle, delimiter='\t')}
        if len(shards) != 22:
            problems.append(f'shard_manifest.tsv must hold 22 datasets, found {len(shards)}')
        for item in file_rows:
            parts = pathlib.PurePosixPath(item['Relative_Path']).parts
            if len(parts) != 4 or parts[:2] != ('expression', 'shards') or parts[2] not in shards:
                continue
            source = shards[parts[2]]
            cells, features = int(source['Cells']), int(source['Features'])
            expected = {'matrix.mtx.gz': f'{features}x{cells}',
                        'features.tsv.gz': f'{features}x3',
                        'barcodes.tsv.gz': f'{cells}x1'}.get(parts[3])
            if item['Dimensions'] != expected:
                problems.append(f'expression shard dimensions differ from shard manifest: {item["Relative_Path"]}')
    # metadata
    metadata = root/'metadata/metadata.tsv.gz'
    if metadata.is_file():
        with gzip.open(metadata, 'rt') as handle:
            columns = handle.readline().rstrip('\n').split('\t')
            cells = sum(1 for _ in handle)
        if len(columns) != 21:
            problems.append(f'metadata.tsv.gz must hold 21 columns, found {len(columns)}')
        if cells != 2602031:
            problems.append(f'metadata.tsv.gz must hold 2,602,031 cells, found {cells}')
        if 'Sample' in columns:
            problems.append('metadata.tsv.gz still carries the redundant Sample column')
        if columns[4:6] != ['Sample_ID', 'Original_Sample_ID']:
            problems.append(f'Original_Sample_ID must follow Sample_ID, found {columns[4:6]}')
        if columns[-1] != 'Original_CellType':
            problems.append(f'last metadata column must be Original_CellType, found {columns[-1]}')
    else:
        problems.append('missing metadata/metadata.tsv.gz')
    # dictionary and feature table
    dictionary = root/'metadata/metadata_dictionary.tsv'
    if dictionary.is_file():
        with open(dictionary, newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if list(rows[0].keys()) != ['Field', 'Definition'] or len(rows) != 21:
            problems.append(f'metadata_dictionary.tsv must hold 21 Field/Definition rows, found {len(rows)}')
    features = root/'metadata/feature_metadata.tsv'
    if features.is_file():
        with open(features, newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if len(rows) != 24659 or len({r['Gene_Key'] for r in rows}) != 24659:
            problems.append(f'feature_metadata.tsv must hold 24,659 unique Gene_Key rows, found {len(rows)}')
    # forbidden wording in text tables and README
    for rel in ['README.md', 'metadata/metadata_dictionary.tsv', 'provenance/dataset_manifest.tsv']:
        path = root/rel
        if path.is_file():
            text = path.read_text(errors='ignore')
            for word in FORBIDDEN_WORDS:
                if word in text:
                    problems.append(f'{rel} contains forbidden wording: {word}')
    if problems:
        print('FAIL package contract')
        for item in problems[:20]:
            print('  -', item)
        return 1
    print(f'PASS package contract{" (strict)" if strict else ""}: manifest hashes, table contracts, forbidden paths and wording')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
