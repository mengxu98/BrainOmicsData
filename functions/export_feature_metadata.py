#!/usr/bin/env python3
"""Copy the released feature table into a package and validate its contract.

Usage: export_feature_metadata.py REPO_DIR PACKAGE_DIR
"""
import csv, pathlib, sys

EXPECTED_COLUMNS = ['Gene_ID', 'Gene_Key', 'Ensembl_IDs', 'Gene_Type']
EXPECTED_ROWS = 24659


def main():
    repo = pathlib.Path(sys.argv[1]).resolve()
    package = pathlib.Path(sys.argv[2]).resolve()
    source = repo/'data/feature_metadata.tsv'
    if not source.is_file():
        raise SystemExit(f'missing feature table: {source}')
    with open(source, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if list(rows[0].keys()) != EXPECTED_COLUMNS:
        raise SystemExit(f'feature table columns differ: {list(rows[0].keys())}')
    # Gene symbol (Gene_ID) can repeat for distinct Ensembl records; Gene_Key is unique.
    if len(rows) != EXPECTED_ROWS or len({r['Gene_Key'] for r in rows}) != EXPECTED_ROWS:
        raise SystemExit(f'feature table must hold {EXPECTED_ROWS} rows with unique Gene_Key, found {len(rows)}')
    target = package/'metadata/feature_metadata.tsv'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(source.read_bytes())
    print(f'export_feature_metadata: {len(rows)} genes -> {target.relative_to(package)}')


if __name__ == '__main__':
    main()
