#!/usr/bin/env python3
"""Fill the reference and licence columns of provenance/dataset_manifest.tsv.

Sources are the repository tables data/source_access_summary.tsv (dataset and
citation fields) and data/source_redistribution_evidence.tsv (processed-input
routes and licence wording).  The evidence table holds the released wording, so
every value is copied verbatim.

Usage: merge_dataset_manifest.py REPO_DIR PACKAGE_DIR
"""
import csv, pathlib, re, sys


def source_first_authors(repo):
    """Parse the single first-author mapping kept in sciencedb/package_metadata.R."""
    text = (repo/'sciencedb/package_metadata.R').read_text()
    block = re.search(r'brainomics_source_first_authors <- function\(\) \{(.*?)\n\}', text, re.S)
    if not block:
        raise SystemExit('cannot locate brainomics_source_first_authors()')
    return dict(re.findall(r'([A-Za-z0-9_]+)\s*=\s*"([^"]+)"', block.group(1)))


def read_tsv(path):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return reader.fieldnames, list(reader)


def write_tsv(path, columns, rows):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, '') for column in columns})


def main():
    repo = pathlib.Path(sys.argv[1]).resolve()
    package = pathlib.Path(sys.argv[2]).resolve()
    access_columns, access = read_tsv(repo/'data/source_access_summary.tsv')
    evidence_columns, evidence = read_tsv(repo/'data/source_redistribution_evidence.tsv')
    authors = source_first_authors(repo)
    access_by_dataset = {r['dataset']: r for r in access}
    evidence_by_dataset = {r['Dataset']: r for r in evidence}
    manifest_path = package/'provenance/dataset_manifest.tsv'
    columns, rows = read_tsv(manifest_path)
    for row in rows:
        dataset = row['Dataset']
        a = access_by_dataset.get(dataset)
        e = evidence_by_dataset.get(dataset)
        if a is None or e is None:
            raise SystemExit(f'source tables lack dataset {dataset}')
        row['Source_Accession'] = a['source_accession']
        row['Source_Repository'] = a['source_repository']
        row['Publication_DOI'] = e['Publication_DOI']
        row['Title'] = a['verified_title']
        row['First_Author'] = authors.get(dataset, '')
        row['Publication_Year'] = a['publication_year']
        row['Source_URL'] = a['source_url']
        row['Repository_Record_URL'] = a['repository_record_url']
        row['Processed_Input_Route'] = e['Processed_Input_Route']
        row['Processed_Input_URLs'] = e['Processed_Input_URLs']
        row['Processed_Record_URL'] = e['Processed_Record_URL']
        row['Data_License'] = e['Data_License']
        row['Derived_Matrix_Scope'] = e['Derived_Matrix_Scope']
        row['Article_License_URL'] = e['Article_License_URL']
    write_tsv(manifest_path, columns, rows)
    print({'rows': len(rows)})


if __name__ == '__main__':
    main()
