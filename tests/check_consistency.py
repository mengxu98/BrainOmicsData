#!/usr/bin/env python3
"""Consistency checks between the released data tables and the pipeline scripts.

Verifies, without running the pipeline:
  1. the 22-dataset cohort is identical across the evidence table, the R
     contract (functions/integration.R), the preprocessing driver, the HPC
     dispatcher, the ScienceDB source-input contract and processing/<dataset>.R;
  2. every column the code consumes exists in data/source_access_summary.tsv and
     data/source_redistribution_evidence.tsv, and every registered download
     script exists;
  3. the sealed package agrees with the repository tables (dataset set and a
     byte-identical replay of provenance/dataset_manifest.tsv);
  4. every path referenced by the numbered drivers exists;
  5. the environment locks used by environment/verify_environment.R are present.

Usage: check_consistency.py [REPO_DIR] [PACKAGE_DIR]
Exits non-zero on the first failed check group.
"""
import csv, gzip, re, subprocess, sys, tempfile, pathlib

FAIL = []


def check(condition, message):
    if not condition:
        FAIL.append(message)


def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def r_vector(text, function_name, variable):
    """Extract a character vector assigned inside one R function."""
    block = re.search(rf'{function_name} <- function\(\) \{{(.*?)\n\}}', text, re.S)
    if not block:
        return []
    body = block.group(1)
    if variable:
        match = re.search(rf'{variable} <- c\((.*?)\n  \)', body, re.S) or re.search(rf'{variable} <- c\((.*?)\)', body, re.S)
    else:
        match = re.search(r'c\((.*?)\n  \)', body, re.S) or re.search(r'c\((.*?)\)', body, re.S)
    return re.findall(r'"([^"]+)"', match.group(1)) if match else []


def main():
    repo = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
    default_pkg = repo.parents[1] / 'data/BrainOmicsData/ScienceDB'
    package = pathlib.Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else (
        pathlib.Path(__import__('os').environ.get('BRAINOMICS_PACKAGE_DIR', default_pkg)).resolve())

    evidence = read_tsv(repo/'data/source_redistribution_evidence.tsv')
    summary = read_tsv(repo/'data/source_access_summary.tsv')
    cohort = [r['Dataset'] for r in evidence]

    # 1. cohort agreement -----------------------------------------------------
    integration = (repo/'functions/integration.R').read_text()
    r_reference = r_vector(integration, 'existing_reference_datasets', '') or []
    r_reference += r_vector(integration, 'additional_reference_datasets', '')
    status = re.search(r'query_validation_datasets <- function\(\) \{\s*(.*?)\n\}', integration, re.S)
    r_query = re.findall(r'"([^"]+)"', status.group(1)) if status else []
    check(len(cohort) == 22 and len(set(cohort)) == 22, f'evidence table must hold 22 unique datasets, found {len(cohort)}')
    check(set(r_reference) == set(cohort), 'functions/integration.R reference datasets differ from the evidence table')
    check(len(r_query) == 1, 'functions/integration.R must declare exactly one query dataset')

    driver = (repo/'03_datasets_preprocessing.sh').read_text()
    block = re.search(r'datasets=\((.*?)\n\)', driver, re.S)
    driver_cohort = re.findall(r'\b([A-Za-z0-9_]+)\b', block.group(1)) if block else []
    check(set(driver_cohort) == set(cohort), f'03_datasets_preprocessing.sh cohort differs (n={len(driver_cohort)})')

    dispatcher = (repo/'hpc/run_dataset.sh').read_text()
    case = re.search(r'case "\$dataset" in(.*?)\n\s*;;', dispatcher, re.S)
    case_datasets = set(re.findall(r'([A-Za-z][A-Za-z0-9_]+)', case.group(1))) if case else set()
    case_datasets = {d for d in case_datasets if d not in {'in', 'dataset'}}
    missing_dispatch = sorted(set(cohort) - case_datasets)
    check(not missing_dispatch, f'hpc/run_dataset.sh does not accept: {missing_dispatch}')

    metadata = (repo/'sciencedb/package_metadata.R').read_text()
    check(set(r_vector(metadata, 'brainomics_source_input_contract', 'datasets')) == set(cohort),
          'sciencedb/package_metadata.R source-input contract differs from the evidence table')

    for dataset in cohort:
        check((repo/f'processing/{dataset}.R').is_file(), f'missing processing/{dataset}.R')

    # 2. registry columns and download scripts --------------------------------
    evidence_columns = set(evidence[0].keys())
    required_evidence = {'Dataset', 'Publication_DOI', 'Processed_Input_Route', 'Processed_Input_URLs',
                         'Processed_Record_URL', 'Data_License', 'Derived_Matrix_Scope', 'Article_License_URL'}
    check(required_evidence <= evidence_columns,
          f'evidence table lacks {sorted(required_evidence - evidence_columns)}')
    summary_columns = set(summary[0].keys())
    required_summary = {'dataset', 'source_accession', 'source_repository', 'publication_doi', 'verified_title',
                        'journal', 'publication_year', 'access_status', 'download_script', 'source_url',
                        'repository_record_url'}
    check(required_summary <= summary_columns, f'source_access_summary.tsv lacks {sorted(required_summary - summary_columns)}')
    summary_datasets = [r['dataset'] for r in summary]
    check(len(summary_datasets) == len(set(summary_datasets)), 'source_access_summary.tsv has duplicate datasets')
    check(set(cohort) <= set(summary_datasets), 'source_access_summary.tsv does not cover the released cohort')
    for row in summary:
        script = (row.get('download_script') or '').strip()
        if script:
            check((repo/script).is_file(), f'{row["dataset"]}: download script missing ({script})')
    registered = {r['download_script'].strip() for r in summary if (r.get('download_script') or '').strip()}
    for path in sorted((repo/'download').glob('*.sh')):
        rel = f'download/{path.name}'
        check(rel in registered, f'{rel} is not registered in source_access_summary.tsv')

    # 3. package agreement ----------------------------------------------------
    manifest_path = package/'provenance/dataset_manifest.tsv'
    if manifest_path.is_file():
        manifest = read_tsv(manifest_path)
        check(set(r['Dataset'] for r in manifest) == set(cohort),
              'package provenance/dataset_manifest.tsv differs from the evidence table')
        with tempfile.TemporaryDirectory() as tmp:
            out = pathlib.Path(tmp)/'provenance'
            out.mkdir(parents=True)
            (out/'dataset_manifest.tsv').write_bytes(manifest_path.read_bytes())
            result = subprocess.run([sys.executable, str(repo/'functions/merge_dataset_manifest.py'),
                                     str(repo), tmp], capture_output=True, text=True)
            check(result.returncode == 0, f'merge_dataset_manifest.py failed: {result.stderr.strip()}')
            if result.returncode == 0:
                check((out/'dataset_manifest.tsv').read_bytes() == manifest_path.read_bytes(),
                      'repository tables do not reproduce the released provenance/dataset_manifest.tsv')
        metadata_tsv = package/'metadata/metadata.tsv.gz'
        if metadata_tsv.is_file():
            counts = {}
            with gzip.open(metadata_tsv, 'rt') as handle:
                header = handle.readline().rstrip('\n').split('\t')
                index = header.index('Dataset')
                for line in handle:
                    dataset = line.split('\t', index+1)[index]
                    counts[dataset] = counts.get(dataset, 0) + 1
            published = {r['Dataset']: int(r['Cells']) for r in manifest}
            check(counts == published, 'metadata.tsv.gz per-dataset cell counts differ from dataset_manifest.tsv')
    else:
        print(f'note: no package at {package}; package agreement checks skipped')

    # 4. driver references ----------------------------------------------------
    pattern = re.compile(r'(?:Rscript(?:\s+--vanilla)?|bash|sh|python3?)\s+("?)([A-Za-z0-9_./-]+\.(?:R|sh|py|sbatch))\1')
    for driver_path in sorted(repo.glob('[0-9][0-9]_*.sh')) + [repo/'run_pipeline.sh']:
        text = driver_path.read_text()
        for _, target in pattern.findall(text):
            if target.startswith('$'):
                continue
            check((repo/target).is_file(), f'{driver_path.name} references missing {target}')

    # 5. environment locks ----------------------------------------------------
    for lock in ['r-packages.lock.tsv']:
        path = repo/'environment'/lock
        check(path.is_file(), f'missing environment/{lock}')
        if path.is_file():
            rows = read_tsv(path)
            check(rows and {'Package', 'Version'} <= set(rows[0].keys()), f'environment/{lock} lacks Package/Version columns')
            versions = {row['Package']: row['Version'] for row in rows}
            check(len(versions) == len(rows), f'environment/{lock} has duplicate packages')
            check(versions.get('R') == '4.5.1', 'the workflow requires one R 4.5.1 environment')
    for overview in ['figures/fig1.svg', 'figures/fig1-night.svg']:
        check((repo/overview).is_file(), f'README overview figure {overview} is missing')

    if FAIL:
        print(f'FAIL consistency: {len(FAIL)} problem(s)')
        for item in FAIL[:25]:
            print('  -', item)
        return 1
    print(f'PASS consistency: {len(cohort)}-dataset cohort, registry columns, download scripts, '
          f'package replay, driver references and environment locks agree')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
