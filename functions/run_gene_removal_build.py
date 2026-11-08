"""Sequential, resumable real matrix construction; uses existing installed R and sparse chunks."""
import csv
import argparse
from collections import defaultdict, Counter
from decimal import Decimal
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import resource
import subprocess
import gzip
import itertools
import re
import io
import zipfile
from contextlib import ExitStack


def iter_rows(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def read_rows(path):
    return list(iter_rows(path))


def write_rows(path, rows, fields):
    partial = path.with_name(path.name+'.tmp')
    with partial.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def gtf_gene_ids(path, gene_hgnc_ids=None):
    result = defaultdict(set)
    with ExitStack() as stack:
        if path.suffix == '.zip':
            archive = stack.enter_context(zipfile.ZipFile(path))
            members = [name for name in archive.namelist() if name.endswith('.gtf')]
            if 'rsem_GRCh38.p2.gtf' in members:
                members = ['rsem_GRCh38.p2.gtf']
            if len(members) != 1:
                raise RuntimeError(f'Ambiguous gene-reference archive: {path}')
            handle = stack.enter_context(io.TextIOWrapper(archive.open(members[0])))
        else:
            handle = stack.enter_context(gzip.open(path, 'rt'))
        for line in handle:
            if line.startswith('#'):
                continue
            columns = line.rstrip('\n').split('\t')
            if len(columns) != 9 or columns[2] not in {'gene', 'transcript'}:
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', columns[8]))
            name = attributes.get('gene_name') or attributes.get('gene_symbol')
            identifier = attributes.get('gene_id')
            if name and identifier:
                identifier = identifier.split('.')[0]
                result[name].add(identifier)
                if gene_hgnc_ids is not None and attributes.get('hgnc_id'):
                    gene_hgnc_ids[identifier].add(attributes['hgnc_id'])
    return result


def gene_reference_context(reference_dir, panel):
    hgnc = {r['hgnc_id']: r for r in read_rows(reference_dir/'hgnc_complete_set.txt')
            if r['status'] == 'Approved'}
    approved = {r['symbol']: key for key, r in hgnc.items()}
    aliases, identifiers = defaultdict(set), defaultdict(set)
    for key, row in hgnc.items():
        for name in (row['prev_symbol']+'|'+row['alias_symbol']).split('|'):
            if name:
                aliases[name].add(key)
        for identifier in row['ensembl_gene_id'].split('|'):
            if identifier:
                identifiers[identifier].add(key)
    for row in panel:
        for identifier in filter(None, row['Ensembl_IDs'].split(';')):
            identifiers[identifier.split('.')[0]].add(row['Gene_Key'])
    entrez = defaultdict(set)
    for row in read_rows(reference_dir/'Homo_sapiens.gene_info.gz'):
        keys = {x[5:] for x in row['dbXrefs'].split('|') if x.startswith('HGNC:HGNC:')}
        keys &= hgnc.keys()
        if len(keys) == 1:
            entrez[row['GeneID']].update(keys)
            for value in row['dbXrefs'].split('|'):
                if value.startswith('Ensembl:ENSG'):
                    identifiers[value[8:]].update(keys)
    history = {}
    for row in iter_rows(reference_dir/'gene_history.gz'):
        if row['#tax_id'] == '9606':
            history[row['Discontinued_GeneID']] = row['GeneID']
    for old in history:
        current, visited = old, set()
        while current in history and current not in visited:
            visited.add(current)
            current = history[current]
        if current in entrez:
            entrez[old].update(entrez[current])
    historical_names = defaultdict(set)
    current_names = {}
    gencode_hgnc = defaultdict(set)
    for version in (17, 19, 24, 27, 32, 38, 43, 50):
        names = gtf_gene_ids(reference_dir/f'gencode.v{version}.annotation.gtf.gz',
                             gencode_hgnc if version == 50 else None)
        if version == 50:
            current_names = names
        for name, ids in names.items():
            for identifier in ids:
                historical_names[identifier].add(name)
        if version == 50:
            for identifier, keys in gencode_hgnc.items():
                supported = keys & hgnc.keys()
                identifiers[identifier].update(supported)
                if len(supported) == 1:
                    key = next(iter(supported))
                    if identifier in hgnc[key]['ensembl_gene_id'].split('|'):
                        identifiers[identifier] = {key}
            for ids in names.values():
                for identifier in ids:
                    if not identifiers[identifier]:
                        identifiers[identifier].add('ENSEMBL:'+identifier)
    official_names = defaultdict(set)
    for name, key in approved.items():
        official_names[name].add(key)
    for name, ids in current_names.items():
        for identifier in ids:
            official_names[name].update(identifiers[identifier])
    for identifier, names in historical_names.items():
        if not identifiers[identifier]:
            keys = set().union(*(official_names[name] for name in names))
            if len(keys) == 1:
                identifiers[identifier].update(keys)
    return dict(approved=approved, aliases=aliases, identifiers=identifiers,
                entrez=entrez, reference_dir=reference_dir, reported={})


def source_identity(row, dataset, context, metadata_ids):
    feature, key = row['Original_Feature_ID'], row['Gene_Key']
    supplied = set(re.findall(r'ENSG\d+', row.get('Source_Ensembl', '') or ''))
    supplied.update(re.findall(r'ENSG\d+', feature))
    supplied.update(metadata_ids.get(feature, set()))
    if dataset == 'AllenM1':
        filename = 'Allen_GRCh38p2_reference.zip'
        if dataset not in context['reported']:
            context['reported'][dataset] = gtf_gene_ids(context['reference_dir']/filename)
        original = context['reported'][dataset].get(feature, set())
        entrez_ids = set()
        for value in original:
            match = re.fullmatch(r'(?:GeneID[:_])?(\d+)', value)
            if match:
                entrez_ids.add(match[1])
        keys = set().union(*(context['entrez'][value] for value in entrez_ids)) if entrez_ids else set()
        complete = bool(entrez_ids) and all(context['entrez'][value] for value in entrez_ids)
        return complete and keys == {key}, {'NCBI:'+value for value in entrez_ids}
    if not supplied and dataset in {'SomaMut', 'HYPOMAP'}:
        filename = 'gencode.v32.annotation.gtf.gz' if dataset == 'SomaMut' else 'Homo_sapiens.GRCh38.98.gtf.gz'
        if dataset not in context['reported']:
            context['reported'][dataset] = gtf_gene_ids(context['reference_dir']/filename)
        supplied.update(context['reported'][dataset].get(feature, set()))
    keys = set().union(*(context['identifiers'][value] for value in supplied)) if supplied else set()
    complete = bool(supplied) and all(context['identifiers'][value] for value in supplied)
    return complete and keys == {key}, supplied


def count_feature_ids(path):
    identifiers = defaultdict(set)
    with path.open() as handle:
        next(handle)
        for line in handle:
            feature = line.split('\t', 1)[0].strip('"\r\n')
            identifier, separator, symbol = feature.partition('|')
            if separator:
                identifiers[symbol.replace('_', '-')].update(re.findall(r'ENSG\d+', identifier))
    if not any(identifiers.values()):
        raise RuntimeError(f'No original Ensembl count-row identifiers found: {path}')
    return identifiers


def prepare_source_policy(run, repo, dataset, context):
    output = run/'run/source_gene_policy'/dataset
    marker = output/'applied_gene_removal_policy.json'
    identity_rule = 'stable_reference_ids_or_historical_official_names'
    if marker.exists():
        return output
    observations = run/'run/source_gene_groups'/dataset
    if not (observations/'SUCCESS.json').is_file():
        subprocess.run([os.environ.get('BRAINOMICS_RSCRIPT', 'Rscript'),
                        str(repo/'functions/build_clean_source_counts.R'),
                        '--describe-source', str(run), dataset], cwd=repo, check=True)
    receipt = json.loads((observations/'SUCCESS.json').read_text())
    if receipt['dataset'] != dataset or not receipt['all_cells_visited_once']:
        raise RuntimeError(f'Incomplete source-gene observations: {dataset}')
    mapping = [r for r in read_rows(run/'code/source_row_audit_map.tsv.gz') if r['Dataset'] == dataset]
    native = read_rows(observations/'native_rows.tsv.gz')
    groups = read_rows(observations/'gene_groups.tsv.gz')
    metadata_ids = defaultdict(set)
    if dataset == 'GSE81475':
        data_root = Path(os.environ.get('BRAINOMICS_DATA_ROOT', str(repo.parents[1]/'data/BrainOmicsData')))
        count_file = data_root/'raw/GSE81475/counts.txt'
        metadata_ids.update(count_feature_ids(count_file))
    feature_file = run/'inputs/processed'/dataset/'features_raw.tsv.gz'
    if feature_file.is_file():
        for row in read_rows(feature_file):
            name = row.get('Original_Feature_ID', '')
            metadata_ids[name].update(re.findall(r'ENSG\d+', ' '.join(row.values())))
    mask = defaultdict(set)
    gene_keys = {(r['Dataset'], r['Original_Feature_ID']): r['Gene_Key'] for r in mapping}
    for row in mapping:
        feature = row['Original_Feature_ID']
        if row['Action'] == 'quarantine':
            mask['*', feature].add('previously_quarantined_source_identity')
        if row['Action'] != 'candidate' or row['Source_Ensembl'] or re.search(r'ENSG\d+', feature):
            continue
        name = row['Existing_Feature']
        if name not in context['approved'] and len(context['aliases'][name]) > 1:
            supported, identifiers = source_identity(row, dataset, context, metadata_ids)
            if not supported:
                mask['*', feature].add('unresolved_source_alias_identity')
    by_group = defaultdict(list)
    multiple = {(r['Layer'], r['Gene_Key']): r for r in groups}
    identity_ids = {}
    for row in native:
        pair = row['Layer'], row['Gene_Key']
        if pair not in multiple:
            continue
        by_group[pair].append(row)
        supported, identifiers = source_identity(row, dataset, context, metadata_ids)
        identity_ids[row['Layer'], row['Original_Feature_ID']] = identifiers
        if not supported:
            mask[row['Layer'], row['Original_Feature_ID']].add('unresolved_source_row_identity')
    umi_sources = {'AllenM1', 'GSE168408', 'GSE217511', 'GSE296073', 'SomaMut'}
    allowed = []
    for pair, group in multiple.items():
        layer, key = pair
        remaining = [r for r in by_group[pair]
                     if ('*', r['Original_Feature_ID']) not in mask
                     and (layer, r['Original_Feature_ID']) not in mask]
        ids = [identity_ids[layer, r['Original_Feature_ID']] for r in remaining]
        disjoint = all(not (a & b) for a, b in itertools.combinations(ids, 2))
        all_supported = len(remaining) == len(by_group[pair]) and all(ids)
        if dataset in umi_sources and all_supported and disjoint and int(group['Positive_Rows']) > 1:
            allowed.append(dict(Dataset=dataset, Layer=layer, Gene_Key=key))
        elif int(group['Coexpressed_Cells']) and sum(int(r['Detected_Cells']) > 0 for r in remaining) > 1:
            for row in remaining:
                mask[layer, row['Original_Feature_ID']].add('unresolved_additivity_in_coexpressed_source_group')
    records = [dict(Dataset=dataset, Layer=layer, Original_Feature_ID=feature,
                    Gene_Key=gene_keys[dataset, feature], Reason=';'.join(sorted(reasons)))
               for (layer, feature), reasons in sorted(mask.items())]
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output/'applied_gene_row_exclusions.tsv', records,
               ['Dataset', 'Layer', 'Original_Feature_ID', 'Gene_Key', 'Reason'])
    write_rows(output/'allowed_count_addition.tsv', allowed, ['Dataset', 'Layer', 'Gene_Key'])
    policy = dict(dataset=dataset,
                  retained_cells=receipt['retained_cells'], cell_membership_modified=False,
                  feature_to_gene_assignments_modified=False,
                  exclusions=len(records), allowed_groups=len(allowed),
                  reason_counts=dict(Counter(reason for reasons in mask.values() for reason in reasons)),
                  identity_rule=identity_rule)
    partial = marker.with_name(marker.name+'.tmp')
    partial.write_text(json.dumps(policy, indent=2)+'\n')
    partial.replace(marker)
    return output

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset')
    parser.add_argument('--prepare-policy', action='store_true')
    args = parser.parse_args()
    run_value = os.environ.get('BRAINOMICS_RUN_ROOT')
    assert run_value, 'BRAINOMICS_RUN_ROOT must name the run directory holding the source inputs and gene-panel policy'
    run = Path(run_value).resolve()
    out = run / 'run/gene_removal_build'
    out.mkdir(parents=True, exist_ok=True)
    lock = (out/'runner.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    repo = Path(__file__).resolve().parents[1]
    script = repo/'functions/build_clean_source_counts.R'
    assert script.is_file(), f'Assembly source is missing: {script}'
    list_file = Path(os.environ.get(
        'BRAINOMICS_DATASET_LIST', repo/'data/source_redistribution_evidence.tsv'))
    assert list_file.is_file(), f'dataset list not found: {list_file}'
    lines = list_file.read_text().splitlines()
    datasets = lines if list_file.suffix == '.txt' else [line.split('	')[0] for line in lines[1:] if line]
    assert len(datasets) == len(set(datasets)) == 22, datasets
    if args.dataset:
        if args.dataset not in datasets:
            raise RuntimeError(f'Unknown reference dataset: {args.dataset}')
        datasets = [args.dataset]
    env = dict(os.environ, BRAINOMICS_REPO_ROOT=str(repo), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    completed = []
    state = dict(pid=os.getpid(), counts_modified=not args.prepare_policy, cell_membership_modified=False,
                 sources=len(datasets), memory_limit_gib=96, sequential_sources=True)


    def status(stage, **extra):
        state.update(state=stage, updated_utc=datetime.now(timezone.utc).isoformat(),
                     completed_sources=len(completed), **extra)
        temp = out/'status.json.tmp'
        temp.write_text(json.dumps(state, indent=2)+'\n'); temp.replace(out/'status.json')


    try:
        reference_context = None
        legacy_policy = all((run/name).is_file() for name in
                            ['applied_gene_row_exclusions.tsv', 'allowed_count_addition.tsv',
                             'applied_gene_removal_policy.json'])
        resource.setrlimit(resource.RLIMIT_AS, (96 * 1024**3, 96 * 1024**3))
        os.nice(19)
        if not legacy_policy:
            status('PREPARING_GENE_REFERENCES', current_source=None)
            data_root = Path(os.environ.get('BRAINOMICS_DATA_ROOT', str(repo.parents[1]/'data/BrainOmicsData')))
            reference_dir = Path(os.environ.get('BRAINOMICS_GENE_REFERENCE_DIR', str(data_root/'raw/gene_references')))
            panel = read_rows(run/'inputs/plan/selected_gene_panel_34387.tsv')
            reference_context = gene_reference_context(reference_dir, panel)
        for ds in datasets:
            destination = out/ds; destination.mkdir(exist_ok=True)
            marker = destination/'SUCCESS.json'
            status('BUILDING_CLEANED_COUNTS', current_source=ds)
            if reference_context is not None:
                prepare_source_policy(run, repo, ds, reference_context)
            if args.prepare_policy:
                continue
            if not marker.exists():
                with (destination/'driver.log').open('a') as log:
                    subprocess.run([os.environ.get('BRAINOMICS_RSCRIPT', 'Rscript'), str(script), str(run), ds],
                                   cwd=repo, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            result = json.loads(marker.read_text())
            assert result['dataset'] == ds
            assert result['all_cells_visited_once'] and result['counts_modified']
            assert result['no_unsupported_positive_addition'] and not result['cell_membership_modified']
            completed.append(result)
        if args.prepare_policy:
            status('SOURCE_GENE_POLICIES_COMPLETE', current_source=None,
                   policy_sources=datasets)
            return
        if args.dataset:
            status('SELECTED_SOURCE_COUNTS_COMPLETE', current_source=None,
                   reference_cells=sum(r['reference_cells'] for r in completed))
            return
        assert sum(r['reference_cells'] for r in completed) == 2602031
        status('SELECTING_FINAL_GENE_PANEL', current_source=None)
        maximum, total = {}, {}
        for ds in datasets:
            with (out/ds/'gene_detection.tsv').open() as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
            assert len(rows) == len({r['Gene_Key'] for r in rows}) == 34387
            for r in rows:
                value = Decimal(r['Detected_Cells'])
                assert value.is_finite() and value >= 0 and value == value.to_integral_value()
                key, n = r['Gene_Key'], int(value)
                maximum[key] = max(maximum.get(key, 0), n)
                total[key] = total.get(key, 0) + n
        with (run/'inputs/plan/selected_gene_panel_34387.tsv').open() as f:
            candidates = list(csv.DictReader(f, delimiter='\t'))
        fields = list(candidates[0]) + ['Cleaned_Max_Source_Detection', 'Cleaned_Total_Detection']
        retained, excluded = [], []
        for row in candidates:
            key = row['Gene_Key']
            row.update(Cleaned_Max_Source_Detection=maximum[key], Cleaned_Total_Detection=total[key])
            (retained if maximum[key] >= 10 else excluded).append(row)
        for name, rows in [('final_gene_panel.tsv', retained), ('genes_removed_after_count_cleaning.tsv', excluded)]:
            with (out/name).open('w') as f:
                w = csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(rows)
        summary = dict(reference_cells=2602031, sources=22, candidate_genes=34387,
                       final_genes=len(retained), removed_genes=len(excluded),
                       zero_count_cells_before_final_panel=sum(r['zero_count_cells'] for r in completed),
                       source_results=completed,
                       matrix_stage='Source count chunks complete; apply final panel to assembly and verify zero-count cells before integration.')
        (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        status('CLEANED_SOURCE_COUNTS_COMPLETE', current_source=None, final_genes=len(retained))
    except Exception as error:
        status('FAILED', error=str(error))
        raise


if __name__ == '__main__':
    main()
