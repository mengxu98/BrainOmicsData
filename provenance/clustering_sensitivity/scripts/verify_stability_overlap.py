"""Independently recompute full-cohort partition metrics from saved overlap counts."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import csv, hashlib, json, math, sys

def read_table(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def verify(directory):
    contract = json.loads((directory / 'contract.json').read_text())
    complete = json.loads((directory / 'COMPLETE.json').read_text())
    summary, = read_table(directory / 'summary.tsv')
    overlaps = read_table(directory / 'overlap.tsv')
    cluster_rows = read_table(directory / 'cluster_stability.tsv')
    expected = {row['Cluster']: int(row['Cells']) for row in
                read_table(Path(__file__).resolve().parents[1] / 'annotations/final_cluster_annotation.tsv')}
    baseline, perturbed, pairs = Counter(), Counter(), {}
    for row in overlaps:
        pair = row['Baseline'], row['Perturbed']
        count = int(row['Freq'])
        assert count > 0 and pair not in pairs
        pairs[pair] = count
        baseline[pair[0]] += count
        perturbed[pair[1]] += count
    assert dict(baseline) == expected
    n = sum(baseline.values())
    assert n == 2602031 == int(summary['Cells']) == contract['cells'] == complete['cells']
    assert int(summary['Clusters']) == len(perturbed)
    assert contract['task'] == complete['task'] == int(summary['Task'])
    assert complete['state'] == 'COMPLETE'
    assert contract['subsampling'] is False and contract['formal_annotation_modified'] is False
    setting, = contract['settings']
    assert setting['seed'] == int(summary['Seed'])
    assert setting['resolution'] == float(summary['Resolution'])
    choose2 = lambda value: value * (value - 1) // 2
    a = sum(map(choose2, baseline.values()))
    b = sum(map(choose2, perturbed.values()))
    observed = sum(map(choose2, pairs.values()))
    expected_pairs = Fraction(a * b, choose2(n))
    ari = float((observed - expected_pairs) / (Fraction(a + b, 2) - expected_pairs))
    assert math.isclose(ari, float(summary['ARI']), abs_tol=1e-12)
    assert math.isclose(ari, complete['ARI'], abs_tol=1e-4)  # jsonlite default precision
    assert len(cluster_rows) == 75 and {r['Baseline'] for r in cluster_rows} == set(baseline)
    for row in cluster_rows:
        key = row['Baseline']
        fragments = {b: count for (a, b), count in pairs.items() if a == key}
        assert int(row['Cells']) == baseline[key] and int(row['Fragments']) == len(fragments)
        largest = max(fragments.values()) / baseline[key]
        jaccard = max(count / (baseline[key] + perturbed[b] - count) for b, count in fragments.items())
        assert math.isclose(largest, float(row['Largest_Fragment_Fraction']), abs_tol=1e-12)
        assert math.isclose(jaccard, float(row['Best_Jaccard']), abs_tol=1e-12)
    return {'state': 'PASS_FULL_OVERLAP_RECOMPUTATION', 'task': contract['task'],
            'cells': n, 'clusters': len(perturbed), 'ARI': ari,
            'graph_md5': contract['graph_md5'], 'baseline_md5': contract['baseline_md5'],
            'input_sha256': {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                             for name in ['contract.json', 'COMPLETE.json', 'summary.tsv',
                                          'overlap.tsv', 'cluster_stability.tsv']},
            'scope': 'Independent exact-integer overlap arithmetic on saved full-cohort counts. Does not independently reload per-cell assignments, rerun clustering, or validate biological labels.'}

if __name__ == '__main__':
    directory = Path(sys.argv[1])
    report = verify(directory)
    (directory / 'independent_overlap_check.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
