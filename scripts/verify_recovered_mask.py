#!/usr/bin/env python3
"""Replay archived SD06 observation-mask counts without BAMs or production imports.

The inputs are byte-pinned by the original SD07 contract. This checks retained
intervals and their CDS intersection. It does not regenerate masks from reads,
verify reference sequence bases, or modify the historical SD13 pending status.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

CONTRACT_SHA256 = '9a6de46c2901b216fe7ecdc4d4d530b5b6ea21df5b7b5869bd865f9478c58184'
PREFIX = 'MQ40_BQ20.non_sparse_142.DP5.GE128.'
GENOME = 'resource/masks/' + PREFIX + 'REFERENCE_ACGT.bed.gz'
CDS = 'resource/masks/' + PREFIX + 'CDS_REFERENCE_ACGT.bed.gz'
UNION = 'resource/tables/CDS_union_intervals.tsv'
SUMMARY = 'resource/tables/mask_summary.tsv'
DEFINITION = 'resource/resource_definition.json'
REQUIRED = [GENOME, CDS, UNION, SUMMARY, DEFINITION]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        require(0 <= start < end, 'Invalid interval')
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def intersect(left, right):
    result = {}
    for chrom in left:
        a, b = left[chrom], right.get(chrom, [])
        i = j = 0
        rows = []
        while i < len(a) and j < len(b):
            start = max(a[i][0], b[j][0])
            end = min(a[i][1], b[j][1])
            if start < end:
                rows.append((start, end))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
        if rows:
            result[chrom] = merge(rows)
    return result


def read_bed(path):
    result = {}
    previous_chrom = None
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        for number, line in enumerate(f, 1):
            fields = line.rstrip('\n').split('\t')
            require(len(fields) == 3, 'BED3 expected at line %s' % number)
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            require(0 <= start < end, 'Invalid BED interval')
            if chrom != previous_chrom:
                require(chrom not in result, 'BED contig block repeated')
                result[chrom] = []
                previous_chrom = chrom
            rows = result[chrom]
            require(not rows or start >= rows[-1][1], 'BED overlap or order error')
            rows.append((start, end))
    require(bool(result), 'Empty BED')
    return result


def interval_counts(intervals):
    return {'retained_bp': sum(b - a for rows in intervals.values() for a, b in rows),
            'bed_intervals': sum(len(rows) for rows in intervals.values())}


def verify_recovered_mask(data, out):
    """Verify six pinned input files; write fresh JSON/TSV evidence; return summary."""
    data, out = Path(data), Path(out)
    require(not out.exists(), 'Output directory already exists; choose a new directory')
    contract_file = data / 'provenance/SD07_contract.json'
    require(sha256(contract_file) == CONTRACT_SHA256, 'SD07 contract identity mismatch')
    contract = json.loads(contract_file.read_text(encoding='utf-8'))
    identities = []
    for rel in REQUIRED:
        source = contract['source_files'][rel]
        p = data / rel
        actual = sha256(p)
        require(p.stat().st_size == source['bytes'], 'Size mismatch: ' + rel)
        require(actual == source['sha256'], 'SHA256 mismatch: ' + rel)
        identities.append({'path': rel, 'bytes': p.stat().st_size, 'sha256': actual,
                           'original_path': source['path'], 'status': 'PASS'})
    definition = json.loads((data / DEFINITION).read_text(encoding='utf-8'))
    require(definition['depth'] == 5 and definition['samples'] == 142,
            'Mask definition mismatch')
    genome, cds = read_bed(data / GENOME), read_bed(data / CDS)
    union_raw, union_rows = {}, 0
    with (data / UNION).open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        require(reader.fieldnames == ['gene', 'chrom', 'start0', 'end0'], 'CDS schema mismatch')
        for row in reader:
            union_raw.setdefault(row['chrom'], []).append((int(row['start0']), int(row['end0'])))
            union_rows += 1
    union = {chrom: merge(rows) for chrom, rows in union_raw.items()}
    overlap = intersect(genome, union)
    require(overlap == cds,
            'Independent genome/CDS intersection differs from archived CDS mask')
    require(union_rows == 42025, 'CDS row count mismatch')
    union_bp = interval_counts(union)['retained_bp']
    require(union_bp == contract['CDS_reference_acgt_bp'], 'CDS union denominator mismatch')
    with (data / SUMMARY).open(newline='', encoding='utf-8') as f:
        source_rows = list(csv.DictReader(f, delimiter='\t'))
    results = []
    for scope, bed_name, intervals, denominator in [
        ('REFERENCE_ACGT', GENOME, genome, contract['reference_acgt_bp']),
        ('CDS_REFERENCE_ACGT', CDS, cds, contract['CDS_reference_acgt_bp'])
    ]:
        matches = [r for r in source_rows if r['bed_file'] == bed_name.removeprefix('resource/')]
        require(len(matches) == 1, 'Missing or repeated mask summary row')
        row = matches[0]
        require((row['profile'], row['source_group'], row['minimum_depth'],
                 row['required_sample_count'], row['samples'], row['scope']) ==
                ('MQ40_BQ20', 'non_sparse_142', '5', '128', '142', scope), 'Summary scope mismatch')
        counts = interval_counts(intervals)
        require(int(row['denominator_bp']) == denominator, 'Denominator record mismatch')
        for key, value in counts.items():
            require(int(row[key]) == value, key + ' differs from original summary')
        fraction = counts['retained_bp'] / denominator
        require(abs(float(row['retained_fraction']) - fraction) < 1e-14, 'Fraction mismatch')
        results.append({'scope': scope, **counts, 'denominator_bp': denominator,
                        'retained_fraction': fraction, 'status': 'PASS'})
    summary = {
        'status': 'ARCHIVED_MASK_REPLAY_PASS',
        'source_run': 'RAMOS_SD06_20261120T184153Z_c5e96b11',
        'source_contract': 'RAMOS_SD07_CERES_VARIANT_OBSERVATION_AUDIT_v0.7.0',
        'source_contract_sha256': CONTRACT_SHA256,
        'profile': 'MQ40_BQ20', 'minimum_depth': 5,
        'required_sample_count': 128, 'samples': 142,
        'counts': results,
        'cds_union': {'input_rows': union_rows, **interval_counts(union)},
        'independent_intersection': {**interval_counts(overlap), 'exact_interval_equality': True},
        'historical_SD13_status_unchanged': 'MASK_AUDIT_PENDING',
        'limits': [
            'Independent interval replay of archived, hash-pinned masks; no BAM processing.',
            'Reference and CDS ACGT labels/denominators inherit original SD06/SD07 provenance. FASTA and GFF contents were not reread.',
            'The unique CDS union length independently matches the recorded CDS ACGT denominator.',
            'At least 128 of the fixed 142 samples support each position; supporting samples can differ by position.',
            'Observation is not evidence of genotype correctness, independent molecules, mapping accuracy, or gene presence/absence.',
            'This replay does not modify the original SD13 record and does not check targeted-site membership.'
        ]
    }
    out.mkdir(parents=True)
    for filename, obj in [('mask_replay_summary.json', summary), ('mask_input_identities.json', identities)]:
        (out / filename).write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')
    with (out / 'mask_replayed_counts.tsv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(results)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    summary = verify_recovered_mask(args.data, args.out)
    print(summary['status'])
    for row in summary['counts']:
        print('%s: %s/%s bp' % (row['scope'], row['retained_bp'], row['denominator_bp']))


if __name__ == '__main__':
    main()
