#!/usr/bin/env python3
"""Draw manuscript figures from frozen resource tables and carrier-group recounts."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

BLUE = '#0072B2'
ORANGE = '#D55E00'
GRAY = '#777777'
EIGHT = {'FSP348', 'FSP350', 'FSP351', 'FSP353', 'FSP357', 'FSP360', 'FSP363', 'FSP369'}


def draw(tables, out):
    if out.exists():
        raise FileExistsError('Choose a new figure output directory: ' + str(out))
    out.mkdir(parents=True)
    def rows(name):
        with (tables / name).open(encoding='utf-8', newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    })
    obs = rows('observation_by_isolate.tsv')
    calls = rows('primary142_isolate_call_counts.tsv')
    summary = json.loads((tables / 'structure_recount.json').read_text(encoding='utf-8'))
    if set(summary['eight_ids']) != EIGHT:
        raise ValueError('Eight-isolate membership differs from the documented grouping')
    if len(obs) != 169 or len(calls) != 142:
        raise ValueError('This figure expects 169 observation records and 142 genotype records')
    if set(x['fsp'] for x in calls if x['descriptive_group'].startswith('eight')) != EIGHT:
        raise ValueError('Isolate summary grouping differs from the documented grouping')
    fig, axes = plt.subplots(2, 2, figsize=(10.3, 7.55), layout='constrained')
    ax = axes.flat
    groups = [
        ([r for r in obs if r['primary'] == '0'], 'Sparse alignments (27)', GRAY, 'x'),
        ([r for r in obs if r['primary'] == '1' and r['fsp'] not in EIGHT], 'Other primary isolates (134)', BLUE, 'o'),
        ([r for r in obs if r['fsp'] in EIGHT], 'High-ALT primary isolates (8)', ORANGE, '^'),
    ]
    for group, label, color, marker in groups:
        ax[0].scatter([float(r['mean_aligned_depth']) for r in group],
                      [100 * float(r['fraction_reference_DP5']) for r in group],
                      s=27, alpha=.75, color=color, marker=marker, label=label,
                      linewidths=.8 if marker == 'x' else 0)
    ax[0].set(xscale='log', xlabel='Mean aligned-base depth',
              ylabel='Reference A/C/G/T positions\nwith depth ≥5 (%)', ylim=(-3, 103))
    ax[0].set_title('a', loc='left', fontweight='bold', fontsize=13)
    ax[0].legend(loc='center right', frameon=False, fontsize=8.1)
    for group, label, color, marker in groups[1:]:
        ax[1].scatter([100 * float(r['fraction_reference_DP5']) for r in group],
                      [100 * float(r['fraction_CDS_DP5']) for r in group],
                      s=29, alpha=.75, color=color, marker=marker, label=label, linewidths=0)
    ax[1].set(xlabel='Reference A/C/G/T positions with depth ≥5 (%)',
              ylabel='Coding-sequence A/C/G/T positions\nwith depth ≥5 (%)',
              xlim=(71, 94), ylim=(83, 99))
    ax[1].set_title('b', loc='left', fontweight='bold', fontsize=13)
    ax[1].legend(loc='lower right', frameon=False, fontsize=8.1)
    ordered = sorted(calls, key=lambda r: (int(r['n_ALT']), int(r['fsp'][3:])))
    for ids, label, color, marker in [(False, 'Other primary isolates (134)', BLUE, 'o'),
                                     (True, 'High-ALT primary isolates (8)', ORANGE, '^')]:
        values = [(i + 1, int(r['n_ALT'])) for i, r in enumerate(ordered) if (r['fsp'] in EIGHT) == ids]
        ax[2].scatter([r[0] for r in values], [r[1] for r in values],
                      color=color, marker=marker, s=22, linewidths=0, label=label, alpha=.8)
    ax[2].set(yscale='log', xlabel='Isolate rank by number of ALT calls',
              ylabel='Retained ALT calls per isolate', xlim=(-2, 146), ylim=(2, 50000))
    ax[2].set_title('c', loc='left', fontweight='bold', fontsize=13)
    ax[2].set_xticks([1, 25, 50, 75, 100, 125, 142])
    ax[2].set_yticks([10, 100, 1000, 10000])
    ax[2].yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{int(v):,}'))
    ax[2].legend(loc='upper left', frameon=False, fontsize=8.1)
    bin_keys = ['1', '2', '3_5', '6_14', '15plus']
    bin_labels = ['1', '2', '3–5', '6–14', '≥15']
    values = [int(summary['MAC_bins'][key]) for key in bin_keys]
    if sum(values) != int(summary['counts']['selected_sites']):
        raise ValueError('Minor-allele counts do not sum to selected-panel sites')
    bars = ax[3].bar(bin_labels, values, color=BLUE, width=.65)
    ax[3].bar_label(bars, labels=[f'{v:,}' for v in values], padding=4, fontsize=9)
    ax[3].set(xlabel='Minor allele count among retained calls',
              ylabel='Selected-panel sites', ylim=(0, max(values) * 1.20))
    ax[3].set_title('d', loc='left', fontweight='bold', fontsize=13)
    ax[3].yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{int(v):,}'))
    for a in axes.flat:
        a.grid(axis='y', color='#DDDDDD', lw=.6)
        a.set_axisbelow(True)
    fig.savefig(out / 'Figure_1_observation.png', dpi=300)
    fig.savefig(out / 'Figure_1_observation.pdf')
    plt.close(fig)

    support = rows('read_support_by_allele.tsv')
    parts = ['SELECTED', 'OUTSIDE_MASK_ONLY', 'LOW_CALL_COUNT_ONLY', 'BOTH']
    labels = ['Selected panel', 'Fails observation mask only',
              'Fails retained-call count only', 'Fails both criteria']
    site_counts = [128, 85, 128, 128]
    profiles = ['MQ20_BQ20_noBAQ', 'MQ40_BQ20_noBAQ', 'MQ40_BQ20_redoBAQ']
    fig, axes = plt.subplots(2, 2, figsize=(10.3, 7.1), layout='constrained')
    for i, (part, label, n, a) in enumerate(zip(parts, labels, site_counts, axes.flat)):
        for allele, color, marker, linestyle, y_text in [
            ('REF', BLUE, 'o', '-', .14), ('ALT', ORANGE, 's', '--', .065)
        ]:
            rows_for_allele = [next(v for v in support if v['partition'] == part and v['profile'] == pr and v['baseline_allele'] == allele) for pr in profiles]
            denoms = {int(v['baseline_retained']) for v in rows_for_allele}
            if len(denoms) != 1:
                raise ValueError('Read-support denominators changed across profiles')
            y = [100 * float(v['same_fraction']) for v in rows_for_allele]
            a.plot(range(3), y, marker=marker, color=color, ls=linestyle,
                   lw=1.6, ms=5, label=allele,
                   markerfacecolor='white' if allele == 'ALT' else color)
            final = rows_for_allele[-1]
            a.text(.035, y_text,
                   f"Final {allele}: {int(final['same']):,}/{int(final['baseline_retained']):,}",
                   transform=a.transAxes, fontsize=9, color=color)
        a.set_title(chr(97 + i) + '  ' + label, loc='left', fontweight='bold', fontsize=10.8, pad=16)
        a.text(0, 1.015, f'{n} sites', transform=a.transAxes, fontsize=8.5, color='#444444')
        a.set(ylim=(0, 107), xlim=(-.15, 2.2), xticks=[0, 1, 2],
              xticklabels=['MQ20\nBAQ off', 'MQ40\nBAQ off', 'MQ40\nBAQ recalculated'])
        a.set_ylabel('Retained calls supported\nas the same allele (%)', fontsize=9)
        a.grid(axis='y', color='#DDDDDD', lw=.6)
        a.set_axisbelow(True)
        if i == 0:
            a.legend(loc='center left', frameon=False, title='Baseline allele', fontsize=9, title_fontsize=9)
    fig.savefig(out / 'Figure_2_allele_support.png', dpi=300)
    fig.savefig(out / 'Figure_2_allele_support.pdf')
    plt.close(fig)
    print('Created two figures in PNG and vector PDF')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    draw(args.tables, args.out)


if __name__ == '__main__':
    main()
