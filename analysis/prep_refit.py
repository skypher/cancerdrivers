"""prep_refit.py -- build per-person input files for refit.cpp from the public data of
Watson et al. (Nature 2026), repository the-blundell-lab/preAML_evolutionary_dynamics, Data_files/.

Per person: sample ages (dx_age_rounded - months_to_diagnosis/12; relative times exact to a
month, absolute offset < 1 yr), the clone structure of Clonal_structure_for_each_sample.csv,
and per variant and timepoint either read counts (SNV/indel: alt, depth), a cell fraction with
s.d. (mCA: phased-caller 95% CI / 3.92, or 0.03 when no CI), or 'U depth' for an undetected
timepoint (depth from UKCTOCS_read_depths_at_variant_positions.csv, else 1500).

Usage: python3 -u prep_refit.py <Data_files dir> <out dir>
"""
import sys, os, csv, collections

if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
    print(__doc__)
    sys.exit(0)
D, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
EXCLUDE = {'C92_041_s2', 'C92_044_s4', 'C92_015_s6', 'C92_054_s3', 'CNTRL_203_s5'}  # as in v16 code

rows = list(csv.DictReader(open(os.path.join(D, 'UKCTOCS_non-germline_variants_calls_SNVs_indels_mCAs_rounded_ages.csv'))))
cs = list(csv.DictReader(open(os.path.join(D, 'Clonal_structure_for_each_sample.csv'))))
dep = {}
for r in csv.DictReader(open(os.path.join(D, 'UKCTOCS_read_depths_at_variant_positions.csv'))):
    dep[(r['sample_timepoint'], r['chromosome'], r['position'])] = int(float(r['total_depth']))
mca_ci = {}
for r in csv.DictReader(open(os.path.join(D, 'Somatic_mCA_calls.csv'))):
    lo, hi = r['95% lower cell fraction (phased caller)'], r['95% upper cell fraction (phased caller)']
    if lo != '' and hi != '':
        mca_ci[(r['Sample name'], r['Chromosome'], r['Start'])] = (float(hi) - float(lo)) / 3.92


def pid_of(sample):
    return '_'.join(sample.split('_')[:2])


def vname(r):
    if r['variant_type'] in ('SNV', 'Deletion', 'Insertion', 'Complex', 'ITD'):
        return r['gene'] + '_' + r['AA_change'], 'S'
    return r['chromosome'] + '_' + r['variant_type'].replace('-', ''), 'M'


times = collections.defaultdict(dict)   # pid -> sample -> age
dxage = {}
calls = collections.defaultdict(dict)   # pid -> (vname) -> sample -> record
vinfo = {}
case_time = {}
for r in rows:   # case sample ages: rounded dx age minus exact months
    if r['age_at_AML_diagnosis'] != '' and r['months_to_diagnosis'] != '':
        case_time[r['sample name']] = (float(r['age_at_AML_diagnosis']) - float(r['months_to_diagnosis']) / 12.0,
                                       float(r['age_at_AML_diagnosis']))
for r in rows:
    s = r['sample name']
    if s in EXCLUDE:
        continue
    p = pid_of(s)
    if s in case_time:
        times[p][s], dxage[p] = case_time[s]
    elif r['matched_sample'] in case_time:   # control: timepoint-matched case sample
        times[p][s], dxage[p] = case_time[r['matched_sample']]
    else:
        times[p][s] = float(r['age_sample_taken']) + 0.5
        dxage.setdefault(p, float(r['age_sample_taken']) + 5.0)
    nm, kind = vname(r)
    vinfo[(p, nm)] = (kind, r['chromosome'], r['start'])
    calls[p].setdefault(nm, {})[s] = r

clones = collections.defaultdict(list)
for c in cs:
    clones[c['Sample name']].append([v.strip() for v in c['Variants in clone'].split(',')])

n_written = 0
for p in sorted(clones):
    ts = sorted(times[p].items(), key=lambda kv: kv[1])
    samples = [k for k, _ in ts]
    lines = []
    lines.append(f"PERSON {p} {dxage.get(p, float('nan')):.3f}")
    lines.append("TIMES " + " ".join(f"{a:.4f}" for _, a in ts))
    allv = []
    for cl in clones[p]:
        for v in cl:
            if v not in allv:
                allv.append(v)
    for v in allv:
        if v in ('X', 'Y', 'Z'):
            lines.append(f"VAR {v} H")   # hidden (missing) driver
            continue
        kind, chrom, start = vinfo[(p, v)]
        obs = []
        for s in samples:
            r = calls[p][v].get(s)
            if kind == 'S':
                if r is not None and r['variant_depth'] not in ('', 'not processed'):
                    obs.append(f"{int(float(r['variant_depth']))}/{int(float(r['total_depth']))}")
                else:
                    d = dep.get((s, chrom, start), 1500)
                    obs.append(f"U{d}")
            else:
                if r is not None and r['VAF (cell fraction for mCAs)'] not in ('', 'not processed'):
                    f = float(r['VAF (cell fraction for mCAs)'])
                    sd = mca_ci.get((s, chrom.rstrip('pq'), start), 0.03)
                    obs.append(f"F{f:.5f}:{max(sd, 0.005):.5f}")
                else:
                    obs.append("NA")
        lines.append(f"VAR {v} {kind} " + " ".join(obs))
    for i, cl in enumerate(clones[p]):
        lines.append("CLONE " + " | ".join(cl))
    open(os.path.join(OUT, p + '.txt'), 'w').write("\n".join(lines) + "\n")
    n_written += 1
print(f"wrote {n_written} person files to {OUT}")
