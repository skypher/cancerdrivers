"""late_obs.py -- observation-matched 'late evolution' classification of simulated AML cases.

The simulation's own classifier (figure5_pipeline.classify) calls a case late when the AML founder's
lineage first exceeds 0.2% cell fraction within 2 years of diagnosis, using the complete trajectory.
The paper classifies real cases from sparse samples and a finite panel: late = the putative
pre-leukaemic clone is detectable only within 1-2 years before diagnosis, or no driver is detected.
This script applies the real observation process to simulated AML cases:
  * sampling: each simulated case gets the sample schedule (years before diagnosis) of a randomly
    drawn real case (47 cases; exact months where the public calls file has them, otherwise the
    midpoint of the sample's timing group in UKCTOCS_samples_processed_minimal_information.csv);
  * detection: a trunk mutation is detected at a sample if its lineage VAF >= threshold;
  * panel: each trunk driver is off-panel (never detectable) independently with probability p;
  * rule: late if the selected lineage is not detected in any sample, or its first detection is at a
    sample taken <= 2 years before diagnosis AND some earlier sample (> 2 yr) did not detect it.
    Rule A selects the first visible trunk driver's lineage, rule B the second visible one. These are
    sensitivity scenarios, not the paper's classifier (which is not specified as an algorithm).
Intervals: case-level bootstrap over the stored simulated AML histories (each case's late
probability averaged over its schedule/panel draws), i.e. for a population of AML histories from the
model, not Clopper-Pearson on reused draws.
Usage: python3 -u late_obs.py <cases.bin> <Data_files dir> [draws per case] [comma-separated off-panel probabilities]"""
import sys, struct, csv, os
if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
import numpy as np
from scipy.stats import beta
src, D = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 20
PS = [float(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else [0.0, 0.1, 0.2, 0.3, 0.5]
rng = np.random.default_rng(20261112)
# ---- real schedules
exact = {}
for r in csv.DictReader(open(os.path.join(D, 'UKCTOCS_non-germline_variants_calls_SNVs_indels_mCAs_rounded_ages.csv'))):
    if r['sample name'].startswith('C92') and r['months_to_diagnosis'] not in ('', None):
        exact[r['sample name']] = float(r['months_to_diagnosis']) / 12.0
mid = {'0-6 mths preDx': 0.25, '3-6 mths preDx': 0.375, '6-9 mths preDx': 0.625, '6-12 mths preDx': 0.75,
       '0-1 yr preDx': 0.5, '>10 yrs preDx': 11.0}
for k in range(1, 10):
    mid[f'{k}-{k+1} yrs preDx'] = k + 0.5
EXCL = {'C92_041_s2', 'C92_044_s4', 'C92_015_s6', 'C92_054_s3'}
sched = {}
for r in csv.DictReader(open(os.path.join(D, 'UKCTOCS_samples_processed_minimal_information.csv'))):
    if r['GROUP'] != 'Case' or r['SAMPLE ID'] in EXCL:
        continue
    sid = r['SAMPLE ID']; pid = '_'.join(sid.split('_')[:2])
    t = exact.get(sid, mid.get(r['Sample timing group']))
    if t is None:
        raise SystemExit('unknown timing group ' + r['Sample timing group'])
    sched.setdefault(pid, []).append(t)
S = [np.array(sorted(v)) for v in sched.values()]
lastsample = np.array([s.min() for s in S])
print(f"real schedules: {len(S)} cases, {sum(len(s) for s in S)} samples; median samples per case {np.median([len(s) for s in S]):.0f}; "
      f"last sample > 2 yr before diagnosis in {np.sum(lastsample > 2)} of {len(S)} cases")
# ---- simulated cases: lineage VAF of each trunk mutation on the full grid
b = open(src, 'rb').read(); off = 0
nc, T, dt = struct.unpack_from('<qid', b, off); off += 20
cases = []
for ci in range(nc):
    sim, dx, ncl = struct.unpack_from('<qdi', b, off); off += 20
    cl = []
    for _ in range(ncl):
        cid, par, nmut = struct.unpack_from('<iii', b, off); off += 12
        muts = struct.unpack_from('<4i', b, off); off += 16
        off += 32 + 16
        aml = b[off]; off += 1
        traj = np.frombuffer(b, dtype=np.float32, count=T, offset=off).astype(float); off += 4 * T
        cl.append((muts[:nmut], aml, traj))
    pop = sum(c[2] for c in cl); pop[pop <= 0] = 1.0
    trunk = [c for c in cl if c[1]][0][0]
    V = np.array([0.5 * sum(c[2] for c in cl if m in c[0]) / pop for m in trunk])  # 4 x T
    cases.append((dx, V))
print(f"simulated AML cases: {len(cases)}")
def ci(k, n):
    return (beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(0.975, k + 1, n - k) if k < n else 1.0)
lo0, hi0 = ci(10, 47); print(f"observed: 10/47 late = 21.3% [{100*lo0:.1f}-{100*hi0:.1f}]")
# Rules (both require an earlier negative sample for a 'late detection'):
#  A: the putative pre-leukaemic clone is the lineage of the first visible trunk driver;
#  B: it is the lineage of the second visible trunk driver.
# Also counted, for comparison with the paper's 'no pathogenic driver detected at all': no visible
# trunk driver detected in any sample (off-trunk clones are not in this file, so this is an upper
# bound on the true no-driver rate).
def boot(xs, B=4000):
    xs = np.asarray(xs); n = len(xs)
    m = np.array([xs[rng.integers(n, size=n)].mean() for _ in range(B)])
    return np.quantile(m, 0.025), np.quantile(m, 0.975)
for p in PS:
    thr = 0.001
    case_late = {'A': [], 'B': []}; case_none = []; case_selmiss = {'A': [], 'B': []}
    for dx, V in cases:
        cnt = {'A': 0, 'B': 0}; sm = {'A': 0, 'B': 0}; none = 0
        for _ in range(R):
            s = S[rng.integers(len(S))]
            ages = dx - s
            ok = ages >= 0
            ss = s[ok]
            idx = np.clip(np.round(ages[ok] / dt).astype(int), 0, T - 1)
            vis = np.where(rng.random(4) >= p)[0]
            detv = [(V[v][idx] >= thr) for v in vis]
            if not any(d.any() for d in detv):
                none += 1
            for rule, k in (('A', 0), ('B', 1)):
                if len(vis) <= k:
                    cnt[rule] += 1; sm[rule] += 1; continue
                det = detv[k]
                if not det.any():
                    cnt[rule] += 1; sm[rule] += 1
                elif ss[det].max() <= 2.0 and (ss > 2.0).any():
                    cnt[rule] += 1
        for rule in 'AB':
            case_late[rule].append(cnt[rule] / R); case_selmiss[rule].append(sm[rule] / R)
        case_none.append(none / R)
    lo, hi = boot(case_none)
    print(f"[P(>=1 of 4 trunk drivers off-panel) = {1-(1-p)**4:.3f}] off-panel prob {p:.4f}: "
          f"no visible trunk driver detected at all {100*np.mean(case_none):.2f}% [{100*lo:.2f}-{100*hi:.2f}]")
    for rule in 'AB':
        lo, hi = boot(case_late[rule])
        print(f"    rule {rule}: late {100*np.mean(case_late[rule]):.2f}% [{100*lo:.2f}-{100*hi:.2f}] "
              f"(selected lineage absent or undetected: {100*np.mean(case_selmiss[rule]):.2f}%)")
