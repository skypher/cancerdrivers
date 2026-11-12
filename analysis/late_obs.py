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
  * rule: late if no detectable trunk mutation is detected in any sample, or the first detection is
    at a sample taken <= 2 years before diagnosis.
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
# Rules (both require an earlier negative sample for a 'late detection', as the authors evidently did:
# real cases C92_065 and C92_066, sampled only in the last 2 years, were not called late):
#  A: the putative pre-leukaemic clone is the lineage of the first visible trunk driver;
#  B: it is the lineage of the second visible trunk driver (an earlier visible driver counts as
#     clonal haematopoiesis, as for real case C92_016, listed 'linear/late').
# late = no detection of that lineage in any sample, or its first detection is <= 2 yr before
# diagnosis and some earlier sample (> 2 yr) did not detect it.
for rule in ('A', 'B'):
    for p in PS:
        thr = 0.001
        late = tot = nodet = 0
        for dx, V in cases:
            for _ in range(R):
                s = S[rng.integers(len(S))]
                ages = dx - s
                ok = ages >= 0
                ss = s[ok]
                idx = np.clip(np.round(ages[ok] / dt).astype(int), 0, T - 1)
                vis = np.where(rng.random(4) >= p)[0]
                k = 0 if rule == 'A' else 1
                tot += 1
                if len(vis) <= k:
                    late += 1; nodet += 1; continue
                det = V[vis[k]][idx] >= thr
                if not det.any():
                    late += 1; nodet += 1
                elif ss[det].max() <= 2.0 and (ss > 2.0).any():
                    late += 1
        lo, hi = ci(late, tot)
        print(f"[P(>=1 of 4 trunk drivers off-panel) = {1-(1-p)**4:.2f}] rule {rule}, off-panel prob {p:.2f}: late {late}/{tot} = {100*late/tot:.2f}% [{100*lo:.2f}-{100*hi:.2f}] "
              f"(of which no detection: {100*nodet/tot:.2f}%)")
