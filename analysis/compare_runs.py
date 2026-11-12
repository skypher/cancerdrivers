"""compare_runs.py -- compare amlsim runs on the outputs the paper reports from its Fig. 5 model:
lifetime risk and age-incidence (vs CRUK female and male-female average), pattern fractions,
PPV of 'largest VAF > 30% and second > 10%' at 40-70, median founding ages of the AML trunk,
conditioned increments and clone fitness.
Usage: python3 -u compare_runs.py <label=prefix> [<label=prefix> ...]   (prefix of *_summary.txt etc.)"""
import sys, csv
if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
import numpy as np
F = np.array([0.9,0.4,0.5,0.7,0.8,1.0,1.3,1.2,1.7,2.1,2.8,3.5,5.5,7.9,10.9,13.6,20.4])
M = np.array([1.1,0.4,0.6,0.7,1.0,0.9,1.2,1.3,1.6,2.1,3.3,4.8,7.5,12.6,18.2,26.7,34.0])
A = 0.5*(F+M)
def med(v):
    v = sorted(v); return v[len(v)//2]
print('CRUK cumulative risk to 85: female %.2f%%, average %.2f%%' % (F.sum()*5/1e3, A.sum()*5/1e3))
print('CRUK incidence 60-64/70-74/80-84 (avg): %.1f / %.1f / %.1f per 100k' % (A[12], A[14], A[16]))
for arg in sys.argv[1:]:
    label, pre = arg.split('=', 1)
    summ = dict(l.split(None, 1) for l in open(pre + '_summary.txt') if l.strip())
    inc = np.array([float(x) for x in summ['incidence_per100k_per_yr'].split()])
    rows = list(csv.DictReader(open(pre + '_cases.tsv'), delimiter='\t'))
    n = len(rows)
    pat = [sum(r['pattern'] == str(k) for r in rows) for k in range(3)]
    b = open(pre + '_ppv.bin', 'rb').read()
    ns = int(np.frombuffer(b[:8], dtype=np.int64)[0])
    dx = np.frombuffer(b[8:8+4*ns], dtype=np.float32)
    rec = np.frombuffer(b[8+4*ns:], dtype=np.float32).reshape(ns, 8)
    ppv = []
    for a in range(4):
        v1, v2 = rec[:, 2*a], rec[:, 2*a+1]
        crit = (v1 >= 0) & (v1 > 0.3) & (v2 > 0.1)
        ppv.append('%d/%d=%.0f%%' % ((crit & (dx >= 0)).sum(), crit.sum(), 100*(crit & (dx >= 0)).sum()/max(crit.sum(), 1)))
    occ = [med([float(r['occ%d' % k]) for r in rows]) for k in range(1, 5)]
    incs = [med([float(r['inc%d' % k]) for r in rows]) for k in range(1, 5)]
    fit = [med([float(r['fit%d' % k]) for r in rows]) for k in range(1, 5)]
    print(f"\n== {label} ({pre}): {n} cases in {ns} lives, lifetime risk {100*n/ns:.3f}%")
    print(f"  incidence 60-64/70-74/80-84: {inc[12]:.1f} / {inc[14]:.1f} / {inc[16]:.1f} per 100k; ratio to CRUK avg {inc[12]/A[12]:.2f} / {inc[14]/A[14]:.2f} / {inc[16]/A[16]:.2f}")
    print(f"  patterns linear/branched/late: {pat[0]}/{pat[1]}/{pat[2]} = {100*pat[0]/n:.0f}%/{100*pat[1]/n:.0f}%/{100*pat[2]/n:.0f}% (paper data 26/11/10 of 47; paper model 642/301/0)")
    print(f"  PPV at 40/50/60/70: {', '.join(ppv)}  (paper: 44/50, 214/292, 564/1152, 719/3313)")
    print(f"  median founding ages drivers 1-4: {', '.join('%.1f' % x for x in occ)}; median dx {med([float(r['dx']) for r in rows]):.1f}")
    print(f"  median conditioned increments: {', '.join('%.3f' % x for x in incs)}; clone fitness: {', '.join('%.3f' % x for x in fit)}")
    print(f"  first driver by 20: {sum(float(r['occ1']) <= 20 for r in rows)/n:.3f}")
