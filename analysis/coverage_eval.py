"""coverage_eval.py -- empirical coverage of the refit profile intervals on simulated cases.
For each trunk clone: does the fitness interval (binomial and quasi) contain the true clone fitness,
and does the t_e interval contain the true founding time (the quantity the paper interprets)?
Also median interval widths, for comparison with the real-data intervals.
Usage: python3 -u coverage_eval.py <cov.tsv> <truth.tsv>"""
import sys, csv
if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
lines = open(sys.argv[1]).read().splitlines()
phi = lines[0]
fit = {(r['person'], r['clone']): r for r in csv.DictReader(lines[1:], delimiter='\t')}
truth = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
def med(v):
    v = sorted(v); return v[len(v)//2] if v else float('nan')
print(phi)
for k in ('1', '2', '3', '4', 'all'):
    rows = [(fit[(t['person'], t['clone'])], t) for t in truth if t['trunk'] == '1' and (k == 'all' or t['nmut'] == k)
            and (t['person'], t['clone']) in fit]
    if not rows: continue
    n = len(rows)
    cs = sum(float(f['s_lo']) <= float(t['fit']) <= float(f['s_hi']) for f, t in rows)
    csq = sum(float(f['s_lo_q']) <= float(t['fit']) <= float(f['s_hi_q']) for f, t in rows)
    ct = sum(float(f['te_lo']) <= float(t['occ']) <= float(f['te_hi']) for f, t in rows)
    ctq = sum(float(f['te_lo_q']) <= float(t['occ']) <= float(f['te_hi_q']) for f, t in rows)
    w = [float(f['te_hi_q']) - float(f['te_lo_q']) for f, _ in rows]
    print(f"trunk clones with {k} drivers (n={n}): fitness coverage binomial {cs/n:.2f}, quasi {csq/n:.2f}; "
          f"founding-time coverage by t_e interval binomial {ct/n:.2f}, quasi {ctq/n:.2f}; median quasi t_e width {med(w):.1f} yr")
