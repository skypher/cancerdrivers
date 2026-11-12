"""refit_eval.py -- compare the re-fit of the public data (refit.cpp) with the published estimates.
Reports: pooled dispersion; likelihood at the published values; interval widths (profile, quasi-
binomial) against the published 95% optimiser ranges; intervals that reach the search bounds
(parameter not identified on that side); published values outside the quasi intervals; and the
first-driver ages behind the Fig. 4b statement.
Usage: python3 -u refit_eval.py <refit_full.tsv> <All_fitness_establishment_times_with_sd_and_CI.csv>"""
import sys, csv, math
if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
lines = open(sys.argv[1]).read().splitlines()
phi = float(lines[0].split('=')[1])
R = list(csv.DictReader(lines[1:], delimiter='\t'))
pub = {}
for r in csv.DictReader(open(sys.argv[2])):
    pub[(r['Sample_name'], ','.join(v.strip() for v in r['Clone'].split(',')))] = r
def rng(sv):
    a, b = sv.split(' - '); return float(a), float(b)
def med(v):
    v = sorted(v); return v[len(v)//2] if v else float('nan')
print(f"pooled Pearson dispersion phi = {phi:.2f}  (binomial model: 1)")
cases = [r for r in R if r['person'].startswith('C92')]
ctrls = [r for r in R if r['person'].startswith('CNTRL')]
for label, rows in (('cases', cases), ('controls', ctrls)):
    persons = {}
    for r in rows: persons[r['person']] = (float(r['nll']), float(r['nll_published']) if r['nll_published'] != 'nan' else float('nan'))
    d = [p - m for m, p in persons.values() if not math.isnan(p)]
    print(f"\n== {label}: {len(persons)} persons, {len(rows)} clones")
    if d:
        print(f"  NLL(published) - NLL(MLE), best common time shift: median {med(d):.1f}; >{1.92*phi:.1f} (quasi 95% for one parameter) in "
              f"{sum(x > 1.92*phi for x in d)}/{len(d)} persons; >{1.92*phi*4:.0f} in {sum(x > 1.92*phi*4 for x in d)}")
    n = len(rows)
    lo_te = sum(r['te_bound'] in ('1', '3') for r in rows); hi_te = sum(r['te_bound'] in ('2', '3') for r in rows)
    lo_s = sum(r['s_bound'] in ('1', '3') for r in rows); hi_s = sum(r['s_bound'] in ('2', '3') for r in rows)
    print(f"  quasi 95% interval for t_e reaches birth (t_e=0): {lo_te}/{n}; reaches diagnosis: {hi_te}/{n}")
    print(f"  quasi 95% interval for s reaches 0.002/yr: {lo_s}/{n}; reaches 15/yr: {hi_s}/{n}")
    wq, wo, ins, insfit, nn = [], [], 0, 0, 0
    for r in rows:
        p = pub.get((r['person'], r['clone']))
        if not p: continue
        try:
            a, b = rng(p['95% optimizer result range establishment age'])
            fa, fb = rng(p['95% optimizer result range fitness'])
        except Exception:
            continue
        nn += 1
        wq.append(float(r['te_hi_q']) - float(r['te_lo_q'])); wo.append(b - a)
        te_pub = float(p['Establishment_time']); s_pub = float(p['Fitness'])/100
        ins += float(r['te_lo_q']) - 1 <= te_pub <= float(r['te_hi_q']) + 1
        insfit += float(r['s_lo_q']) <= s_pub <= float(r['s_hi_q'])
    if nn:
        print(f"  width of 95% interval for t_e: profile (quasi) median {med(wq):.1f} yr vs published optimiser range median {med(wo):.1f} yr "
              f"({sum(x > y for x, y in zip(wq, wo))}/{nn} wider)")
        print(f"  published t_e inside the quasi interval (+-1 yr for the age rounding): {ins}/{nn}; published s inside: {insfit}/{nn}")
# Fig. 4b: earliest clone per case
first = {}
for r in cases:
    if r['person'] not in first or float(r['te']) < float(first[r['person']]['te']): first[r['person']] = r
n = len(first)
sure = sum(float(r['te_hi_q']) <= 20 for r in first.values())
maybe = sum(float(r['te_lo_q']) <= 20 for r in first.values())
mle = sum(float(r['te']) <= 20 for r in first.values())
print(f"\nFig. 4b (earliest clone per case, {n} cases; paper: about two thirds of 47 by age 20):")
print(f"  MLE t_e <= 20: {mle}; whole quasi 95% interval <= 20: {sure}; interval reaches below 20: {maybe}")
for r in sorted(first.values(), key=lambda r: float(r['te'])):
    print(f"  {r['person']} {r['clone'][:45]:45s} s {float(r['s']):.3f} [{float(r['s_lo_q']):.3f},{float(r['s_hi_q']):.3f}]  "
          f"t_e {float(r['te']):5.1f} [{float(r['te_lo_q']):5.1f},{float(r['te_hi_q']):5.1f}] bounds s{r['s_bound']} t{r['te_bound']}")
