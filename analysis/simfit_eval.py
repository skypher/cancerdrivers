"""simfit_eval.py -- score the paper's estimator (refit MLE, eqs. 1-2, binomial likelihood) on
simulated AML cases with known truth.
Inputs: refit output on simulated persons, truth.tsv from sim2refit.py.
Reports, for the first driver of the AML trunk: error of the reported age t_e against the true
founding time, compared with Proposition 4; the analogue of the paper's Fig. 4b statement; fitness
recovery by number of drivers; and founders already swept at the first sample.
Usage: python3 -u simfit_eval.py <simfit.tsv> <truth.tsv> <inputdir>"""
import sys, csv, math, os
if len(sys.argv) < 4 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
fit = {}
for r in csv.DictReader((l for l in open(sys.argv[1]) if not l.startswith('#')), delimiter='\t'):
    fit[(r['person'], r['clone'])] = r
truth = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
indir = sys.argv[3]
G = 0.5772156649
def q(F, h=0.1):
    m = math.exp(F*h); x = 2*(m-1)/m
    for _ in range(100):
        g = 1-x-math.exp(-m*x); dg = -1+m*math.exp(-m*x); x -= g/dg
    return x
def med(v):
    v = sorted(v); n = len(v)
    return v[n//2] if n else float('nan')
def mean(v): return sum(v)/len(v) if v else float('nan')
def sd(v):
    m = mean(v); return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1)) if len(v) > 1 else float('nan')
# first lineage VAF at first sample, to flag founders swept before sampling
def first_vaf(person, mname):
    for line in open(os.path.join(indir, person + '.txt')):
        if line.startswith('VAR ' + mname + ' '):
            tok = line.split()[3]
            if tok.startswith('U'): return 0.0
            k, n = tok.split('/'); return int(k)/int(n)
    return float('nan')
err, pred, s_hat, s_true, swept_err, notswept_err = [], [], [], [], [], []
te_le20 = occ_le20 = n1 = 0
persons_first = {}
for t in truth:
    key = (t['person'], t['clone'])
    if key not in fit: continue
    f = fit[key]
    s, te, occ, F = float(f['s']), float(f['te']), float(t['occ']), float(t['fit'])
    if t['trunk'] == '1' and t['nmut'] == '1':
        n1 += 1
        e = te - occ
        err.append(e)
        kappa = math.log(F/q(F))  # N_inf = N
        pred.append((G - kappa + math.log(1+F))/F)
        te_le20 += te <= 20; occ_le20 += occ <= 20
        v = first_vaf(t['person'], t['clone'])
        (swept_err if v > 0.4 else notswept_err).append((e, s, F))
    persons_first.setdefault(t['person'], []).append((te, occ))
print(f"AML-trunk first drivers fitted: {n1}")
print(f"  reported age minus true founding age: mean {mean(err):+.1f} yr, median {med(err):+.1f}, s.d. {sd(err):.1f}")
print(f"  Proposition 4 (true fitness, X ~ Exp(1)): predicted mean {mean(pred):+.1f} yr")
print(f"  reported t_e <= 20: {te_le20/n1:.3f}; true founding <= 20: {occ_le20/n1:.3f}")
if notswept_err:
    print(f"  founders not swept at first sample (n={len(notswept_err)}): error mean {mean([e for e,_,_ in notswept_err]):+.1f}, s.d. {sd([e for e,_,_ in notswept_err]):.1f}; "
          f"median s_hat/F {med([s/F for _,s,F in notswept_err]):.2f}")
if swept_err:
    print(f"  founders swept at first sample, VAF>0.4 (n={len(swept_err)}): error mean {mean([e for e,_,_ in swept_err]):+.1f}, s.d. {sd([e for e,_,_ in swept_err]):.1f}; "
          f"median s_hat/F {med([s/F for _,s,F in swept_err]):.2f}; IQR of s_hat/F "
          f"{sorted([s/F for _,s,F in swept_err])[len(swept_err)//4]:.2f}-{sorted([s/F for _,s,F in swept_err])[3*len(swept_err)//4]:.2f}")
# Fig. 4b rule: earliest fitted clone per person
e20 = sum(1 for p, v in persons_first.items() if min(v)[0] <= 20)
t20 = sum(1 for p, v in persons_first.items() if min(o for _, o in v) <= 20)
print(f"Fig. 4b rule (earliest fitted clone per person, {len(persons_first)} persons): reported <=20 {e20/len(persons_first):.3f}; true earliest founding <=20 {t20/len(persons_first):.3f}")
# fitness recovery by number of drivers (trunk clones)
for k in ('1', '2', '3', '4'):
    rat = [float(fit[(t['person'], t['clone'])]['s'])/float(t['fit']) for t in truth
           if t['trunk'] == '1' and t['nmut'] == k and (t['person'], t['clone']) in fit]
    sh = [float(fit[(t['person'], t['clone'])]['s']) for t in truth if t['trunk'] == '1' and t['nmut'] == k and (t['person'], t['clone']) in fit]
    tf = [float(t['fit']) for t in truth if t['trunk'] == '1' and t['nmut'] == k and (t['person'], t['clone']) in fit]
    if rat:
        print(f"trunk clones with {k} drivers (n={len(rat)}): median fitted s {med(sh):.3f}, median true F {med(tf):.3f}, median ratio {med(rat):.2f}, "
              f"IQR {sorted(rat)[len(rat)//4]:.2f}-{sorted(rat)[3*len(rat)//4]:.2f}")
