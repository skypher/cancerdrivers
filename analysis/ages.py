"""ages.py -- first-driver acquisition ages implied by the published (s, t_e) estimates
(All_fitness_establishment_times_with_sd_and_CI.csv) once the random establishment factor X
and the population scale N*tau are accounted for.

Model: an established clone of fitness s founded at t_mut has n(t) ~ X exp(s (t - t_mut)) / (2 s tau)
cells, X ~ Exp(1).  The paper's fit uses size (exp(s(t-t0))-1)/s against N_inf = 1e5, and reports
t_e = t0 + log(1+s)/s.  Hence t_mut = t0 + (log X + kappa)/s with kappa = log(N_inf s / (q N)).
Prior on t_mut: uniform on [0, inf) (constant first-driver supply), so the posterior of
L = log X is its Gumbel-min law truncated to t_mut >= 0.
Usage: python3 -u ages.py <csv>"""
import sys, csv, math
if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
rows = list(csv.DictReader(open(sys.argv[1])))
by = {}
for r in rows:
    by.setdefault(r['Sample_name'], []).append(r)
def FL(l):  # CDF of log X, X ~ Exp(1)
    if l > 50: return 1.0
    return 1.0 - math.exp(-math.exp(l))
def p_le(age, t0, s, kappa):
    # P(t_mut <= age | t0, s), t_mut = t0 + (L + kappa)/s, truncated to t_mut >= 0
    lo = FL(-s*t0 - kappa)
    hi = FL(s*(age - t0) - kappa)
    return max(0.0, (hi - lo) / (1.0 - lo)) if lo < 1 else 1.0
def mean_tmut(t0, s, kappa, n=4000):
    # posterior mean by quadrature in L
    num = den = 0.0
    lmin = -s*t0 - kappa
    for i in range(n):
        l = -12 + 16*(i+0.5)/n
        if l < lmin: continue
        w = math.exp(l - math.exp(l))
        num += w*(t0 + (l + kappa)/s); den += w
    return num/den
Ninf = 1e5
for label, sel in (('earliest clone with >=1 driver (Fig. 4b rule)', 'any'),
                   ('first clone of the AML lineage', 'lineage')):
    recs = []
    for k in sorted(by):
        if not k.startswith('C92'): continue
        cl = by[k]
        if sel == 'any':
            first = min(cl, key=lambda r: float(r['Establishment_time']))
        else:
            fittest = max(cl, key=lambda r: float(r['Fitness']))
            aml = set(v.strip() for v in fittest['Clone'].split(','))
            lin = [r for r in cl if r['Clone'].split(',')[-1].strip() in aml]
            first = min(lin, key=lambda r: float(r['Establishment_time']))
        s = float(first['Fitness'])/100.0
        te = float(first['Establishment_time'])
        t0 = float(first['Tau_time'])
        recs.append((k, s, te, t0))
    n = len(recs)
    print(f"\n== {label}: {n} cases with a fitted clone (paper's denominator for Fig. 4b is 47 cases)")
    print(f"  as reported: t_e <= 20 in {sum(r[2] <= 20 for r in recs)}/{n};  t_e pinned near birth (t_e <= 1.5): {sum(r[2] <= 1.5 for r in recs)}")
    # kappa = log(N_inf s / (q N)); q ~ s*tau in the convention of ref 39 (Watson 2020, Suppl. Note 2:
    # self-renewal and differentiation each at rate 1/tau), q ~ 2 s C/N with C the coalescent N_e*g
    # (refs 38, 41). The code has tau = 2 dt (ref 39 convention), C = N dt.
    for label2, kappa in (('code (ref-39 Ntau = 2e4, coalescent 1e4)', math.log(5.0)),
                          ('stated Ntau = 1e5, ref-39 convention', 0.0),
                          ('stated Ntau = 1e5, coalescent convention', math.log(0.5))):
        Ntau = label2
        e20 = sum(p_le(20, r[3], r[1], kappa) for r in recs)
        shift = sum(mean_tmut(r[3], r[1], kappa) - r[2] for r in recs)/n
        print(f"  {Ntau} (kappa = {kappa:+.3f}): expected number with true t_mut <= 20: {e20:.1f}/{n}; "
              f"mean (posterior mean t_mut - reported t_e) = {shift:+.1f} yr")
    # spread
    sds = [1.2825/r[1] for r in recs]
    sds.sort()
    print(f"  irreducible s.d. of t_mut given (s, t0), pi/sqrt(6)/s: median {sds[n//2]:.1f} yr, range {sds[0]:.1f}-{sds[-1]:.1f} yr")
