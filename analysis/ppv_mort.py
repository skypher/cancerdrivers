"""ppv_mort.py -- PPV of 'largest VAF > 30% and second VAF > 10%' for AML with competing mortality.
Uses an amlsim *_ppv.bin produced with --years 100 (diagnosis ages up to 100). For a person meeting
the criterion at age a who is diagnosed at age T, the AML counts with weight S(a -> T), the
probability of surviving other causes from a to T, under Gompertz mortality mu(x) = A exp(B x).
A is set so that life expectancy at 70 equals E70 (default 16.5 years, roughly UK women; an
assumption, not fitted to a life table). Reports the PPV for several B.
Usage: python3 -u ppv_mort.py <y100_ppv.bin> [E70]"""
import sys, math
if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
import numpy as np
E70 = float(sys.argv[2]) if len(sys.argv) > 2 else 16.5
b = open(sys.argv[1], 'rb').read()
n = int(np.frombuffer(b[:8], dtype=np.int64)[0])
dx = np.frombuffer(b[8:8+4*n], dtype=np.float32).astype(float)
rec = np.frombuffer(b[8+4*n:], dtype=np.float32).reshape(n, 8)
def surv(A, B, x, y):
    return np.exp(-(A/B)*(np.exp(B*y) - np.exp(B*x)))
def e_at(A, B, x):
    t = np.linspace(0, 80, 8001)
    return np.trapz(surv(A, B, x, x + t), t)
for B in (0.085, 0.10, 0.115):
    lo, hi = 1e-8, 1e-1
    for _ in range(200):   # bisection on log A for e(70) = E70
        mid = math.sqrt(lo*hi)
        if e_at(mid, B, 70.0) > E70: lo = mid
        else: hi = mid
    A = math.sqrt(lo*hi)
    out = []
    for k, age in enumerate([40, 50, 60, 70]):
        v1, v2 = rec[:, 2*k], rec[:, 2*k+1]
        crit = (v1 >= 0) & (v1 > 0.3) & (v2 > 0.1)
        case = crit & (dx >= 0)
        w = surv(A, B, age, np.where(case, dx, age))
        out.append(f"age {age}: {w[case].sum()/crit.sum()*100:.0f}% (horizon 85: {((case) & (dx < 85)).sum()/crit.sum()*100:.0f}%, 100: {case.sum()/crit.sum()*100:.0f}%)")
    print(f"Gompertz B={B}, A={A:.3g} (e70={e_at(A, B, 70):.1f}, e40={e_at(A, B, 40):.1f}): " + '; '.join(out))
