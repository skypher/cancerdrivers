"""firstdriver_law.py -- check of the first-driver age law
    f(u | diagnosis at T) = h'(T-u) / h(T),   E[u | T] = H(T)/h(T),   P(u <= x | T) = 1 - h(T-x)/h(T)
(h = age-specific incidence, H its integral), valid when first drivers are founded at a constant rate
and lineages from different founders evolve independently and time-homogeneously.
(1) simulation: h estimated from the simulated diagnosis ages (smoothed), compared with the
    simulated first-driver ages of the AML trunks;  (2) the same formula applied to CRUK incidence.
Usage: python3 -u firstdriver_law.py <base_cases.tsv>"""
import sys, csv, math
if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
dx = [float(r['dx']) for r in rows]; u1 = [float(r['occ1']) for r in rows]
n = len(dx)
# smoothed incidence: kernel density of dx (Gaussian, bw 2.5 yr, reflected at 85); only its shape matters
bw = 2.5
def h(T):
    if T <= 0: return 0.0
    # reflection at the age-85 cut-off removes the boundary bias of the kernel estimate
    return sum(math.exp(-0.5*((T-d)/bw)**2) + math.exp(-0.5*((T-(170.0-d))/bw)**2) for d in dx)
grid = [0.25*i for i in range(0, 341)]
hg = [h(t) for t in grid]
def Hcum(T):
    s = 0.0
    for i in range(1, len(grid)):
        if grid[i] > T: break
        s += 0.25*0.5*(hg[i-1]+hg[i])
    return s
pred20 = sum(1 - h(T-20)/h(T) for T in dx)/n
predmean = sum(Hcum(T)/h(T) for T in dx)/n
print(f"simulation ({n} AML cases): P(first driver <= 20): law {pred20:.3f}, simulated {sum(x<=20 for x in u1)/n:.3f}")
print(f"simulation: mean first-driver age: law {predmean:.2f}, simulated {sum(u1)/n:.2f}")
# per-quartile of dx
qs = sorted(dx)
for lo, hi in [(0, qs[n//4]), (qs[n//4], qs[n//2]), (qs[n//2], qs[3*n//4]), (qs[3*n//4], 100)]:
    sel = [i for i in range(n) if lo <= dx[i] < hi]
    if not sel: continue
    p = sum(1 - h(dx[i]-20)/h(dx[i]) for i in sel)/len(sel)
    m = sum(Hcum(dx[i])/h(dx[i]) for i in sel)/len(sel)
    print(f"  dx in [{lo:.1f},{hi:.1f}): n={len(sel)} P(u<=20) law {p:.3f} sim {sum(u1[i]<=20 for i in sel)/len(sel):.3f}; "
          f"mean u law {m:.1f} sim {sum(u1[i] for i in sel)/len(sel):.1f}")
# CRUK incidence (per 100k per yr), 5-yr bands centred 2.5..82.5 (from figure5_pipeline.py)
F = [0.9,0.4,0.5,0.7,0.8,1.0,1.3,1.2,1.7,2.1,2.8,3.5,5.5,7.9,10.9,13.6,20.4]
Mm = [1.1,0.4,0.6,0.7,1.0,0.9,1.2,1.3,1.6,2.1,3.3,4.8,7.5,12.6,18.2,26.7,34.0]
cent = [2.5+5*i for i in range(17)]
def interp(arr, T):
    if T <= cent[0]: return arr[0]
    if T >= cent[-1]: return arr[-1]
    i = int((T-2.5)//5); w = (T-cent[i])/5
    return arr[i]*(1-w)+arr[i+1]*w
for name, arr in (('female', F), ('male', Mm)):
    for start in (0, 30):
        def hc(T): return interp(arr, T) if T >= start else 0.0
        out = []
        for T in (60, 72, 80):
            Hs = sum(hc(start + 0.05*k)*0.05 for k in range(int((T-start)/0.05)))
            out.append(f"T={T}: P(u<=20) {1-hc(T-20)/hc(T):.2f}, E[u] {Hs/hc(T):.1f}")
        print(f"CRUK {name}, incidence counted from age {start}: " + '; '.join(out))
