"""trunkstats.py -- statistics of the AML trunk in simulated cases (amlsim --dump-cases):
frequency of the parent clone (exclusive) and of the parent lineage when each trunk driver arises,
and the mean fitness at that moment.
Usage: python3 -u trunkstats.py <cases.bin>"""
import sys, struct
if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
import numpy as np
b = open(sys.argv[1], 'rb').read(); off = 0
nc, T, dt = struct.unpack_from('<qid', b, off); off += 20
res = {k: [] for k in (2, 3, 4)}
for ci in range(nc):
    sim, dx, ncl = struct.unpack_from('<qdi', b, off); off += 20
    cl = []
    for _ in range(ncl):
        cid, par, nmut = struct.unpack_from('<iii', b, off); off += 12
        muts = struct.unpack_from('<4i', b, off); off += 16
        off += 32
        occ, fit = struct.unpack_from('<dd', b, off); off += 16
        aml = b[off]; off += 1
        traj = np.frombuffer(b, dtype=np.float32, count=T, offset=off).astype(float); off += 4 * T
        cl.append((cid, nmut, muts[:nmut], occ, fit, aml, traj))
    pop = sum(c[6] for c in cl)
    xbar = sum(c[4] * c[6] for c in cl) / pop
    A = [c for c in cl if c[5]][0]
    byid = {c[0]: c for c in cl}
    for k in (2, 3, 4):
        d = byid.get(A[2][k - 1]); p = byid.get(A[2][k - 2])
        if d is None or p is None: continue
        t = int(round(d[3] / dt))
        lin = sum(c[6][t] for c in cl if A[2][k - 2] in c[2]) / pop[t]
        res[k].append((p[6][t] / pop[t], lin, xbar[t], p[4], d[4]))
for k in (2, 3, 4):
    a = np.array(res[k])
    q = lambda v: '%.3f [%.3f, %.3f]' % (np.median(v), np.quantile(v, .25), np.quantile(v, .75))
    print(f"driver {k} arises: parent clone freq {q(a[:,0])}; parent lineage freq {q(a[:,1])}; "
          f"mean fitness {q(a[:,2])}; parent F {q(a[:,3])}; P(parent clone freq>0.5) {np.mean(a[:,0]>0.5):.2f}")
