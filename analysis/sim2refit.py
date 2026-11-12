"""sim2refit.py -- turn simulated AML cases (amlsim --dump-cases) into refit input files, with
simulated sequencing, so the paper's estimator can be scored against the true mutation times.

Design: annual samples at ages dx-NS, ..., dx-1 (NS = 8 by default); VAF of a mutation =
0.5 * (cells of clones carrying it) / (all cells); alt reads ~ Binomial(depth, VAF); a variant is
'detected' if its true VAF >= 0.001 at some sample; an observation with fewer than 3 alt reads is
written as undetected (censored, as refit treats it).  The clone structure given to the estimator
is the true tree of the detected clones (best case for the estimator).
Also writes truth.tsv: person, clone name, nmut, true occurrence time, true fitness, on-AML-trunk.

Usage: python3 -u sim2refit.py <cases.bin> <outdir> [--samples NS] [--depth D] [--seed S] [--max M]
"""
import sys, os, struct
if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
import numpy as np
src, out = sys.argv[1], sys.argv[2]
args = sys.argv[3:]
NS = int(args[args.index('--samples') + 1]) if '--samples' in args else 8
DEPTH = int(args[args.index('--depth') + 1]) if '--depth' in args else 1500
SEED = int(args[args.index('--seed') + 1]) if '--seed' in args else 1
MAXC = int(args[args.index('--max') + 1]) if '--max' in args else 10**9
rng = np.random.default_rng(SEED)
os.makedirs(out, exist_ok=True)
b = open(src, 'rb').read()
off = 0
nc, T, dt = struct.unpack_from('<qid', b, off); off += 8 + 4 + 8
truth = open(os.path.join(out, 'truth.tsv'), 'w')
truth.write('person\tclone\tnmut\tocc\tfit\ttrunk\tdx\n')
written = 0
for ci in range(nc):
    sim, dx, ncl = struct.unpack_from('<qdi', b, off); off += 8 + 8 + 4
    cl = []
    for _ in range(ncl):
        cid, par, nmut = struct.unpack_from('<iii', b, off); off += 12
        muts = struct.unpack_from('<4i', b, off); off += 16
        incs = struct.unpack_from('<4d', b, off); off += 32
        occ, fit = struct.unpack_from('<dd', b, off); off += 16
        aml = b[off]; off += 1
        traj = np.frombuffer(b, dtype=np.float32, count=T, offset=off).astype(float); off += 4 * T
        cl.append(dict(id=cid, parent=par, nmut=nmut, muts=muts[:nmut], occ=occ, fit=fit, aml=aml, traj=traj))
    if written >= MAXC:
        continue
    times = [dx - NS + j for j in range(NS)]
    if times[0] < 1:
        continue
    idx = [int(round(t / dt)) for t in times]
    pop = sum(c['traj'] for c in cl)[idx]
    trunk = set()
    for c in cl:
        if c['aml']:
            trunk = set(c['muts'])
            break
    # lineage VAF of each mutation
    allm = sorted({m for c in cl for m in c['muts']})
    vaf = {m: 0.5 * sum(c['traj'][idx] for c in cl if m in c['muts']) / pop for m in allm}
    det = [m for m in allm if np.max(vaf[m]) >= 0.001]
    clones = [c for c in cl if c['nmut'] >= 1 and c['muts'][-1] in det]
    if not clones:
        continue
    pid = 'S%06d' % sim
    lines = [f'PERSON {pid} {dx:.3f}', 'TIMES ' + ' '.join(f'{t:.4f}' for t in times)]
    for m in det:
        obs = []
        for j in range(NS):
            k = int(rng.binomial(DEPTH, min(max(vaf[m][j], 0.0), 1.0)))
            obs.append(f'{k}/{DEPTH}' if k >= 3 else f'U{DEPTH}')
        lines.append(f'VAR m{m} S ' + ' '.join(obs))
    for c in clones:
        name = ' | '.join(f'm{m}' for m in c['muts'])
        lines.append('CLONE ' + name)
        truth.write(f"{pid}\t{','.join('m%d' % m for m in c['muts'])}\t{c['nmut']}\t{c['occ']:.2f}\t{c['fit']:.4f}\t{int(c['muts'][-1] in trunk)}\t{dx:.1f}\n")
    open(os.path.join(out, pid + '.txt'), 'w').write('\n'.join(lines) + '\n')
    written += 1
truth.close()
print(f'{written} simulated persons written to {out} (of {nc} cases)')
