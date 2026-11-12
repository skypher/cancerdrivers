"""ppv.py -- PPV of 'sweep + second driver' and OR of 'sweep' at ages 40..70 from amlsim *_ppv.bin.
Usage: python3 -u ppv.py <prefix>_ppv.bin"""
import sys
if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
    print(__doc__); sys.exit(0)
import numpy as np
b = open(sys.argv[1], 'rb').read()
n = int(np.frombuffer(b[:8], dtype=np.int64)[0])
dx = np.frombuffer(b[8:8+4*n], dtype=np.float32)
rec = np.frombuffer(b[8+4*n:], dtype=np.float32).reshape(n, 8)
case = dx >= 0
for a, age in enumerate([40, 50, 60, 70]):
    v1, v2 = rec[:, 2*a], rec[:, 2*a+1]
    alive = v1 >= 0
    crit = alive & (v1 > 0.3) & (v2 > 0.1)
    sw = alive & (v1 > 0.3)
    nc, ncr = int((crit & case).sum()), int(crit.sum())
    a_, b_ = int((sw & case).sum()), int((sw & ~case).sum())
    c_, d_ = int((alive & ~sw & case).sum()), int((alive & ~sw & ~case).sum())
    OR = (a_*d_)/(b_*c_) if b_*c_ > 0 else float('nan')
    print(f"age {age}: PPV sweep+second {nc}/{ncr} = {100*nc/max(ncr,1):.1f}%   sweep OR {OR:.1f} (a={a_} b={b_} c={c_} d={d_})")
