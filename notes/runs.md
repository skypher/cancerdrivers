# Reproduction commands

Build (clang only):

    cd sim && mkdir -p bin
    for p in amlsim refit theory_branch theory_renewal theory_reduced wlaw; do
      clang++ -O3 -march=native -fopenmp=libomp -std=c++17 $p.cpp -o bin/$p; done

Data: `git clone https://github.com/the-blundell-lab/preAML_evolutionary_dynamics sources/code/preAML_evolutionary_dynamics`
(commit 5a53fe7). All programs print usage with `--help`.

| result | command (from the repo root) | time |
|---|---|---|
| reproduction of Fig. 5 numbers | `sim/bin/amlsim --sims 200200 --threads 8 --out runs/repro/base --dump-cases` then `python3 -u analysis/ppv.py runs/repro/base_ppv.bin` | 2 min |
| scale tests | `amlsim --sims 100000 --dt 0.05`; `--N 5e5 --dt 0.02`; `--sims 4000 --N 1e6` | 2-7 min each |
| horizon 100 | `amlsim --sims 200200 --years 100 --out runs/horizon/y100` | 3 min |
| refit of public data | `python3 -u analysis/prep_refit.py <Data_files> runs/refit/input`; `sim/bin/refit runs/refit/full/refit_full.tsv runs/refit/input/*.txt --threads 6 --published runs/refit/published.tsv`; `python3 -u analysis/refit_eval.py ...` | 6 min |
| estimator on simulated cases | `python3 -u analysis/sim2refit.py runs/repro/base_cases.bin runs/simfit/input`; `sim/bin/refit runs/simfit/simfit.tsv runs/simfit/input/S*.txt --no-profiles --threads 8`; `python3 -u analysis/simfit_eval.py runs/simfit/simfit.tsv runs/simfit/input/truth.tsv runs/simfit/input` | 25 min |
| published estimates and Prop. 4 | `python3 -u analysis/ages.py <Data_files>/All_fitness_establishment_times_with_sd_and_CI.csv` | 1 s |
| first-driver law | `python3 -u analysis/firstdriver_law.py runs/repro/base_cases.tsv` | 5 s |
| trunk statistics | `python3 -u analysis/trunkstats.py runs/repro/base_cases.bin` | 20 s |
| branching recursion | `sim/bin/theory_branch` | 3 s |
| reduced model | `sim/bin/theory_reduced --samples 3000000 --background 0` (and `1`) | 2 min each |
| establishment factor | `sim/bin/wlaw 20000` | 1 min |

`runs/refit/published.tsv` is written from the published CSV with columns person, clone, s
(fraction per year), t0 (`Tau_time`).
