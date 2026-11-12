# Statistics re-check of Watson et al., Nature 2026

Paper: Watson, Poon, MacGregor, Fonseca, Apostolidou, Damaskou, Boniface, Hackett,
Gentry-Maharaj, Menon, Blundell, *Evolutionary dynamics in the years preceding acute myeloid
leukaemia*, Nature, online 2026-11-11, doi:10.1038/s41586-026-11269-x.
Code and public data: github.com/the-blundell-lab/preAML_evolutionary_dynamics, commit 5a53fe7
(`Data_files/`). Controlled-access data (EGA) were not used; ages in the public tables are rounded
down to whole years (relative timing within a person is exact to a month).

Companion: `theory/theory.md` (propositions referred to below as Prop. n).

## Summary

1. **The simulation does not run at the stated population scale.** The paper states
   $N\tau=10^5$; `sim_fast.py` resamples every 0.1 yr, which makes $N\tau=N\cdot dt=10^4$
   (Corollary 3 of the theory note). At the stated $10^5$ the same model gives a 16% lifetime AML
   risk instead of 0.47%. Halving the time step at fixed $N$ gives 0.15%. The published fit to
   UK incidence holds only at the code's implicit $N\tau=10^4$.
2. **Acquisition ages have an irreducible spread and a scale-dependent offset.** The fitted
   establishment time differs from the true founding time by $-(\log X+\kappa)/s$, with $X$ the
   unobservable establishment luck (s.d. $1.28/s$, 13 years at $s=0.1$/yr) and
   $\kappa=\log(10^5/(2N\tau))$ (Prop. 4). The sign of the offset flips between $N\tau=10^4$ and
   $10^5$.
3. **The reported uncertainties are not confidence intervals.** They are the spread of optimiser
   end points over 25 seeds within ±20% of the best negative log-likelihood. On the public data,
   proper 95% profile-likelihood intervals for reported ages, corrected for the measured
   overdispersion (Pearson dispersion 3.1), have median width 51 years, against 7 years for the
   published ranges. 144 of 156 case clones are wider.
4. **The "first driver by age 20 in two thirds of cases" statement is not supported case by
   case, and it is what the incidence curve predicts anyway.** In the re-fit, 19 of 42 earliest
   clones sit on the constraint "founded at or after birth" (their likelihood prefers a founding
   before birth). Only 3 of 42 have a 95% interval entirely below age 20. Separately, the
   first-driver age law (Prop. 5) gives 74% by age 20 for a woman diagnosed at 72 from UK
   incidence alone, for any time-homogeneous multistage model.
5. **Founders that had swept before sampling are not identified by their trajectories.** When the
   founder's VAF is flat near 0.5 throughout (for example C92_019, C92_023, C92_024, C92_034), its
   fitness and age are set by the birth constraint: the smallest fitness that completes a sweep by
   the first sample when starting at birth.
5b. **Under the paper's own model the estimator does not recover first-driver ages.** On 928
   simulated AML cases (8 annual samples before diagnosis, true clone tree given to the
   estimator), 95% of first drivers had swept before the first sample. Reported ages are on average
   15 years later than the truth (s.d. 24). The estimator reports "first driver by 20" for 57% of
   cases where the truth is 94%. Fitted fitness exceeds the truth by a median factor of 1.5 (1-, 3-
   and 4-driver clones) to 2.2 (2-driver clones).
6. **The geometric fitness ladder does not test the model.** In the simulation the conditioned
   increments sit at the input scales $a_k=s_1r^{k-1}$, so the ladder reproduces the assumed
   $r=2.5$ (theory note §5).
7. **The model cannot produce late evolution, by a factor of about $4\times10^5$** (theory note §6).
8. **The PPV decline with age is mostly the 85-year horizon** of the simulation (§4 below).

Details follow.

## 1. Population scale of the simulation

`sim_fast.select_population` sets $n\leftarrow\mathrm{Poisson}(n\,e^{(F-\bar x)dt})$ with $dt=0.1$
yr and $N=10^5$. This is a Galton-Watson chain with generation time 0.1 yr: the drift variance per
cell per year is $1/dt=10$, which in the Methods' birth-death language is a symmetric division
every $\tau=0.1$ yr, so $N\tau=10^4$ years. The Methods and Discussion state $N\tau=10^5$ and
$\theta=N\tau\mu\approx1$; the code has $\theta=0.1$.

Lifetime risk to 85 in the C++ port (`sim/amlsim.cpp`; it reproduces the published 943/200,200
cases, 642/301/0 patterns and PPVs within Monte Carlo error): $N=10^5$, $dt=0.1$: 0.464%;
$N=5\times10^5$, $dt=0.02$ (same $N\,dt$): 0.536%; $N=10^5$, $dt=0.05$: 0.146%; $N=10^6$, $dt=0.1$
(the stated $N\tau=10^5$): 16.0%.

Either the stated parameter is wrong, or the published simulations do not implement it. The
incidence fit, the pattern fractions and the PPVs of Fig. 5 all come from the $N\tau=10^4$ model.

## 2. Acquisition ages

### 2a. What the fitted age measures (Prop. 4)

Reported "occurrence time" is $t_e=t_0+\log(1+s)/s$, the age at which the fitted deterministic
clone size equals one cell (Figure_4.ipynb uses `establishment_time`). Under the paper's own
stochastic model:
- spread: s.d. of $t_e-t_m$ is $1.28/s$, independent of sequencing depth or sample count;
- offset: at the code's $N\tau=10^4$ reported ages are too early by
  $(\log5-\gamma-\log(1+s))/s$ (9.4 years at $s=0.1$, 2.6 at $s=0.3$); at the stated $N\tau=10^5$
  too late by $(\gamma+\log2+\log(1+s))/s$ (13.7 and 5.1 years).

Applied to the published estimates (`analysis/ages.py`; earliest clone per case, the Fig. 4b
rule): reported $t_e\le20$ in 32 of 42 fitted cases (32/47 = 68%, the paper's "two thirds").
Expected number with true founding $\le20$: 22.3 at $N\tau=10^4$, 28.7 at $3\times10^4$, 32.4 at
$10^5$. Ten of the 42 reported ages are at $t_e\le1.5$, i.e. at the birth constraint.

### 2b. The estimator on simulated patients with known truth

Design (`analysis/sim2refit.py`, `analysis/simfit_eval.py`): the 928 AML cases of the reproduction
run; 8 annual samples at ages dx-8, ..., dx-1; read depth 1500; alt reads binomial on the true
VAF; variants with true VAF $\ge0.001$ at some sample are included, observations with fewer than
3 alt reads are censored. The estimator is the maximum-likelihood fit of §2c, given the true
clone tree of the detected clones, which favours the estimator relative to the real analysis.

Results (`runs/simfit/eval.txt`):
- In 874 of 923 cases the first driver of the AML lineage had VAF > 0.4 at the first sample: it
  had swept before sampling. Its trajectory carries no information on its age or fitness (§3).
- Reported age minus true founding age, first driver: mean +14.9 years, median +7.9, s.d. 23.7.
  Not-swept founders (49): +7.7 ± 18.3. Prop. 4 alone, at the true fitness, predicts -6 years:
  the fitted fitness is biased as well (next items), so the error is dominated by identifiability,
  not by the establishment factor.
- Reported $t_e\le20$: 55% (AML-lineage founders), 57% by the Fig. 4b rule (earliest clone per
  person). True founding $\le20$: 93.5% and 94%. On the real data the estimator gives 68% (32/47).
- Fitted fitness over true fitness, median (IQR): 1-driver 1.52 (0.09-6.44), 2-driver 2.19
  (0.99-5.39), 3-driver 1.51 (1.02-3.08), 4-driver 1.52 (1.07-2.65). Absolute fitness of a clone
  inside a swept founder is the founder's (unidentified) fitness plus an identified difference,
  so the bias propagates up the ladder.
- Pearson dispersion on the simulated read counts: 1.20. On the real data: 3.06. Real VAF series
  fluctuate about 2.5 times more (in variance) than the model plus binomial sequencing implies.

Caveats: the clean MLE differs from the authors' penalised fit; the real design has fewer samples
(median 5), an inferred tree and missing drivers. Each of these makes recovery harder, not easier.

### 2c. Confidence intervals on the real data

`sim/refit.cpp` re-implements the paper's model (Methods eqs. 1-2, the VAF weights of the v16
code for heterozygous SNVs, homozygous JAK2/9p and TET2/4q, and mCA cell fractions) with a
binomial likelihood on read counts, censored likelihood for undetected timepoints and Gaussian
mCA cell fractions. Constraints: $t_0\ge0$, daughter not founded before parent, $0.002\le s\le15$.
It does not impose the authors' penalty that a daughter be fitter than its parent. Starting points
include the published estimates. Profile intervals are computed for $s$ and $t_e$ of every clone,
with the binomial threshold (1.92) and the quasi-binomial threshold $1.92\phi$.

Results (`runs/refit/full/eval.txt`, 42 cases, 33 controls, 222 clones):
- Pooled Pearson dispersion of the detected read counts: $\phi=3.06$ (1 under the binomial
  model). Year-to-year variation is three times the sampling noise the likelihood assumes, so
  binomial intervals are too narrow by a factor $\sqrt{3}$ even before identifiability problems.
- Width of the 95% interval for $t_e$: median 51 years (cases) and 53 years (controls), against
  published optimiser ranges of 7 and 4.4 years.
- The intervals reach birth for 14 of 156 case clones; the fitness interval reaches the upper
  search bound 15/yr for 59 of 156 (fast late clones whose growth is only bounded below).
- Published values lie inside the quasi intervals for 149/156 ages and 131/156 fitnesses. The
  point estimates are compatible with the data; their stated precision is not.
- In 15 of 42 cases the published parameters are worse than the maximum-likelihood fit by more
  than the quasi 95% threshold for one parameter (at the best common time shift for the age
  rounding), in 7 by more than four times that threshold. Their penalties and heuristics (one-
  sided penalties for undetected timepoints, ordering penalties, manual starting values in
  `Clonal_structure_for_each_sample.csv`) move the fit away from the stated likelihood.

### 2d. The earliest clone per case (Fig. 4b)

Maximum-likelihood $t_e\le20$ in 26 of 42; 19 of these sit at the constraint $t_0=0$ ($t_e\approx1$).
The whole quasi 95% interval lies below 20 in 3 cases (C92_023, C92_024, C92_086); it reaches
below 20 in 40. Two of the three (C92_023, C92_024) are founders at VAF about 0.47 in every
sample, whose age is fixed by the birth constraint (§3).

## 3. Founders swept before sampling

For a founder at fixation in every sample the trajectory is flat. Its own fitness and age are not
determined by the data. Through the normalisation they enter only as (i) a lower bound on fitness
— to have swept by the first sample starting no earlier than birth, $e^{s a_1}/s\gtrsim10^5$ at the
first sampling age $a_1$, i.e. $s\gtrsim0.2$/yr for $a_1\approx60$ — and (ii) the difference
$s_{\rm daughter}-s_{\rm founder}$ that governs daughters growing inside it. The optimiser then
lands at the corner $t_0=0$, $s\approx0.2$: in the re-fit C92_022 (0.21, $t_e$ 0.9), C92_023
(0.20, 0.9), C92_024 (0.21, 0.9), and C92_019 (0.18, 0.9). The narrow fitness intervals of some of
these founders (C92_024: 0.208-0.277) come from the birth constraint, not from the data.
Because daughters' absolute fitness equals the founder's plus the identified difference, the
absolute fitness of every clone growing inside a swept founder inherits this indeterminacy.

## 4. Smaller points

- `figure5_pipeline.py` uses observed pattern counts 26/12/9; the paper reports 26/11/10.
  The `fig5e_pipeline.py` docstring quotes OR = 98.8 for Fig. 5f; the paper reports 107.
- Fig. 5c compares with the average of male and female CRUK incidence; the cohort is women
  (female cumulative risk to 85 is lower).
- Fig. 5f counts any later myeloid leukaemia (ICD-10 C92, which includes CML) and hotspot sweeps
  that include JAK2 V617F; the Discussion reads it as AML risk.
- The PPVs of Fig. 5e depend on the simulation horizon (AML before 85). Following the same
  simulated lives to 100 raises the PPV at 70 from 21% to 62% and at 60 from 46% to 80% (theory
  note §7). The decline with age in the main text is mostly this cut-off.

## Files and reproduction

- `sim/amlsim.cpp`: port of `sim_fast.py`. `runs/repro/` (published parameters), `runs/ndt/`
  (scale tests), `runs/horizon/`.
- `analysis/prep_refit.py` (public tables to refit input), `sim/refit.cpp`,
  `analysis/refit_eval.py`; `runs/refit/full/`.
- `analysis/sim2refit.py`, `analysis/simfit_eval.py`; `runs/simfit/`.
- `analysis/ages.py` (Prop. 4 applied to the published estimates).
- The authors' Python code was read but not executed.
