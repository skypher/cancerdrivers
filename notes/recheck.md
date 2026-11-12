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
2. **Acquisition ages have an irreducible spread and a scale-dependent offset.** Under the
   exponential-establishment approximation, the reported age $t_e$ differs from the true founding
   time by $-(\log X+\kappa-\log(1+s))/s$, with $X$ the unobservable establishment luck (s.d.
   $1.28/s$ for an unselected clone, 13 years at $s=0.1$/yr) and $\kappa=\log(10^5/(2N\tau))$
   (Prop. 4). The sign of the offset flips between $N\tau=10^4$ and $10^5$. In cases selected on
   their trajectories both offset and spread change.
3. **The reported uncertainties are not confidence intervals.** They are the spread of optimiser
   end points over 25 seeds within ±20% of the best negative log-likelihood. On the public data,
   approximate 95% profile-likelihood intervals for reported ages, with the quasi-binomial
   correction for the measured overdispersion (Pearson dispersion 3.1), have median width 51 years,
   against 7 years for the published ranges; 144 of 156 case clones are wider. These intervals are
   themselves approximate: on simulated data with known truth their fitness intervals cover the
   true fitness for only 68% of single-driver clones (nominal 95%), and for 93%, 91% and 97% of
   2-, 3- and 4-driver clones (§2c). Read as intervals for the true founding time, a different estimand from the fitted
   $t_e$, they cover 52-86%; that is a diagnostic of using them for that purpose, not a coverage
   audit of $t_e$.
4. **The "first driver by age 20 in two thirds of cases" statement is not supported case by
   case, and the incidence curve predicts it anyway.** In the re-fit, the maximum-likelihood
   founding of 19 of 42 earliest clones sits on the constraint "founded at or after birth"; for 33
   of 42 the approximate 95% interval reaches that constraint. Only 3 of 42 have an interval
   entirely below age 20. Separately, the first-driver age law (Prop. 5) gives 74% by age 20 for a
   woman diagnosed at 72 from UK incidence alone, for any model satisfying its assumptions
   (constant supply of established founders, independent identically distributed delays).
5. **Founders that had swept before sampling are weakly identified by their trajectories.** When
   the founder's VAF is flat near 0.5 throughout (for example C92_019, C92_023, C92_024, C92_034),
   the data say little about its fitness and age; the fit lands on the birth constraint, at about
   the smallest fitness that completes a sweep by the first sample when starting at birth.
5b. **Under the paper's own model, the estimator's first-driver ages show an observed median
   error of +0.2 years in the matched simulation, with wide errors, and it misses early
   founding.** Simulated AML cases matched to the real
   stratum (founder not swept at the first sample, 272 cases): founding-age error median +0.2
   years, s.d. 18; reported "first driver by 20" 70% against a true 84%. Swept founders (95% of
   simulated cases, 11 of 42 real ones; 874 simulated): error mean +15 years, reported 56%
   against a true 95%.
   Reweighted to the real mix, the estimator would report about 66% (the paper finds 68%) where
   the model's truth is about 87%. So the paper's 68% is what its model plus its estimator would
   produce; in the model it understates early founding, and single-case ages carry an 18-year s.d.
   The median error +0.2 years is observed, not a proof of unbiasedness. These are statements
   about the paper's model, not a calibration of the real data (§2b, §2b2).
6. **In simulation the fitness ladder is mostly recovered, except for swept founders.** In the
   matched simulated stratum the estimator recovers clone fitness to within a median factor of
   1.0-1.2 at every driver number. This is specific to the paper's model and does not validate the
   real ladder. For a founder that swept before sampling it assigns roughly the smallest
   fitness that could sweep from birth (about 0.13 per year) whatever the truth: from an $r=1$
   truth (additive ladder, 3.9-fold rise) such cases give a 19.8-fold fitted rise. In the real data
   this concerns the 11 of 42 cases with swept founders (first-sample VAF > 0.4) and, overlapping
   with them in 6 cases, the 19 earliest clones fitted at the birth constraint (fitted fitness
   0.05-0.21 per year); those founders' fitted values should be discounted. Separately, in the published model the conditioned
   increments sit at the input scales $a_k=s_1r^{k-1}$, so the model's agreement with Fig. 4c
   reproduces the assumed $r=2.5$ rather than testing it (§3b, theory note §5).
7. **The model essentially never produces late evolution.** No late case among 928 simulated AML
   cases (published: 0 of 943): two-sided 95% Clopper-Pearson upper end 0.40% of AML cases. The
   observed fraction is 10/47, two-sided 95% interval 10.7-35.7%. The separate interval endpoints
   differ by about 27-fold; a jointly calibrated conservative statement (Bonferroni, each
   interval at 97.5%) gives a ratio of at least 20 with 95% confidence. This assumes independent
   cases and comparable classification of simulated and real cases. A first-moment heuristic puts
   the gap near $4\times10^5$; that figure is an estimate, not a bound (theory §6).
8. **The PPV decline with age stands.** The 85-year cut-off is a crude proxy for death from other
   causes; under an explicit competing-mortality scenario the PPVs are close to the published ones
   (§4 below). An earlier draft of this note called the decline an artefact of the cut-off; that
   was wrong.

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
stochastic model, with the exponential-establishment approximation A1, for a clone of known
fitness $s$ fitted in its exponential phase and not selected on its trajectory:
- spread: s.d. of $t_e-t_m$ is $1.28/s$, independent of sequencing depth or sample count;
- offset: at the code's $N\tau=10^4$ reported ages are too early by
  $(\log5-\gamma-\log(1+s))/s$ (9.4 years at $s=0.1$, 2.6 at $s=0.3$); at the stated $N\tau=10^5$
  too late by $(\gamma+\log2+\log(1+s))/s$ (13.7 and 5.1 years).

Real cases are selected (they are AML cases, observed in a window) and $s$ is estimated, so they do
not inherit this error law; the simulation tests in §2b and §2b2 show both offset and spread
changing under that selection.

Applied to the published estimates (`analysis/ages.py`; earliest clone per case, the Fig. 4b
rule): reported $t_e\le20$ in 32 of 42 fitted cases (32/47 = 68%, the paper's "two thirds").
Expected number with true founding $\le20$: 22.3 at $N\tau=10^4$, 28.7 at $3\times10^4$, 32.4 at
$10^5$. Ten of the 42 reported ages are at $t_e\le1.5$, i.e. at the birth constraint.

### 2b. The estimator on simulated patients with known truth

Design (`analysis/sim2refit.py`, `analysis/simfit_eval.py`): the 928 AML cases of the reproduction
run; 8 annual samples at ages dx-8, ..., dx-1; read depth 1500; alt reads binomial on the true
VAF; variants with true VAF $\ge0.001$ at some sample are included, observations with fewer than
3 alt reads are censored. The estimator is the maximum-likelihood fit of §2c, given the true
clone tree of the detected clones. That helps the estimator. The case mix does not: in 95% of
these simulated cases the founder had swept before the first sample, against 11 of 42 real cases
(earliest clone's variant at VAF > 0.4 in the first sample). The pooled numbers below therefore
measure mainly the failure on swept founders. Results for the stratum that dominates the real data
(founder not swept) are in §2b2.

Results (`runs/simfit/eval.txt`):
- In 874 of 923 cases the first driver of the AML lineage had VAF > 0.4 at the first sample: it
  had swept before sampling. Its trajectory carries little information on its age or fitness (§3).
- Reported age minus true founding age, first driver: mean +14.9 years, median +7.9, s.d. 23.7.
  Not-swept founders (49): +7.7 ± 18.3. Prop. 4 alone, at the true fitness, predicts -6 years:
  the fitted fitness is biased as well (next items), so the error is dominated by identifiability,
  not by the establishment factor.
- Reported $t_e\le20$: 55% (AML-lineage founders), 57% by the Fig. 4b rule (earliest clone per
  person). True founding $\le20$: 93.5% and 94%. On the real data the estimator gives 68% (32/47).
- Swept stratum alone (874 cases): founding-age error mean +15.4 years, median +8.0, s.d. 23.9;
  Fig. 4b rule reported 56.3%, true 94.7%.
- Fitted fitness over true fitness, median (IQR): 1-driver 1.52 (0.09-6.44), 2-driver 2.19
  (0.99-5.39), 3-driver 1.51 (1.02-3.08), 4-driver 1.52 (1.07-2.65). Absolute fitness of a clone
  inside a swept founder is the founder's (unidentified) fitness plus an identified difference,
  so the bias propagates up the ladder.
- Pearson dispersion on the simulated read counts: 1.20. On the real data: 3.06. Real VAF series
  fluctuate about 2.5 times more (in variance) than the model plus binomial sequencing implies.

Caveats: the clean MLE differs from the authors' penalised fit; the real design has fewer samples
(median 5), an inferred tree and missing drivers. Each of these makes recovery harder. The case mix
(95% swept founders, against 26% in the real data) makes the pooled simulation numbers worse than
what the real cases would show; see §2b2.

### 2b2. Matched stratum: founders not swept at the first sample

Design as in §2b, on a larger run of the published model (1.2 million lives, 5606 AML cases,
`runs/matched/`), keeping only the 272 cases whose AML founder had VAF <= 0.4 at the first sample
(`sim2refit.py --notswept`).

Results (`runs/matched/ns/eval.txt`):
- Founding age of the AML-lineage founder, reported minus true: median +0.2 years, mean +5.0, s.d.
  18.2.
- Reported $t_e\le20$: 62% (AML-lineage founders) and 70% (Fig. 4b rule) against true 80% and 84%.
- Fitted over true fitness, median (IQR): 1-driver 1.20 (0.95-2.00), 2-driver 1.04 (0.96-1.22),
  3-driver 1.02 (0.99-1.07), 4-driver 1.10 (0.98-1.67). Medians fitted/true: 0.165/0.128,
  0.468/0.424, 1.449/1.368, 4.283/3.976.
- Pearson dispersion 1.14 (real data 3.06).

Mixture at the real proportions (28 not swept, 11 swept among the 39 cases with an observed
founder variant), using the swept stratum of §2b (Fig. 4b rule): estimator reports "by 20" for
66.1%, truth 87.2%. The paper's 68% is in line with what its model would produce through its
estimator. This is arithmetic on simulated strata, not a calibration of the real data.

Prop. 4's offset ($-9.6$ years at these fitness values) does not show up here because Prop. 4 is
the average over the establishment factor $X$ for a clone with a given founding time. Selecting
AML cases whose founder is still small eight years before diagnosis selects small $X$, which
moves the fitted time later, and changes the spread too.

### 2b3. A mismatch between the model and the data

In the published model 95% of AML founders have swept (VAF > 0.4) by eight years before
diagnosis; the paper itself reports that about 99% of simulated AMLs had a full sweep more than
5 years before diagnosis. In the real data the earliest clone had swept by the first sample in 11
of 42 cases, and the paper counts sweeps in "at least a third" of the 47. The model produces early
sweeps far more often than the data show.

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
  model). Year-to-year variation is three times the sampling noise the likelihood assumes. Under
  regular, local conditions this would make binomial intervals too narrow by about $\sqrt{3}$;
  near boundaries and under weak identification the factor is not reliable.
- Width of the 95% interval for $t_e$: median 51 years (cases) and 53 years (controls), against
  published optimiser ranges of 7 and 4.4 years.
- The age interval reaches the birth constraint ($t_0=0$, within 0.3 yr) for 77 of 156 case
  clones, and for 33 of the 42 earliest clones; 19 of those 42 have their maximum-likelihood fit
  on it. (An earlier draft reported 14 of 156 from a flag that could not detect this constraint.)
  The fitness interval reaches the upper search bound 15/yr for 59 of 156 (fast late clones whose
  growth is only bounded below).
- Coverage of these intervals on simulated cases with known truth (`runs/coverage/eval.txt`; 150
  cases of the matched stratum of §2b2, the paper's model, nominal 95%, simulated dispersion 1.16):
  fitness intervals contain the true clone fitness for 68% of single-driver clones and 93%, 91%,
  97% of 2-, 3-, 4-driver clones: undercoverage for single-driver fitness in this design. As a
  diagnostic of reading the $t_e$ intervals as intervals for the true founding time (a different
  estimand from the fitted deterministic $t_e$), they contain it for 86%, 83%, 52%, 59% of 1- to
  4-driver clones (median widths 40, 21, 1.5, 2.9 years); for fast clones the interval is narrower
  than the gap between $t_e$ and the founding time (Prop. 4). This is not a coverage audit of the
  $t_e$ parameter itself. A calibrated interval could differ in centre and shape as well as width,
  so these results do not imply an ordering of widths.
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

For a founder at fixation in every sample the trajectory is flat. Its own fitness and age are only
weakly determined by the data. Through the normalisation they enter mainly as (i) a lower bound on fitness
— to have swept by the first sample starting no earlier than birth, $e^{s a_1}/s\gtrsim10^5$ at the
first sampling age $a_1$, i.e. $s\gtrsim0.2$/yr for $a_1\approx60$ — and (ii) the difference
$s_{\rm daughter}-s_{\rm founder}$ that governs daughters growing inside it. The optimiser then
lands at the corner $t_0=0$, $s\approx0.2$: in the re-fit C92_022 (0.21, $t_e$ 0.9), C92_023
(0.20, 0.9), C92_024 (0.21, 0.9), and C92_019 (0.18, 0.9). For some of these founders the
fitness interval is narrow (C92_024: 0.208-0.277); its lower end is set by the birth constraint,
its upper end by the data through the daughters and the normalisation.
Because daughters' absolute fitness equals the founder's plus the identified difference, the
absolute fitness of every clone growing inside a swept founder inherits this indeterminacy.

## 3b. Does the estimator manufacture the fitness ladder?

Test (`runs/null/`): simulate the model with $r=1$ (every driver's increment from the same
distribution, scale $s_1=0.5$, chosen so the lifetime risk is 0.51%), take 400 AML cases, and run
the same estimator and design as in §2b.

| drivers in clone | true conditioned fitness (median) | fitted (median) |
|---|---|---|
| 1 | 0.55 | 0.13 |
| 2 | 1.05 | 0.58 |
| 3 | 1.60 | 1.94 |
| 4 | 2.15 | 2.48 |

The truth is an additive ladder (3.9-fold from 1 to 4 drivers). The estimator reports a 19.8-fold
rise, steeper than geometric over the first two steps. All 400 founders had VAF > 0.4 at the first
sample. Their fitted fitness (median 0.13, a quarter of the truth) is the birth-constraint floor of
§3. One of the 400 fits failed (S062166, chi-square per observation about 6 × 10⁵); excluding it
changes no median by more than 0.004. With the published $r=2.5$ (§2b) the estimator's ladder
(0.23, 1.24, 2.33, 5.85) is also steeper at the first step than the truth (0.18, 0.57, 1.54, 3.93).

So swept founders produce a steep first step whether or not drivers escalate. The test cannot be
repeated on not-swept founders under $r=1$ (none of the 400 was), but in the published model the
estimator recovers fitness well for not-swept founders (§2b2). The bias therefore concerns the
real cases with swept founders (11 of 42 by first-sample VAF) and the 19 of 42 earliest clones
fitted at the birth constraint (§2d; fitted fitness 0.05-0.21; 6 cases are in both groups). This
is a statement about the paper's model in simulation; it does not validate the real ladder.

## 4. Smaller points

- `figure5_pipeline.py` uses observed pattern counts 26/12/9; the paper reports 26/11/10.
  The `fig5e_pipeline.py` docstring quotes OR = 98.8 for Fig. 5f; the paper reports 107.
- Fig. 5c compares with the average of male and female CRUK incidence; the cohort is women
  (female cumulative risk to 85 is lower).
- Fig. 5f counts any later myeloid leukaemia (ICD-10 C92, which includes CML) and hotspot sweeps
  that include JAK2 V617F; the Discussion reads it as AML risk.
- The PPVs of Fig. 5e count AML before 85 in a model without other causes of death. Following the
  same simulated lives to 100 raises the PPV at 70 from 21% to 62%, but that is the figure for a
  woman who cannot die of anything else. Weighting each AML by survival from the screening age to
  diagnosis, $S(T)/S(a)$, with other-cause death independent of the clones (Gompertz, life
  expectancy 16.5 years at 70, an assumed scenario), the PPVs are 81%, 64%, 44%, 28% at 40-70,
  close to the published 89%, 72%, 46%, 21% (theory note §7). AML after 100 is not counted (at most
  4-5 points). Under this scenario the decline with age stands.

## Files and reproduction

- `sim/amlsim.cpp`: port of `sim_fast.py`. `runs/repro/` (published parameters), `runs/ndt/`
  (scale tests), `runs/horizon/`.
- `analysis/prep_refit.py` (public tables to refit input), `sim/refit.cpp`,
  `analysis/refit_eval.py`; `runs/refit/full/`.
- `analysis/sim2refit.py`, `analysis/simfit_eval.py`; `runs/simfit/`.
- `analysis/ages.py` (Prop. 4 applied to the published estimates).
- The authors' Python code was read but not executed.
