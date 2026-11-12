# Statistics re-check of Watson et al., Nature 2026

Paper: Watson, Poon, MacGregor, Fonseca, Apostolidou, Damaskou, Boniface, Hackett,
Gentry-Maharaj, Menon, Blundell, *Evolutionary dynamics in the years preceding acute myeloid
leukaemia*, Nature, online 2026-11-11, doi:10.1038/s41586-026-11269-x.
Code and public data: github.com/the-blundell-lab/preAML_evolutionary_dynamics, commit 5a53fe7
(`Data_files/`). Controlled-access data (EGA) were not used; ages in the public tables are rounded
down to whole years (relative timing within a person is exact to a month).

Companion: `theory/theory.md` (propositions referred to below as Prop. n).

## Summary

Revised 2026-11-12 after an independent ChatGPT Pro review (`notes/pro_review_01.md`, job
job_f9cc4981, gpt-6-pro); verdicts noted per item.

1. **The stated simulation parameters do not reproduce the published results; the effective ones
   do** (Pro: holds with qualification; re-tuning added after the review). The paper states
   $N\tau=10^5$ (from its refs 38, 39, 41) and $U=10^{-5}$. `sim_fast.py` resamples every 0.1 yr;
   matched to the explicit birth-death rates of ref 39 (Watson 2020, same group) that is $\tau=0.2$ yr
   and $N\tau=2\times10^4$, five times below the stated value. (An earlier version also claimed a
   tenfold gap in "coalescent units"; withdrawn.) At the stated scale with $U=10^{-5}$ unchanged the
   model gives a 6.6% lifetime risk against 0.47% observed. Re-tuned at the stated scale
   ($U=3.8\times10^{-6}$, 200,200 lives) it gives 0.50% and approximately preserves the simulated
   incidence shape, the pattern fractions (66/34/0%, published 69/31/0%), the median driver ages (5.4, 39.8, 61.9, 74.2
   against 6.0, 39.1, 61.1, 74.1) and the conditioned fitness ladder; the sweep PPVs change (40-70:
   24/24, 133/142, 480/783, 770/2984, i.e. 100% [85.75-100], 94%, 61%, 26%, against 89%, 72%, 46%, 21%). Re-tuning the
   fitness scale instead ($s_1=0.1165$) also fits incidence but lowers absolute conditioned fitness
   by 21-28% (adjacent ratios approximately similar), gives more branched cases (57/43/0%) and higher PPVs (60: 82%, 70: 39%). So this is
   an internal inconsistency between the stated and the simulated parameters, which matters for the
   reported PPVs, possibly for the absolute fitness scale, and for anyone re-using $(N\tau, U)$; it
   does not by itself overturn the listed outputs (theory §8; two one-parameter re-tunings only).
   Because the code's drift scale is $N\tau=2N\,dt$, the time step is a model parameter, not a
   refinable numerical step: halving it at fixed $N$ halves $N\tau$ and gives 0.15%.
2. **Fitted establishment times carry an irreducible spread and a scale-dependent offset, for a
   surrogate time** (Pro: the earlier wording failed; repaired). Under the exponential-establishment
   approximation, the own-size time $\tilde t_e=t_0+\log(1+s)/s$ differs from the true founding time
   by $-(\log X+\kappa-\log(1+s))/s$, with s.d. $1.28/s$ for a clone fitted in its exponential phase
   and not selected on its trajectory, and $\kappa\approx\log(10^5/(N\tau))$ in ref-39 units (theory
   Prop. 4). The paper's reported establishment time is a different quantity (first time the
   normalised frequency reaches $10^{-5}$); the two agree only when other clones are negligible.
3. **The published uncertainty ranges have no confidence-interval justification** (Pro: holds with
   qualification). They are the spread of optimiser end points over 25 seeds within ±20% of the best
   negative log-likelihood; seeds do not generate a sampling distribution, and a 20% tolerance does
   not define coverage (it even depends on additive constants in the likelihood). An alternative
   re-fit of the public data (§2c) exposes large uncertainty and instability: approximate profile
   intervals for the own-size age have median width 52 years, and different restarts give very
   different point estimates at numerically equal likelihood. These are not calibrated 95%
   intervals: they target a modified objective and estimand and undercover in simulation.
4. **Individual first-driver ages are weakly identified; the aggregate "by age 20" figure needs its
   own analysis** (Pro: holds with qualification). In the re-fit, 19 of 42 earliest clones have their
   maximum-likelihood founding on the birth constraint and only 3 have an interval entirely below
   20. That shows weak identification of single cases, not that the population proportion is wrong:
   many uncertain cases can still inform a proportion. Separately, Prop. 5 shows that for models
   with constant founder supply and independent delays the incidence curve itself predicts the
   founding-age distribution (74% by 20 for a woman diagnosed at 72 under UK female rates); agreement
   with that is support for the model class, not discriminating support for the paper's mechanism.
5. **Founders whose lineage is near fixation before sampling are weakly identified** (Pro: holds with
   qualification). A first-driver VAF near 0.5 means its mutation-bearing lineage dominates, possibly
   through descendants; the trajectories then constrain fitness differences much better than
   absolute fitness and founding time, and the fit becomes sensitive to constraints (often landing
   on the birth boundary).
5b. **Simulated estimator tests (alternative re-fit, not the authors' full estimator)** (Pro: holds
   with qualification). In simulated AML cases whose founder lineage had not swept at the first
   sample (272 cases), founding-age errors had observed median within 0.2 years of zero and s.d. 18
   years, and the estimator reported "first driver by 20" for 70% against a true 84% (the same with
   the repaired re-fit); with the full simulated case
   mix (95% swept) it reported 57% against 94%. These describe the tested designs; they are not a
   calibration of the paper's 68% (32/47), and the real cohort's mix differs from the model's.
6. **Fitness ladder** (Pro: holds with qualification). In the non-swept simulated stratum the
   estimator recovers clone fitness with median ratio 1.0-1.2, though with wide spread for 1- and
   4-driver clones. In one experiment with no built-in escalation ($r=1$) and all founders swept,
   the fitted ladder was strongly distorted (19.8-fold fitted against 3.9-fold true). That is a
   sensitivity warning for swept founders, not a general floor: under the published model the swept
   founders' fitted median was 0.23 against 0.18 true. In the simulation the geometric escalation is
   imposed by the input distributions (scales $s_1r^{k-1}$), so the model's agreement with Fig. 4c
   reflects that assumption.
7. **Late evolution: a gap remains under comparable classification, but it is far smaller than first
   stated** (Pro: holds with qualification; observation-matched check added 2026-11-12, §2e). Against
   the simulation's own classifier (complete lineage trajectories) the model gives 0 late cases in
   928, about 27 times below the observed 10/47 (21%, 95% interval 10.7-35.7%). That comparison is not
   like for like. Applying the real observation process to the simulated cases (each case sampled on
   a real case's schedule, detection at VAF 0.1%, the paper's rule: no detection, or first detection
   within 2 years after an earlier negative sample) gives 0 late cases with a complete panel:
   sampling and detection limits alone do not create late cases. Drivers outside the panel do. With
   the off-panel rate calibrated to the real data (12-15% of drivers, matching the 40-48% of cases in
   which the authors inferred a missing driver or found none), the model gives 0.15-0.23% late cases
   if the pre-leukaemic clone is the first driver's lineage, or 2.6-4.2% under the authors' apparent
   rule that an earlier clone can be background (their case C92_016, listed "linear/late", was
   detected 6 years before diagnosis). The gap is then roughly 3- to 8-fold under the authors' rule,
   not 27-fold. It disappears only if about 30% of drivers are off-panel, which would imply about 70%
   of cases with a missing driver against 40% observed. Excluding C92_016 the observed fraction is
   9/47 (19%, 9.1-33.3%).
8. **The PPV decline with age stands** (Pro: holds). With explicit competing mortality the PPVs are
   close to the published ones (§4); an earlier draft calling the decline an artefact of the age-85
   cut-off was wrong. This shows robustness within the model and scenarios, not calibration of
   real-world PPVs.
9. **Possible issues in the authors' inference code** (found by the Pro review; effects not
   assessed). In v16 `likelihood_clones`, `establishment_age` is not reset between clones, so a clone
   whose frequency never reaches $10^{-5}$ would inherit the previous clone's value (whether any
   clone takes that path needs a code-path check); reported ages use `np.linspace(0, T*dt, T)` while
   trajectories use `i*dt`, which differ by at most one step $h$; in `sim_fast.py` diagnosis is
   recorded as `t*dt` after the update at index $t$ (effect not checked).

Details follow.

## 1. Population scale of the simulation

`sim_fast.select_population` sets $n\leftarrow\mathrm{Poisson}(n\,e^{(F-\bar x)dt})$ with $dt=0.1$
yr and $N=10^5$. This is a Galton-Watson chain with generation time 0.1 yr: the drift variance per
cell per year is $1/dt=10$ and a new clone of fitness $s$ survives with probability $\approx2s\,dt$.
In the birth-death model of ref 39 (Watson et al. 2020, Supplementary Note 2 of bioRxiv 569566:
self-renewal and differentiation each at rate $1/\tau$ per cell) the variance is $2/\tau$ per year
and survival $\approx s\tau$, so the code has $\tau=0.2$ yr and $N\tau=2\times10^4$. In the
phylodynamic (coalescent) units of refs 38 and 41 (effective population size times generation
time; pairwise coalescence rate $1/(N\,dt)$ per year for the Poisson update) the code has
$N\tau=10^4$. The Methods and Discussion state $N\tau=10^5$ and $\theta=N\tau\mu\approx1$; the code
has $N\tau\mu=0.2$ (ref-39 units) or 0.1 (coalescent units).

Lifetime risk to 85 in the C++ port (`sim/amlsim.cpp`; it reproduces the published 943/200,200
cases, 642/301/0 patterns and PPVs within Monte Carlo error): $N=10^5$, $dt=0.1$: 0.464%;
$N=5\times10^5$, $dt=0.02$ (same $N\,dt$): 0.536%; $N=10^5$, $dt=0.05$: 0.146%; $N=5\times10^5$,
$dt=0.1$ (the stated $N\tau=10^5$ in ref-39 units): 6.6% (20,000 lives); $N=10^6$, $dt=0.1$ (the
stated value in coalescent units): 16.0% (4,000 lives).

Re-tuned at the stated scale ($N=5\times10^5$, $dt=0.1$, $U=3.8\times10^{-6}$; `runs/retune/`,
`runs/retune/compare_U.txt`), the model reproduces the published incidence, patterns, driver timing
and conditioned fitness; only the sweep PPVs shift (summary item 1). Either the stated parameters
are mislabelled, or the published simulations do not implement them. The
incidence fit, the pattern fractions and the PPVs of Fig. 5 all come from the model at the code's
implicit scale. In this code the time step sets the drift scale ($N\tau=2N\,dt$), so it is a model
parameter rather than a refinable numerical step.

## 2. Acquisition ages

### 2a. What the fitted age measures (Prop. 4)

Reported "occurrence time" is $t_e=t_0+\log(1+s)/s$, the age at which the fitted deterministic
clone size equals one cell (Figure_4.ipynb uses `establishment_time`). Under the paper's own
stochastic model, with the exponential-establishment approximation A1, for a clone of known
fitness $s$ fitted in its exponential phase and not selected on its trajectory:
- spread: s.d. of $t_e-t_m$ is $1.28/s$, independent of sequencing depth or sample count;
- offset: at the code's scale reported ages are too early by $(\log5-\gamma-\log(1+s))/s$ (9.4
  years at $s=0.1$, 2.6 at $s=0.3$); at the stated $N\tau=10^5$ too late by $(\gamma+\log(1+s))/s$
  in ref-39 units (6.7 and 2.8 years) or $(\gamma+\log2+\log(1+s))/s$ in coalescent units (13.7
  and 5.1 years).

Real cases are selected (they are AML cases, observed in a window) and $s$ is estimated, so they do
not inherit this error law; the simulation tests in §2b and §2b2 show both offset and spread
changing under that selection.

Applied to the published estimates (`analysis/ages.py`; earliest clone per case, the Fig. 4b
rule): reported $t_e\le20$ in 32 of 42 fitted cases (32/47 = 68%, the paper's "two thirds").
Expected number with true founding $\le20$ under the simplified posterior of theory §2 (constant
founding prior, fixed fitted $s$): 22.3 at the code's scale, 30.7 at the stated $N\tau=10^5$ in
ref-39 units. These rest on that simplified prior and observation model and are not calibrated
estimates for the real cohort. Ten of the 42 reported ages are at $t_e\le1.5$, i.e. at the birth constraint.

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

The simulated fits in §2b and §3b used refit as it was before the 2026-11-12 repair (for these
inputs the only difference is the former exponent cap, whose effect can propagate through the
common denominator); they were not rerun. The matched-stratum test (§2b2) and the coverage check
(§2c) were rerun with the repaired refit, with unchanged results (`runs/matched/ns5/`,
`runs/coverage5/`).

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
- Pearson dispersion on the simulated read counts: 1.20. On the real data: 3.14. Real VAF series
  show about 2.5 times more residual variance than the fitted model plus binomial sequencing
  implies; this can come from biological variation, measurement differences or mean-model
  misspecification, and does not by itself measure biological year-to-year variance.

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
- Pearson dispersion 1.14 (real data 3.14).
- Rerun with the repaired re-fit (`runs/matched/ns5/eval.txt`): median error -0.2 years, mean +4.7,
  s.d. 17.9; Fig. 4b rule reported 70.2% (true 84.2%); fitness ratios 1.18, 1.04, 1.02, 1.10. The
  repair does not change these results.

A standardised exercise: weighting the simulated strata by the real proportions (28 not swept,
11 swept among the 39 cases with an observed founder variant) gives "by 20" reported 66.1%, true
87.2%. This is arithmetic on simulated strata with a mix the model itself does not produce; it is
not a calibration of the paper's 68% (32/47, or 32/42 = 76% among fitted cases), and the observed
fraction among the same 39 cases is not computed here.

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

`sim/refit.cpp` re-implements the core of the paper's model: Methods eqs. 1-2; the VAF weights of
the v16 `likelihood_clones` in its branch order (a clone containing JAK2 V617F takes the JAK2
branch, with JAK2 weight 1 only alongside chr9p CN-LOH, and never the chr4q/TET2 branch; otherwise
chr4q CN-LOH gives weight 1 to TET2 variants; mCAs weight 1; others 1/2); the v16 `Upper_freq`
selection (0.4 for C92_034, 1 otherwise), implemented as excluding detected observations with
alt/depth (or mCA fraction) at or above `Upper_freq`; and a binomial likelihood on read counts.
This is not a literal copy of v16's condition $0<\text{freq}<$ `Upper_freq` on the raw VAF field
(refit does not exclude zero values and filters on alt/depth), but on the public inputs the two
select the same observations: there are no detected rows with zero alternate reads and no row where
the alt/depth decision differs from the VAF-field decision, including C92_034 (independent audit by
yang2, 2026-11-12). The Jander code review (2026-11-12) found the earlier weight rule out of
branch order; no real clone contains both JAK2 V617F and chr4q CN-LOH, so that change alters no
fit, whereas adding the omitted `Upper_freq` filter changes C92_034.

Deliberate differences from the authors' objective: their penalties are omitted (refit has its
own soft constraints: $t_0\ge0$, daughter not founded before parent, $0.002\le s\le15$; no
requirement that a daughter be fitter than its parent); undetected timepoints enter as a censored
$P(K\le2)$ term; mCA cell fractions are Gaussian with the phased-caller s.d. (v16: binomial at
depth 1800); sample times are exact rather than on a 0.1-yr grid; probabilities are clipped to
$[10^{-12},1-10^{-12}]$ (v16 returns a fixed penalty at 0 or 1); read depths are the public
`total_depth` (v16 reads the DCS depth of the all-positions files, under controlled access). The
reported $t_e$ is where the clone's own deterministic size is one cell; v16's `establishment_time`
is where its normalised frequency first reaches $10^{-5}$, which is later when other clones are
large (equal for clones founded while the population is about $10^5$). Clone sizes are normalised
in log space, with no cap on the exponent. Starting points include the published estimates.
Profile intervals are computed for $s$ and $t_e$ of every clone, with the binomial threshold
(1.92) and the quasi-binomial threshold $1.92\phi$.

Results (`runs/refit/full5/eval.txt`, 42 cases, 33 controls, 222 clones; repaired code of
2026-11-12, replacing `runs/refit/full/`):
- Pooled Pearson dispersion of the detected read counts: $\phi=3.14$ (1 under the binomial
  model): residual variation is three times what the binomial likelihood assumes (from biology,
  measurement or mean-model misspecification; not separable here). Under
  regular, local conditions this would make binomial intervals too narrow by about $\sqrt{3}$;
  near boundaries and under weak identification the factor is not reliable.
- Width of the 95% interval for $t_e$: median 52 years (cases) and 53 years (controls), against
  published optimiser ranges of 7 and 4.4 years.
- The age interval reaches the birth constraint ($t_0=0$, within 0.3 yr) for 74 of 156 case
  clones, and for 31 of the 42 earliest clones; 19 of those 42 have their maximum-likelihood fit
  on it. (An earlier draft reported 14 of 156 from a flag that could not detect this constraint.)
  The fitness interval reaches the upper search bound 15/yr for 57 of 156 (fast late clones whose
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
  so these results do not imply an ordering of widths. Rerun with the repaired re-fit
  (`runs/coverage5/eval.txt`): fitness 67%, 93%, 91%, 97%; founding time 87%, 83%, 52%, 59%.
- Flat ridges: rerunning the fit (repaired code, new random restarts) moved the point estimates of
  16 persons other than C92_034 while their maximum log-likelihood agreed with the previous run to
  three decimals (two others improved by 0.07 and 0.89, where the old exponent cap had been
  active). Example: C92_011's fitness 7.5 to 1.5 per year at a likelihood that agrees to 0.001.
  The likelihood is numerically indistinguishable at that precision between such end points (not
  shown to be exactly flat, nor flat along the whole path between them), so the point estimate is
  largely wherever the optimiser stops; the intervals are the meaningful summary.
- Published values lie inside the quasi intervals for 151/156 ages and 137/156 fitnesses. The
  point estimates are compatible with the data; their stated precision is not.
- Joint test of the published parameter vector against the re-fit's optimum (threshold
  $\chi^2_{0.95}(2K)/2\cdot\phi$ for $2K$ parameters, at the best common time shift for the age
  rounding): outside the joint region for 6 of 42 cases and 6 of 33 controls. (An earlier version
  used a one-parameter threshold and reported 15 of 42; that was the wrong reference.) Since the
  re-fit's objective differs from the authors' (penalties, censoring, mCA likelihood, depths), these
  differences cannot be attributed to the authors' optimisation alone.

### 2d. The earliest clone per case (Fig. 4b)

Maximum-likelihood $t_e\le20$ in 25 of 42; 19 of these sit at the constraint $t_0=0$ ($t_e\approx1$).
The whole quasi 95% interval lies below 20 in 3 cases (C92_023, C92_024, C92_086); it reaches
below 20 in 39. A birth-boundary maximum is not positive evidence of founding at birth, but a lower
constraint does not by itself explain an upper interval end below 20 either: for C92_023 and
C92_024 daughter trajectories and other constraints may carry real information.

### 2e. Late evolution under comparable classification

The paper's real-data rule (Methods): late if the putative pre-leukaemic clone is detectable only
within 1-2 years before diagnosis, or no driver is detected. The ten real late cases (authors'
`Figure_2.ipynb` list): five with no detection at all (C92_015, C92_043, C92_044, C92_048, C92_063;
for three of them the last sample was 3.5-5.5 years before diagnosis), four first detected within 2
years after earlier negative samples (C92_047, C92_052, C92_055, C92_072), and C92_016, detected at
every sample from 6.25 years and listed as "linear/late" (its late NPM1 clone is the presumed
pre-leukaemic clone, with an earlier DNMT3A clone treated as background). The two real cases sampled
only in the final 2 years (C92_065, C92_066) were not called late, so a late detection evidently
requires an earlier negative sample.

Simulated version (`analysis/late_obs.py`, `runs/late_obs2.txt`, `runs/late_obs3.txt`): each of the
928 simulated AML cases receives the sample schedule of a randomly drawn real case (40 draws per
case; exact months where the public calls file has them, else timing-group midpoints); a trunk
mutation is detected when its lineage VAF is at least 0.001; each trunk driver is off-panel
independently with probability $p$; late = no detection, or first detection within 2 years after an
earlier negative sample. Rule A takes the first visible trunk driver's lineage as the pre-leukaemic
clone, rule B the second (an earlier visible driver treated as background, as for C92_016).

| off-panel probability $p$ | P(at least 1 of 4 drivers off-panel) | late, rule A | late, rule B |
|---|---|---|---|
| 0 | 0 | 0.00% | 0.00% |
| 0.12 | 0.40 | 0.15% | 2.56% (no detection 1.80%) |
| 0.15 | 0.48 | 0.23% | 4.17% (no detection 2.98%) |
| 0.3 | 0.76 | 2.0% | 16.5% (no detection 13.0%) |
| 0.5 | 0.94 | 9.8% | 43% |

Observed: 10/47 = 21.3% (10.7-35.7%); no detection 5/47 (3.5-23.1%); late detection after a
negative sample 4/47 (2.4-20.4%). The real data show a missing driver (inferred X/Y/Z, 14 cases) or no
driver (5 cases) in 19 of 47 cases (40%), which calibrates $p\approx0.12$-$0.15$ under the
independence assumption.

Conclusions. (i) The simulation's own classifier and the paper's real-data rule are not comparable,
so the 27-fold figure overstates the gap. (ii) Sparse sampling and the detection limit do not by
themselves create late cases in the model. (iii) At the calibrated off-panel rate the model's late
fraction is 0.15-4.2% depending on the definition, still below the observed interval; matching it
would need about 30% of drivers off-panel, inconsistent with the 40% of cases showing a missing
driver. (iv) Independent per-driver panel misses are a crude model: a subset of AMLs with largely
invisible routes (the paper's own suggestion: off-panel drivers or single ultra-rare events) lies
outside the published model and could close the gap. Limits: depth noise and calling thresholds are
not modelled beyond the 0.1% VAF cut; draws reuse the 928 cases, so the effective sample is 928.

## 3. Founders swept before sampling

When a founder's variant is near fixation in every sample, its mutation-bearing lineage dominates,
possibly through descendants rather than the one-driver genotype itself. The trajectory then
constrains fitness differences much better than absolute fitness and founding time. If the
one-driver genotype itself swept, its own parameters enter mainly as (i) a lower bound on fitness
— to have swept by the first sample starting no earlier than birth, $e^{s a_1}/s\gtrsim10^5$ at the
first sampling age $a_1$, i.e. $s\gtrsim0.2$/yr for $a_1\approx60$ — and (ii) the difference
$s_{\rm daughter}-s_{\rm founder}$ that governs daughters growing inside it. In the re-fit the
optimiser often lands near the corner $t_0=0$, $s\approx0.2$ (it need not; nearly equivalent
points exist), for example C92_022 (0.21, $t_e$ 0.9), C92_023
(0.20, 0.9), C92_024 (0.21, 0.9), and C92_019 (0.18, 0.9). For some of these founders the
fitness interval is narrow (C92_024: 0.208-0.277); its lower end is set by the birth constraint,
its upper end by the data through the daughters and the normalisation.
To the extent daughters' absolute fitness is the founder's plus a better-identified difference,
their absolute fitness inherits this weak identification.

## 3b. Sensitivity of the fitted ladder to swept founders (an r = 1 experiment)

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
sample. Their fitted fitness (median 0.13, a quarter of the truth) sits near the birth-constraint
boundary in this experiment (§3). One of the 400 fits failed (S062166, chi-square per observation about 6 × 10⁵); excluding it
changes no median by more than 0.004. With the published $r=2.5$ (§2b) the estimator's ladder
(0.23, 1.24, 2.33, 5.85) is also steeper at the first step than the truth (0.18, 0.57, 1.54, 3.93).

In this experiment swept founders produced a steep fitted first step without built-in
escalation. It is a sensitivity warning, not a general floor: under the published model the swept
founders' fitted median was 0.23 against 0.18 true (§2b). No non-swept founder occurred among these
400 $r=1$ cases, so that stratum was not tested under $r=1$; in the published model the estimator
recovers fitness well for non-swept founders (§2b2). The concern applies to real cases with swept
founder lineages (11 of 42 by first-sample VAF) and the 19 of 42 earliest clones fitted at the
birth constraint (6 in both groups); a boundary estimate is not by itself a reason to discard it.

## 4. Smaller points

- `figure5_pipeline.py` uses observed pattern counts 26/12/9; the paper reports 26/11/10.
  The `fig5e_pipeline.py` docstring quotes OR = 98.8 for Fig. 5f; the paper reports 107.
- Fig. 5c compares with the average of male and female CRUK incidence while the cohort is women
  (female cumulative risk to 85 is lower); whether that is an error depends on the intended target
  population. This and the other points in this section are audit questions; the relevant pipelines
  and clinical definitions are not reproduced here.
- Fig. 5f counts any later myeloid leukaemia (ICD-10 C92, which includes CML) and hotspot sweeps
  that include JAK2 V617F; the Discussion reads it as AML risk.
- The PPVs of Fig. 5e count AML before 85 in a model without other causes of death. Following the
  same simulated lives to 100 raises the PPV at 70 from 21% to 62%, but that is the figure for a
  woman who cannot die of anything else. Weighting each AML by survival from the screening age to
  diagnosis, $S(T)/S(a)$, with other-cause death independent of the clones (Gompertz, life
  expectancy 16.5 years at 70, an assumed scenario), the PPVs are 81%, 64%, 44%, 28% at 40-70,
  close to the published 89%, 72%, 46%, 21% (theory note §7). AML after 100 is not counted (below
  2.6 points in each scenario). Under this scenario the decline with age stands.

## Files and reproduction

- `sim/amlsim.cpp`: port of `sim_fast.py`. `runs/repro/` (published parameters), `runs/ndt/`
  (scale tests), `runs/horizon/`.
- `analysis/prep_refit.py` (public tables to refit input), `sim/refit.cpp`,
  `analysis/refit_eval.py`; `runs/refit/full/`.
- `analysis/sim2refit.py`, `analysis/simfit_eval.py`; `runs/simfit/`.
- `analysis/ages.py` (Prop. 4 applied to the published estimates).
- The authors' Python code was read but not executed.
