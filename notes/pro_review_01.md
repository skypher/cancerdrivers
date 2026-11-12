## Overall assessment

The notes contain several sound criticisms, but also two important implementation/convention errors and some conclusions stronger than their evidence. In particular, **the fivefold population-scale discrepancy is supported; the separate tenfold convention is not established by the supplied excerpts.** The acquisition-time calculation also does **not** describe the inference code’s reported establishment time in general.

I take the simulation and refit numbers as reported results, not as independently verified computations.

## 1. NOTE B: summary items

### 1. Population scale — **HOLDS WITH QUALIFICATION**

Source 1’s neutral Poisson resampling gives variance accumulation \(n/h\) per year. Source 3’s explicit birth–death equations give \(2n/\tau\). Matching these gives
\[
\tau=2h,\qquad N\tau=2Nh=20{,}000.
\]
Thus the **fivefold discrepancy relative to \(100{,}000\)** is supported.

The additional claim that refs 38 and 41 necessarily use \(N\tau=Nh\) is not established. Calling a parameter “coalescent” does not fix its factor-of-two normalization; see Corollary 3 below.

**Correction:** retain the fivefold discrepancy and the reported 6.6% sensitivity result for the specific choice \(N=500{,}000,h=0.1\). Do not present 16% as a second uniquely determined “correct-scale” prediction. Neither run proves that incidence could not be fitted after retuning other parameters.

### 2. Acquisition-age error law — **FAILS as written**

Two corrections are necessary.

First, the summary writes
\[
\kappa=\log\!\frac{10^5}{2N\tau}.
\]
Under ref 39’s explicit convention, the correct expression is
\[
\kappa\simeq\log\!\frac{10^5}{N\tau}.
\]
The factor 2 belongs only in an expression using \(C=Nh\): \(\kappa\simeq\log[10^5/(2C)]\).

Second, Source 2 reports the first **frequency** crossing of \(10^{-5}\), not the time when the clone’s own unnormalized size reaches one. NOTE B §2c acknowledges this, contradicting §2a and the summary’s identification with Proposition 4.

The \(1.28/s\) spread remains a valid approximate result for the stipulated, survival-conditioned, exponential-phase model. It is not a universal uncertainty for the paper’s reported ages, especially for daughters growing inside large ancestors.

### 3. Uncertainty intervals — **HOLDS WITH QUALIFICATION**

The central criticism is correct. Optimizer seeds explore an objective function; they do not generate a sampling distribution. A tolerance of 20% of the best negative log-likelihood does not define 95% coverage and even depends on additive constants in the likelihood.

The reported 52-year profile widths concern a **modified objective and partly different age estimand**, however. They are not established replacements for the authors’ intervals. Scalar quasi-binomial adjustment does not generally cure correlated residuals, structural misspecification, boundaries or weak identification. The reported fitness undercoverage reinforces that limitation.

**Correction:** say the published ranges lack a confidence-interval justification, and that the alternative refit exposes substantial uncertainty and instability—not that calibrated 95% intervals have been demonstrated to be 52 years wide.

### 4. “First driver by 20” — **HOLDS WITH QUALIFICATION**

The boundary pile-up and broad profiles support weak **individual-case identification under the refit**.

But “only three individual intervals lie entirely below 20” does not invalidate a population-level proportion. Many individually uncertain classifications can collectively inform a proportion. Conversely, a birth-boundary MLE alone is not positive evidence of acquisition at birth.

Proposition 5 gives a conditional model-class prediction, not a model-free inference from incidence. Its prediction at diagnosis age 72 is also not automatically the prediction averaged over the study’s diagnosis ages.

**Correction:** distinguish individual identification, estimation of an aggregate proportion, and the incidence-conditioned prediction. Do not treat one as a substitute for a statistical analysis of another.

### 5. Swept founders — **HOLDS WITH QUALIFICATION**

Near fixation, trajectories can identify fitness differences much better than absolute fitnesses and founding times. That concern is sound.

However, a first mutation’s VAF near 0.5 means its **mutation-bearing lineage** is near fixation; descendants may account for most of that lineage. It does not establish that the original one-driver genotype itself swept. Consequently, the standalone-founder sweep-floor calculation need not apply.

Nor must an optimizer choose the birth boundary: it may choose another nearly equivalent point, and daughter trajectories can impose additional constraints.

**Correction:** describe weak identification and sensitivity to constraints, not a universal rule assigning the minimum possible founder fitness.

### 5b. Simulated acquisition-age errors — **HOLDS WITH QUALIFICATION**

The reported errors and missed early founding follow for the tested simulation/refit designs. The acknowledgement that median error \(+0.2\) is not unbiasedness is correct.

Three qualifications are essential:

* These are results for the **alternative refit**, not a demonstrated replay of the authors’ complete estimator.
* The full simulated case mix gives about **57%** by the Fig. 4b rule. **66.1%** appears only after imposing a real-data stratum mixture that the original model does not reproduce.
* That mixture uses **39 cases**, whereas the paper’s 68% is **32/47**. Even among all fitted cases, \(32/42=76.2\%\). The appropriate observed fraction among the same 39 eligible cases is not supplied.

**Correction:** retain the standardized simulation exercise, but remove the claim that it establishes “the paper’s 68% is what its model plus its estimator would produce.” The 18-year standard deviation belongs to the specified nonswept simulation stratum, not to real single-case ages generally.

### 6. Fitness ladder — **HOLDS WITH QUALIFICATION**

The reported median ratios support reasonably accurate **central recovery in the nonswept simulated stratum**. They do not show accurate recovery for every clone; the first- and fourth-driver interquartile ranges remain substantial.

The universal claim that swept founders receive about \(0.13/\mathrm{yr}\) “whatever the truth” is unsupported and contradicted by the other reported experiment: NOTE B §3b gives a published-model founder median of 0.23 fitted versus 0.18 true.

A birth-boundary estimate also does not, by itself, establish that its fitness should be discarded.

**Correction:** retain the demonstrated distortion in the particular \(r=1\), swept-founder experiment. Replace the universal floor claim with a sensitivity/identifiability warning. The criticism that geometric escalation is imposed in the simulator is sound, subject to the qualifications under NOTE A §5.

### 7. Late evolution — **HOLDS WITH QUALIFICATION**

The discrepancy between 0/928 and 10/47 is strong **conditional on comparable classification and sampling populations**. The numerical binomial bounds are consistent with the reported counts.

There is a useful confidence correction: the approximately **27-fold ratio already provides a conservative 95% lower confidence bound**. The lower endpoint of a two-sided 95% interval fails on its relevant side with probability at most 2.5%; the model’s upper endpoint likewise. Bonferroni therefore gives joint coverage at least 95% for their ratio. The 20-fold construction is valid but unnecessarily conservative.

The \(4\times10^5\) figure is not supported as an estimated biological discrepancy; see H4.

### 8. PPV decline — **HOLDS**

The reported horizon-100 and mortality-weighted results support withdrawing the claim that the decline is merely an age-85 cutoff artefact.

This establishes robustness within the reported model/scenarios, not calibration of real-world PPVs or validation of the chosen mortality model.

## 2. NOTE A: propositions, approximations and sections

### Proposition 1 — **HOLDS**

The conditional integrating-factor identity is correct, as is cancellation of the first-order term in the expected total.

There are two presentation corrections. In the proof, “the post-selection total is” must read **“the conditional expected post-selection total is.”** Also,
\[
E[N_{t+h}\mid\mathbf n]
=N_t\sum_i f_i e^{h(F_i-\bar F)}
\ge N_t
\]
by Jensen’s inequality. Selection therefore causes positive expected total-population inflation whenever fitness varies, in addition to mutation adding cells.

The deterministic skeleton is a conditional-mean/mean-field construction, not the exact unconditional mean trajectory of the stochastic system. The simulation does not maintain a fixed census population.

### Proposition 2 — **HOLDS under its stated constant-\(g\) assumption**

Poisson offspring with mean \(e^{gh}\), the extinction equation, the first-order survival expansion, and the martingale conclusions are correct. Finite variance supplies the required integrability.

In the actual competitive model, \(g=F-\bar F(t)\) changes. Thus “survival” here is survival of an auxiliary constant-environment branching process; its use as establishment probability in the full model is an approximation.

### Approximation A1 — **HOLDS WITH QUALIFICATION**

The linear-fractional result and the near-critical Poisson approximation are appropriate. The note correctly states that matching two logarithmic moments does not validate tails.

At larger \(gh\), the reported deviations are material. Neither the exponential approximation nor the constant-environment assumption should be treated as validated for rare-event calculations merely because these moment checks look reasonable.

### Corollary 3 — **HOLDS WITH QUALIFICATION**

The three early-stage scaling calculations are correct to leading order:
\[
NUq(g)\sim2NhUg,\qquad
n/N\sim Xe^{gt}/(2gNh),\qquad
Unq(g')\sim UXe^{gt}g'/g.
\]
They do not prove exact invariance of the complete finite-step model.

The convention discussion needs revision. Source 3 fixes ref 39’s normalization through its explicit rates. The excerpts from refs 38 and 41 do **not** specify that pairwise coalescence equals \(1/(N\tau)\). For example, a continuous-time binary self-renewal model with self-renewal rate \(1/\tau\) has pairwise coalescence approximately \(2/(N\tau)\). Matching that to \(1/(Nh)\) again gives \(\tau=2h\).

**Required correction:** remove the assertion that “coalescent units” necessarily imply the separate tenfold discrepancy. Ref 41’s wording “time between symmetric self-renewal cell divisions” particularly does not justify that assertion.

Also:

* Finite-\(h\) survival probabilities retain separate \(h\)-dependence.
* The fixed-\(Nh\) numerical difference is not shown to arise specifically from late-clone discretization rather than the collection of finite-step effects.
* Source 3 writes \(\theta=2N\tau\mu\) for a **haploid** mutation rate. Comparing this with NOTE B’s \(\theta=N\tau\mu\) requires an explicit haploid-versus-total-cell mutation-rate convention.

### Proposition 4 — **FAILS as a proposition about the reported establishment time**

Its exponential-amplitude matching algebra is correct. Its identification of the reported time is not.

Source 2 uses:
```python
if freq >= 1/(1.0*10**5):
    establishment_age = age
```
Writing \(S_i\) for unnormalized sizes and \(N_{\rm inf}=10^5\), the actual crossing condition is
\[
\frac{S_i(t)}{N_{\rm inf}+\sum_jS_j(t)}
\ge\frac1{N_{\rm inf}},
\]
equivalently
\[
S_i(t)\ge
\frac{N_{\rm inf}+\sum_{j\ne i}S_j(t)}
     {N_{\rm inf}-1}.
\]
It is only approximately \(S_i=1\) when other clones are negligible.

**Corrected proposition:** define a surrogate own-size time
\[
\widetilde t_e=t_0+\frac{\log(1+s)}s.
\]
For a clone fitted in the stipulated exponential phase against a constant background, the displayed Gumbel error law holds approximately for \(\widetilde t_e\), with
\[
\kappa=\log\frac{N_{\rm inf}s}{q(s)N}
\simeq\log\frac{N_{\rm inf}}{N\tau}
\]
in ref 39’s convention.

The stated numerical offsets are consistent with that surrogate calculation. They are not general corrections to v16’s reported times. For daughters, one must also distinguish absolute clone fitness from growth relative to the expanding background.

Finally, “irreducible” means irreducible from the stipulated later exponential-phase observations. Earlier observations or additional information about founding can reduce that uncertainty.

### Proposition 5 — **HOLDS WITH QUALIFICATION**

For absolutely continuous \(D\), the Poisson marking proof, hazard, conditional density, mean and CDF are correct.

The claimed extension to atomic delays needs a correction. The conditional delay is distributed according to \(dD(v)/D(T)\), with \(u=T-v\). Consequently,
\[
P(u\le x\mid T)
=1-\frac{D((T-x)^-)}{D(T)}
=1-\frac{h((T-x)^-)}{h(T)}.
\]
The left limit matters. For deterministic delay \(d\), all conditional mass is at \(u=T-d\); the printed CDF gives the wrong value at that point.

The interpretation also needs narrowing. The assumptions imply a nondecreasing hazard \(h(T)=\Lambda D(T)\); they are substantive restrictions, not consequences of age-homogeneous simulation rules alone.

Agreement with independently measured founding ages can support this **model class**. It cannot distinguish models within the class already calibrated to the same incidence curve. Thus “not independent support” should be replaced by “not discriminating support for this specific mechanism.”

### Approximation A2 — **HOLDS WITH QUALIFICATION**

The exponential-mixture recursion is correct for the specified independent branching construction **after an explicit terminal-size proxy is added**.

As written,
\[
G_4(a\mid F)=e^{-q(F)Ne^{-Fa}}
\]
calculates \(P(n_4(a)>N)\). In a pure-branching population retaining wild type and all ancestors, the actual half-population criterion is
\[
n_4>N+\sum_{j\ne4}n_j.
\]
For a fixed *total* population \(N\), it would instead be \(n_4>N/2\). These are different terminal events.

Other corrections:

* Waiting-time independence of \(h\) is only leading-order, after substituting \(q(F)\simeq2Fh\).
* The exact log-logistic variable in the displayed formula is \(e^{F\tau_2}-1\), not \(e^{F\tau_2}\).
* The tenfold incidence discrepancy does not isolate saturation as its cause; multiple approximations and the terminal criterion differ.

### Approximation A3 — **HOLDS WITH QUALIFICATION**

The reported 0.371% versus 0.464% supports “approximately 20% low.” The inert-competitor experiment supports only the tested, small additional effect; the note appropriately excludes other competition from that conclusion.

“Exact mean dynamics of Proposition 1” should be “deterministic conditional-mean update” or “mean-field dynamics.”

The sweep-time formula
\[
[\log(Nq(s))-\log X]/s
\]
is the approximate time to a **half-frequency crossing** in an appropriate two-type dominance approximation, not a generic fixation time.

The quoted 21-, 12- and 4.2-year stage durations need an explicit derivation. In particular, the expectation in \(UN\,E[q(s)]\) is over mutations **at birth**. Substituting conditioned AML-trunk medians does not establish that arrival rate or the conditional waiting-time distribution.

### Heuristic H4 — **FAILS as a quantitative estimate of the late-event probability**

The arithmetic of the displayed product is correct. The identification of that product with \(\Lambda_{\rm late}\) is not justified.

Multiplying stagewise \(M/F\) occupancy estimates does not enforce the joint event that the **entire lineage** stays below the detection threshold and then reaches diagnosis within two years. The displayed product contains no explicit two-year completion constraint. It also omits alternative paths involving transient, ultimately extinct genotypes.

Moreover, later mutations can occur **after** the lineage crosses the detection threshold and still complete rapidly. That possibility exists even with deterministic growth between random mutation events; unbounded Poisson size fluctuations are not the only objection to the claimed necessity.

**Correction:** retain this only as an illustrative calculation for a restricted path, with unknown approximation error. Delete the \(4\times10^5\) quantitative gap and the claim that a second channel of that magnitude is required. Zero events in 200,200 lives cannot validate a proposed probability near \(10^{-9}\).

### Section 5: conditioned fitness ladder — **HOLDS WITH QUALIFICATION**

The simulator imposes geometrically increasing increment scales. Even before conditioning, expected cumulative fitness is proportional to
\[
1+r+\cdots+r^{k-1}
=\frac{r^k-1}{r-1}.
\]
Its successive ratios approach \(r\); they are not exactly \(r\), especially at the first step.

The reported \(r=1\) experiment shows that AML conditioning does **not necessarily** create a geometric ladder. It does not prove that conditioning can never create an apparently geometric pattern under other distributions or sampling schemes.

**Correction:** describe the simulation’s geometric escalation as imposed rather than independently generated. Quantitative agreement can still test aspects of the model; it is not inherently illegitimate to calibrate a model containing that assumption. “At the input scales” should also be “near the input scales”—the first reported median is \(0.176\), not \(0.16\).

### Section 6: late evolution — **HOLDS WITH QUALIFICATION**

The empirical simulation-versus-data discrepancy survives; the H4-based mechanism and magnitude claims do not.

Comparability must be demonstrated through an observation model. Source 1’s classifier uses the complete first-driver lineage trajectory and its first-ever threshold crossing. Real classification uses limited samples and a finite panel. An unobserved initiating driver can create an apparent late pattern without requiring a separate biological initiation mechanism.

Also, \(M=200\) is approximate: the classifier uses 0.2% of a changing reconstructed population, not a fixed 200-cell threshold.

### Section 7: PPVs — **HOLDS WITH QUALIFICATION**

Under independent competing mortality, the weighting \(S(T)/S(a)\) is correct. The omitted-tail bound is also valid:
\[
P(\text{weighted AML after }100\mid E_a)
\le
\frac{S(100)}{S(a)}
P(T>100\mid E_a).
\]

But \(k\) and \(f\) alone do not specify a completion probability: fitness and other clones matter. Even in the reduction, a decreasing remaining horizon does not mathematically force the **pooled** PPV to fall; the eligible population’s composition also changes. The table establishes the decline, while the asserted enrichment for slow clones needs separate evidence.

At age 70, the table shows 62% diagnosed by 100, with about **40.4 percentage points** diagnosed between 85 and 100. Most progressors by 100 progress after 85; it does not establish that most eligible people progress after 85.

## 3. Additional errors and unsupported interpretations

### Unnoticed implementation issues in Source 2

`establishment_age` is not reset inside the clone loop. If a later clone never crosses the threshold, it can inherit the previous clone’s value; the `try/except` will not detect this. Initialize a per-clone fallback before searching. The excerpt establishes the potential bug, not its effect on published results.

There is also a time-grid inconsistency: trajectories use `i*dt`, but reported ages use `np.linspace(0, T*dt, T)`. These grids differ.

Source 1 likewise performs selection at index \(t\) and records diagnosis as `t*dt`; diagnosis labels run from 0 to 84.9 despite 850 full updates. These are small clock issues, but the notes should acknowledge them.

### Some refit comparisons use the wrong inferential reference

NOTE B §2c compares a complete published parameter vector with the optimum using a **one-parameter** likelihood threshold. That is not a calibrated joint test of a multiparameter vector. Use the appropriate joint test or profile a specified coordinate.

Different penalties, censoring, depths and objectives also prevent attributing objective differences solely to faulty optimization by the authors. Manual starting values do not themselves change the objective.

Pearson dispersion 3.14 measures residual disagreement relative to the fitted binomial model. It does not uniquely establish threefold biological year-to-year variance; measurement differences and mean-model misspecification can contribute.

Finally, simulation coverage and bias results used the pre-repair exponent cap. Its effects can propagate through the common denominator to other clones. Those results are not automatically an audit of the repaired refit.

### Boundary estimates are not automatically uninformative

NOTE B §2d partly dismisses two intervals lying entirely below 20 because their MLEs hit the birth constraint. A **lower** age constraint does not by itself explain an **upper** profile endpoint below 20. Daughter data and other restrictions may supply information.

The appropriate criticism is sensitivity to assumptions and constraints—not that every boundary estimate lacks evidential value.

### Posterior calculations require additional assumptions

NOTE A §2’s truncated log-exponential posterior is mathematically correct under its stipulated constant founding-time prior and simplified observation equation.

The daughter Gamma result requires an exponential **intensity prior**. For the first daughter event, the prior also contains
\[
\exp\!\left[-\int^{t_m}\lambda(v)\,dv\right],
\]
and an AML-ancestral daughter has additional selection. Parent-birth and observed-presence times impose further support restrictions.

Consequently, NOTE B §2a’s expected counts 22.3, 30.7 and 32.4 are not calibrated estimates of the real cohort’s true early-founding count without a justified prior, sampling model and treatment of fitness uncertainty.

### Finite samples and smaller paper criticisms

“No simulated diagnosis before 40” is not a structural zero of the model. Likewise, no nonswept founder among 400 \(r=1\) cases does not make that stratum impossible to study; it was absent from that run.

The 95%-versus-approximately-26% swept-founder discrepancy is worth investigating, but should be tested using the real sampling schedules, detection process and comparable founder definitions.

The count, docstring, sex-standardization and diagnostic-code discrepancies in NOTE B §4 are reasonable audit questions. The relevant pipelines and clinical definitions are not reproduced here, so they are not independently established by Sources 1–3. A male/female-average population benchmark is not automatically erroneous merely because the sequenced cohort consists of women; the intended target population needs clarification.

## Verdict to put to the authors

The strongest supported questions concern the **fivefold drift/population-scale mismatch**, the absence of a confidence-interval justification for optimizer ranges, weak identification of some swept-founder parameters, and the model’s failure to reproduce late patterns under its present classifier. The excessive frequency of early sweeps also warrants an observation-matched comparison.

Do **not** present the tenfold convention, a uniquely corrected 16% incidence, calibrated 52-year intervals, a real-cohort truth of 87% early founding, or a \(400{,}000\)-fold missing biological channel as established conclusions. The blanket fitness-ladder-artifact criticism is unsupported, and withdrawing the claim that the PPV decline is merely a cutoff artefact was justified.