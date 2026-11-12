# Draft email to the corresponding authors (not sent)

To: Caroline J. Watson <cw672@cam.ac.uk>, Jamie R. Blundell <jrb75@cam.ac.uk> (corresponding authors,
per the article's correspondence section)
Subject: Questions on the simulation parameters and uncertainty estimates in "Evolutionary dynamics in the years preceding acute myeloid leukaemia" (Nature, 2026)

Dear Dr Blundell and Dr Watson,

Thank you for publishing the code and the public data alongside your paper; it made a close reading
possible. While working through the modelling we came across a few points we would be grateful if
you could clarify. We may well have missed something, and we would welcome a correction.

Summary

- Our checks do not by themselves overturn the qualitative picture, but they leave uncertainty about
  individual clone ages and fitnesses and about some quantitative predictions of the Fig. 5 model.
- Simulation scale: matching, to leading order, the drift and establishment probability of
  sim_fast.py to the birth-death convention of your 2020 Science paper suggests an effective
  N*tau of about 2x10^4, five times below the stated 10^5. We would like to check this reading of the
  intended normalisation with you. In our implementation, the stated values (N*tau = 10^5,
  U = 10^-5) give a lifetime AML risk of about 6.6%.
- Uncertainty: the ranges in Supplementary Table 9 reflect optimiser variability across seeds rather
  than sampling uncertainty, and some parameters, notably those of founders whose lineage had swept
  before the first sample, appear weakly identified.
- Code: two possible issues in the inference code, whose effect on the published values we have not
  assessed.
- No-driver cases: in an illustrative sensitivity analysis the model rarely produces cases with no
  detected driver, whereas the cohort has 5 of 47.

Details

1. Population scale of the Fig. 5 simulation. The paper states N*tau ≈ 10^5 (refs 38, 39, 41) and
   U ≈ 10^-5 per year. sim_fast.py resamples every clone as Poisson(n·exp((F − x̄)dt)) with N = 10^5 and
   dt = 0.1 yr. Matching, to leading order, its drift variance and small-selection establishment
   probability to the birth-death rates in your Science 2020 paper (Supplementary Note 2:
   self-renewal and differentiation each at rate 1/τ) gives τ ≈ 2·dt = 0.2 yr, i.e. N*tau ≈ 2×10^4.
   Is that the intended normalisation? Our C++ implementation of sim_fast.py approximately reproduces
   the reported summaries we compared (case counts, pattern fractions, PPVs). With it, running at
   N*tau = 10^5 and U = 10^-5 gives a lifetime AML risk of about 6.6% rather than 0.47%. In one
   re-tuning (U = 3.8×10^-6 at N*tau = 10^5) the model approximately restored the listed simulated
   summaries (incidence, patterns, driver timing, fitness ladder), while the sweep PPVs changed in
   that run (for example 61% instead of 46% at age 60). We tried only one-parameter re-tunings, so
   these are illustrations rather than a recalibration. A related observation: because the update
   draws a new Poisson generation every dt, changing dt at fixed N changes the drift scale (halving
   dt gave 0.15% lifetime risk), so dt acts as a model parameter rather than a numerical step.

2. Uncertainty of the fitness and establishment-time estimates. Supplementary Table 9 reports the
   spread of optimiser end points across 25 seeds within ±20% of the best log-likelihood. As far as we
   can see this measures optimiser variability rather than sampling uncertainty, and the 20% window
   does not correspond to a coverage level. As a sensitivity diagnostic we re-fitted the public
   trajectories with an alternative fit that keeps the deterministic growth structure but changes
   the penalties, the observation likelihood and the establishment-time definition (the time at which
   a clone's own deterministic size reaches one cell). Its approximate profile intervals for that age
   had median width about 52 years. These are not calibrated confidence intervals for your reported
   establishment-time estimand. In that fit, different optimiser end points sometimes had objective
   values agreeing to the displayed 0.001 precision (for example C92_011), and parameters of founders
   whose lineage was near fixation before the first sample were only weakly constrained, the
   trajectories informing fitness differences more than absolute fitness and founding time.

3. Two possible code issues in Code_for_inferring_acquisition_age_and_fitness_v16.py
   (likelihood_clones): establishment_age is not reset between clones, so a clone whose frequency
   never reaches 1e-5 would inherit the previous clone's value; and reported ages use
   np.linspace(0, T*dt, T) while trajectories use i*dt (differing by up to one step). We have not
   checked whether either affects the published values.

4. Late-evolution and no-driver cases. In an illustrative sensitivity analysis that samples
   simulated AML cases on the cohort's real schedules, assumes four trunk drivers each missed by the
   panel independently, and ignores off-trunk detections, the estimated rate of cases with no
   detected trunk driver is about 1.6% (case-bootstrap 95% interval 1.4–1.7%) even if 30% of drivers
   are off-panel; this is an estimated upper envelope for the model's no-panel-driver rate under
   those assumptions. The cohort has 5 of 47 such cases. This seems consistent with your suggestion
   that some AMLs follow routes outside the panel or the four-hit model, though it does not identify
   that explanation (correlated misses, other driver counts and detection differences remain
   possible).

Possible checks and their purpose

1. Clarify the simulation parameters (N*tau, U and the role of dt): state the effective scale if it is
   the intended one, or re-run Fig. 5 at the stated parameters with a recalibration. In our U
   re-tuning the PPVs changed (61% instead of 46% at age 60) while the listed summaries were
   approximately preserved; how they change under your intended recalibration remains to be seen.
   Consistent parameters would also let readers reuse the model.
2. Report likelihood-profile or bootstrap intervals for clone fitness and establishment time, and
   flag weakly identified founders. This would give an interpretable uncertainty assessment for your
   estimand (its width and shape would need calibration), and the share of first drivers acquired by
   age 20 would rest on a population-level analysis rather than per-case point estimates.
3. Check the two code paths in item 3, which would determine whether the reported establishment
   times are affected.
4. Compare the evolutionary patterns between model and cohort under matched observation (real
   sampling schedules, panel and detection threshold, with the classification rule made explicit).
   This would make that comparison more directly interpretable. A route for no-driver and late cases
   is one modelling option, not something this analysis establishes.

Our notes and code are at https://github.com/skypher/cancerdrivers (notes/recheck.md,
theory/theory.md).

With best regards,
Leslie Polzer
