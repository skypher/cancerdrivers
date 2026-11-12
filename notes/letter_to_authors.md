# Draft email to the corresponding authors (not sent)

To: Caroline J. Watson <cw672@cam.ac.uk>, Jamie R. Blundell <jrb75@cam.ac.uk> (corresponding authors,
per the article's correspondence section)
Subject: Questions on the simulation parameters and uncertainty estimates in "Evolutionary dynamics in the years preceding acute myeloid leukaemia" (Nature, 2026)

Dear Dr Blundell and Dr Watson,

Thank you for publishing the code and the public data alongside your paper; it made a close reading
possible. While working through the modelling we came across a few points we would be grateful if
you could clarify. We may well have missed something, and we would welcome a correction.

1. Population scale of the Fig. 5 simulation. The paper states Nτ ≈ 10^5 (refs 38, 39, 41) and
   U ≈ 10^-5 per year. sim_fast.py resamples every clone as Poisson(n·exp((F − x̄)dt)) with N = 10^5 and
   dt = 0.1 yr. Matching its drift variance and establishment probability to the birth-death rates in
   your Science 2020 paper (Supplementary Note 2: self-renewal and differentiation each at rate 1/τ)
   gives τ = 2·dt = 0.2 yr, i.e. Nτ = 2×10^4. In our C++ port of sim_fast.py (which reproduces your
   published counts, patterns and PPVs within Monte Carlo error), running at Nτ = 10^5 with
   U = 10^-5 gives a lifetime AML risk of about 6.6% rather than 0.47%. Re-tuning U to 3.8×10^-6 at
   Nτ = 10^5 restores the incidence, patterns, driver timing and fitness ladder, but changes the sweep
   PPVs (for example 61% instead of 46% at age 60). Is our reading of the intended convention right,
   and if so, which parameter values should readers take as the model's?
   A related observation: because the update draws a new Poisson generation every dt, changing dt at
   fixed N changes the drift scale (halving dt gives 0.15% lifetime risk), so dt acts as a model
   parameter rather than a numerical step.

2. Uncertainty of the fitness and establishment-time estimates. Supplementary Table 9 reports the
   spread of optimiser end points across 25 seeds within ±20% of the best log-likelihood. As far as we
   can see this measures optimiser variability rather than sampling uncertainty. In a re-fit of the
   public trajectories with the same clone-growth model (binomial likelihood, profile likelihood,
   allowing for the observed overdispersion), intervals for establishment times were much wider
   (median about 50 years). For several clones, quite different parameter values gave the same
   maximum likelihood to 0.001 (for example C92_011). This applies particularly to founders whose
   lineage was near fixation before the first sample, where the trajectories constrain fitness
   differences much better than absolute fitness and founding time. Would you consider reporting
   likelihood-based or bootstrap intervals for these estimates?

3. Two small points in Code_for_inferring_acquisition_age_and_fitness_v16.py (likelihood_clones):
   establishment_age is not reset between clones, so a clone whose frequency never reaches 1e-5
   would inherit the previous clone's value; and reported ages use np.linspace(0, T*dt, T) while
   trajectories use i*dt (differing by up to one step). We have not checked whether either affects
   the published values.

4. Late-evolution and no-driver cases. In an illustrative sensitivity analysis that samples
   simulated AML cases on the cohort's real schedules and treats drivers as missed independently by
   the panel, the estimated rate of cases with no detected trunk driver is about 1.6% even if 30% of
   drivers are off-panel; ignoring off-trunk detections makes this an estimated upper envelope for
   the model's no-panel-driver rate. The cohort has 5 of 47 such cases. This seems consistent with your suggestion that some AMLs
   follow routes outside the panel or the four-hit model, though other explanations remain possible,
   and we wondered whether you had considered adding such a route.

Our notes and code are at https://github.com/skypher/cancerdrivers (notes/recheck.md,
theory/theory.md).

With best regards,
Leslie Polzer
