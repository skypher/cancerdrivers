# Plan: statistics re-check and theory for Watson et al., Nature 2026

Paper: Watson, Poon, MacGregor, ..., Menon, Blundell, "Evolutionary dynamics in the
years preceding acute myeloid leukaemia", Nature, published online 2026-11-11,
doi:10.1038/s41586-026-11269-x. Local copies in `sources/nature/`, code in
`sources/code/` (cloned 2026-11-12, upstream `the-blundell-lab/preAML_evolutionary_dynamics`).

User request (2026-11-12): do both (1) the statistics re-check and (2) the analytic
theory of the paper's simulation model.

## Part 1: statistics re-check

1a. Simulation model as implemented (`sim_fast.py`) versus as stated (Methods: Nτ = 1e5,
    continuous-time birth-death). The code is a Poisson Galton-Watson chain with
    generation dt = 0.1 yr, so the drift scale is N·dt = 1e4, not 1e5. Check the
    consequences (incidence, timing) and the dt dependence.
1b. Acquisition ages. The fitted t0 (and the reported "establishment time", where the
    deterministic size reaches one cell) confounds the true mutation time with the
    random establishment factor X. Quantify the bias and spread under the authors' own
    model: simulate, fit their deterministic formula, compare with the true times.
    Recompute the Fig. 4b claims ("first driver by age 20 in two thirds of cases").
1c. Uncertainty. Their intervals are the spread of optimizer end points over 25 seeds
    within ±20% of the best negative log-likelihood. Re-fit the public trajectories with
    the same deterministic model, binomial and beta-binomial likelihoods, and give
    profile-likelihood intervals; measure overdispersion.
1d. The geometric increase of fitness with driver number: test how much of it the
    monotonicity penalty and conditioning on AML produce by themselves (null model r = 1).
1e. Smaller items: code constants versus paper (26/12/9 vs 26/11/10), UK Biobank C92
    outcome includes non-AML myeloid leukaemias.

## Part 2: theory of the simulation model

2a. Exact statement of the model as implemented.
2b. Branching-process recursion for P(AML by age T) with exponential establishment
    factors; compare with simulation.
2c. Conditional law of the AML lineage (fitness and timing by stage); which features
    come from conditioning alone.
2d. Late evolution: bound on the probability that the first driver becomes detectable
    within 2 years of diagnosis; what model change would give about 20%.
2e. Positive predictive value of sweeps from the recursion.

## Rules in force

- Computations in C++ (clang), progress lines, runs under `runs/`, at most 8 threads
  (machine shared, load about 24 of 27 cores at 08:00 Taiwan time).
- Claims labelled: derived (with proof), numerical, or conjecture.
- Reviews: Jander and yang2 on the shared ntfy topic before calling anything final.
