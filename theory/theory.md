# Theory of the clone model in Watson et al., Nature 2026

Subject: the stochastic model used for Fig. 5 and Extended Data Fig. 10 of Watson, Poon,
MacGregor, ..., Menon, Blundell, *Evolutionary dynamics in the years preceding acute myeloid
leukaemia*, Nature (online 2026-11-11), doi:10.1038/s41586-026-11269-x, as implemented in
`sim_fast.py` of github.com/the-blundell-lab/preAML_evolutionary_dynamics (commit 5a53fe7).

Status labels. **Proposition** = proved here under the stated assumptions. **Approximation** =
derived by an explicitly named approximation, then compared with simulation. **Numerical** =
simulation only. Simulations use the C++ port `sim/amlsim.cpp`, which reproduces the published
numbers (below).

## 0. The model as implemented

Time step $h=\mathrm{dt}=0.1$ yr, $T=850$ steps (85 years). A person starts with one wild-type
clone of $N=10^5$ cells and fitness $0$. A clone is a set of cells with the same drivers; clone
$i$ has size $n_i\in\mathbb N$, fitness $F_i$ (sum of its drivers' increments) and $k_i\le 4$
drivers. Each step:

1. *Mutation.* Every clone with $k_i<4$, $n_i>0$ produces $\mathrm{Poisson}(n_i U h)$ new clones
   ($U=10^{-5}$ per cell per year). A new clone starts with one cell, $k_i+1$ drivers and fitness
   $F_i+s$, where $s$ is drawn from the table of `mutation_fitness_generator` for scale
   $a_{k_i}=s_1 r^{k_i}$ ($s_1=0.16$, $r=2.5$): a quantised stretched exponential
   $\rho_a(s)\propto e^{-(s/a)^3}$ with 18 atoms at $0.1a,0.2a,\dots,1.8a$.
2. *Mean fitness* $\bar x=\sum_i F_i n_i/\sum_i n_i$.
3. *Selection and drift.* $n_i\leftarrow \mathrm{Poisson}\big(n_i e^{(F_i-\bar x)h}\big)$.
4. Clones that are extinct, childless and never reached 10 cells are deleted.
5. *AML* when some 4-driver clone exceeds half of all cells.

Reproduction (200,200 lives, `runs/repro/`): 928 AML cases (0.464%; published 943, 0.471%),
patterns linear/branched/late 637/291/0 (published 642/301/0), PPV of "largest VAF > 30% and
second VAF > 10%" at ages 40/50/60/70: 50/56, 217/303, 553/1204, 733/3422 (published 44/50,
214/292, 564/1152, 719/3313). All differences are within Monte Carlo error.

## 1. Mean dynamics and establishment

**Proposition 1 (selection step in expectation).** Let $n_i$ be the sizes after the mutation step
of a time step (so they include the clones born in that step), let $\bar x$ be the mean fitness
computed from them, frozen during the selection step, and let $\tilde n_i(t)=n_i(t)\exp\big(\int_0^t\bar x\big)$
with $\bar x$ piecewise constant. Conditionally on this post-mutation state,
$E[\tilde n_i(t+h)]=\tilde n_i(t)\,e^{F_i h}$ for every clone, and the selection step conserves the
post-mutation total to first order: $E[\sum_i n_i(t+h)]=\sum_i n_i\,(1+O(h^2\,\mathrm{Var}_f F))$.
(The total over a whole step also gains the newborn mutants of the mutation step.)

*Proof.* $E[n_i(t+h)\mid\cdot]=n_i e^{(F_i-\bar x)h}$; multiply by $e^{\int_0^{t+h}\bar x}$. The
conditional expected post-selection total is $\sum_i n_i e^{(F_i-\bar x)h}
=\sum_i n_i\big(1+(F_i-\bar x)h+\tfrac12(F_i-\bar x)^2h^2+\dots\big)$, and the first-order term
vanishes by the definition of $\bar x$. $\square$

By Jensen's inequality the conditional expected total is at least the post-mutation total, so
selection inflates the population in expectation whenever fitness varies (besides the newborns);
the simulation does not hold a fixed census size. The conditional-mean (mean-field) skeleton
therefore has the structure of the paper's inference model (Methods,
eqs. 1-2): wild type of constant unnormalised size, clones growing exponentially, frequencies by
normalisation. The two differ near a clone's birth (the paper's $S(t)=(e^{s(t-t_0)}-1)/s$ against
a random establishment factor times $e^{st}$) and in drift.

**Proposition 2 (establishment).** A clone born with one cell whose fitness exceeds the mean by
a constant $g>0$ is a Galton-Watson process with $\mathrm{Poisson}(m)$ offspring, $m=e^{gh}$, and
generation time $h$. It survives with probability $q(g)$, the root in $(0,1)$ of $1-q=e^{-mq}$
(the positive root; $q=0$ is the other), and $q(g)=2gh\,(1+O(gh))$. Almost surely
$n(t)e^{-gt}\to W$ (unconditional limit) with $E[W]=1$, $\{W>0\}=\{\text{survival}\}$, hence
$E[W\mid\text{survival}]=1/q$. In the competitive model $g=F-\bar x(t)$ changes with time, so using
$q$ as the establishment probability there is an approximation (survival of an auxiliary
constant-environment process).

*Proof.* Standard Galton-Watson theory: the extinction probability is the smallest fixed point of
the offspring generating function $e^{m(z-1)}$, and $n_t/m^{t/h}$ is a non-negative martingale with
mean 1 that converges to a limit $W$ with $\{W>0\}=\{\text{survival}\}$ because the offspring law
has finite variance (Kesten-Stigum). Expanding $1-q=e^{-mq}$ gives $q=2(m-1)/m+O((m-1)^2)$. $\square$

**Approximation A1 (exponential establishment factor).** On survival, $qW\approx X$ with
$X\sim\mathrm{Exp}(1)$. This is exact for linear-fractional offspring laws and holds in the
near-critical limit $gh\to 0$. The check below compares two moments of $\log(qW)$ only; it does not
validate the tails. *Numerical check* (`sim/wlaw.cpp`, 20,000 lineages per value):
mean and s.d. of $\log(qW)$ against $\log X$ ($-0.577$, $1.283$): $g=0.16$: $-0.63$, $1.27$;
$g=0.5$: $-0.56$, $1.23$; $g=1$: $-0.54$, $1.23$; $g=2.5$: $-0.51$, $1.21$; $g=4$: $-0.44$, $1.10$.
The deviations at larger $gh$ are material; A1 is not validated for rare-event calculations.

**Corollary 3 (the population scale is $N h$, not $N$; under Proposition 2 and A1).** For clones
that are small compared with $N$ (so $\bar x\approx0$):
- first drivers are established at rate $\Lambda_1=NU\,E[q(s)]\approx 2\,(Nh)\,U\,E[s]$;
- an established clone of fitness $g$ has $n(t)\approx X e^{gt}/q(g)\approx X e^{gt}/(2gh)$
  cells, so its frequency is $X e^{gt}/(2g\,Nh)$;
- it founds established daughters of fitness $g'$ at rate $U\,n(t)\,q(g')\approx U X e^{gt}\,g'/g$,
  which does not depend on $h$.

These three rates, and the frequency threshold for AML, depend on $N$ and $h$ only through $Nh$
to leading order; finite-$h$ survival probabilities keep a separate dependence on $h$. That the
model's outputs inherit the scaling is not proved here; the table below tests it numerically.

*What $N\tau$ the code corresponds to.* The paper takes $N\tau\approx10^5$ from its refs 38, 39 and
41. Ref 39 (Watson et al., Science 2020; bioRxiv 569566, Supplementary Note 2) fixes the
normalisation by explicit per-cell rates: symmetric self-renewal at $b=\lambda r(1+\sigma/2)$ and
symmetric differentiation at $d=\lambda r(1-\sigma/2)$, with $1/\tau=\lambda r$ and $\sigma$ the
dimensionless selection coefficient per $\tau$ (written $s$ there). The per-year growth rate is
$g=b-d=\sigma/\tau$, the drift variance is $(b+d)n=2n/\tau$ per year, and a new clone survives with
probability $q=1-d/b=\sigma/(1+\sigma/2)=g\tau/(1+g\tau/2)\approx g\tau$. The Poisson update has
variance $n e^{gh}/h\approx n/h$ per year and survival $q\approx2gh$, so matching the leading-order
variance and survival gives $\tau=2h=0.2$ yr and $N\tau=2\times10^4$, five times below the stated
value; the stated value corresponds to $Nh=5\times10^4$. Because $\tau=2h$, changing $h$ at fixed
$N$ changes the physical drift scale $N\tau$: in this code $h$ is a model parameter, not a numerical
step that can be refined at fixed model. (That paper's text calls $1/\tau$ "the number of symmetric divisions per year",
while its rates give $2/\tau$; the value used here is the one fixed by its equations, $c=1$ and
$\tilde n=(e^{sT}-1)/s$.) Refs 38 and 41 estimate $N\tau$ phylodynamically from coalescence rates;
the quoted excerpts do not fix the factor-of-two normalisation, and a binary self-renewal model
with self-renewal rate $1/\tau$ has pairwise coalescence rate about $2/(N\tau)$, which matched to the
Poisson update's $1/(Nh)$ again gives $\tau=2h$. (An earlier version of this note claimed a separate
tenfold gap in "coalescent units"; that was not established and is withdrawn.) Ref 39 writes
$\theta=2N\tau\mu$ with a haploid mutation rate, whereas the paper's Discussion writes
$\theta=N\tau\mu$; comparing them needs an explicit mutation-rate convention.

*Numerical* (`runs/ndt/`): lifetime AML risk to age 85

| $N$ | $h$ | $Nh$ | risk |
|---|---|---|---|
| $10^5$ | 0.1 | $10^4$ | 0.464% (published setting) |
| $5\times10^5$ | 0.02 | $10^4$ | 0.536% |
| $10^5$ | 0.05 | $5\times10^3$ | 0.146% |
| $5\times10^5$ | 0.1 | $5\times10^4$ | 6.6% (stated $N\tau$ in ref-39 units, other parameters unchanged; 20,000 lives) |
| $10^6$ | 0.1 | $10^5$ | 16.0% (sensitivity run; 4,000 lives) |

At fixed $Nh$ the risk changes by 15% over a fivefold change of $h$ (finite-step effects; not
attributed to a specific mechanism). At the stated $N\tau=10^5$ in ref-39 units, with all other
parameters unchanged, the model gives a lifetime risk of 6.6%, 14 times the observed 0.48%.
Whether the model retuned at the stated scale reproduces incidence and the paper's other outputs
is examined in §8.

## 2. Acquisition ages

**Proposition 4 (the own-size establishment time; under A1).** Fit the paper's deterministic
size $S(t)=(e^{s(t-t_0)}-1)/s$ against a population of $N_{\rm inf}=10^5$ to a clone that in
truth was founded at $t_m$ and grows at rate $s$ relative to the background, with
$n(t)\approx X e^{s(t-t_m)}/q(s)$ in a population of $N$ cells. Matching the exponential phases
gives
$$t_0=t_m-\frac{1}{s}\Big[\log X+\log\frac{s N_{\rm inf}}{q(s)N}\Big]
\;\approx\;t_m-\frac{1}{s}\Big[\log X+\log\frac{N_{\rm inf}}{N\tau}\Big]$$
with $q(s)\approx s\tau$ in the convention of ref 39 (in the code $\tau=2h$, so $N\tau=2\times10^4$),
and define the own-size time (deterministic size reaches one cell)
$\tilde t_e=t_0+\log(1+s)/s$. (The paper's reported `establishment_time` is a different quantity:
the first time the clone's normalised frequency $S_i/(N_{\rm inf}+\sum_jS_j)$ reaches $1/N_{\rm inf}$,
i.e. $S_i\ge(N_{\rm inf}+\sum_{j\ne i}S_j)/(N_{\rm inf}-1)$; the two agree only when other clones are
negligible at that time, and differ for daughters growing inside large ancestors.) Since $-\log X$ is Gumbel distributed with mean $\gamma=0.577$ and
variance $\pi^2/6$:
$$E[\tilde t_e-t_m]=\frac{\gamma-\kappa+\log(1+s)}{s},\qquad \mathrm{sd}(\tilde t_e-t_m)=\frac{\pi}{\sqrt6\,s}\approx\frac{1.28}{s},
\qquad \kappa=\log\frac{N_{\rm inf}\,s}{q(s)\,N}\approx\log\frac{N_{\rm inf}}{N\tau}.$$

*Proof.* Equate $e^{-st_0}/(sN_{\rm inf})$ with $X e^{-st_m}/(q(s)N)$; use $q(s)\approx s\tau$
(Proposition 2: $q\approx2sh$, and $\tau=2h$). $\square$

Consequences.
- Proposition 4 averages over $X$ for a clone with a given founding time. In samples selected on
  the trajectory (for example AML cases whose founder is still small years before diagnosis) the
  selection acts on $X$; the simulation tests in `notes/recheck.md` §2b2 show this.
- For an unselected clone observed only in its later exponential phase, the spread $1.28/s$ does not
  shrink with more such data: it is the luck of the clone's first few divisions (13 years at
  $s=0.1$/yr). Earlier observations or other information on founding could reduce it. After
  selection on the trajectory both the offset and the spread change, so $1.28/s$ is not a universal
  error for the paper's reported ages.
- The offset of $\tilde t_e$ depends on the population scale. In the code's own model
  ($\kappa=\log 5$) it is on average $(\log 5-\gamma-\log(1+s))/s$ too early: 9.4 years at $s=0.1$,
  2.6 years at $s=0.3$. At the stated $N\tau=10^5$ in ref-39 units ($\kappa=0$) it is too late by
  $(\gamma+\log(1+s))/s$: 6.7 and 2.8 years. These are statements about the surrogate $\tilde t_e$,
  not general corrections to the paper's reported times.
- Posterior of the true founding age given $(s,t_0)$, with constant first-driver supply (uniform
  prior on $t_m\ge0$) and the simplified observation equation above: $t_m=t_0+(L+\kappa)/s$ with
  $L=\log X$ truncated to $t_m\ge 0$. For a daughter founded inside a parent of fitness $s_p$, taking
  an exponential intensity prior $\propto e^{s_p t_m}$ turns $e^{L}$ into a $\mathrm{Gamma}(1+s_p/s)$
  variable; for the first daughter event the prior also carries a survival factor
  $\exp[-\int^{t_m}\lambda]$, and selection on AML ancestry and parent-birth times restrict it further.
  The expected counts in `notes/recheck.md` §2a rest on these simplifications and are not calibrated
  estimates for the real cohort.

## 3. When was the first driver acquired?

**Proposition 5 (first-driver age law).** Assume (A) founders of AML lineages are established at
the points of a Poisson process of constant rate $\Lambda$ on $[0,\infty)$, and (B) the lineage
founded at $u$ causes diagnosis at $u+\theta$, where the $\theta$'s are i.i.d., independent of the
founding times, with defective distribution function $D$ ($D(0)=0$), absolutely continuous with
density $D'$. (For general $D$ the delay given diagnosis at $T$ has law $dD(v)/D(T)$ with $u=T-v$,
and $P(u\le x\mid T)=1-h((T-x)^-)/h(T)$ with a left limit.) Let diagnosis be the first completion. Then the diagnosis hazard is $h(T)=\Lambda D(T)$, and given diagnosis at age $T$ the
founding age $u$ of the diagnosing lineage has density
$$f(u\mid T)=\frac{h'(T-u)}{h(T)},\quad 0\le u\le T;\qquad
E[u\mid T]=\frac{H(T)}{h(T)},\qquad P(u\le x\mid T)=1-\frac{h(T-x)}{h(T)},$$
with $H(T)=\int_0^T h$.

*Proof.* By the marking theorem the pairs $(u_i,u_i+\theta_i)$ form a Poisson process on
$\{0\le u\le c\}$ with intensity $\Lambda\,du\,dD(c-u)$. The number of completions by $T$ is
Poisson with mean $m(T)=\Lambda\int_0^T D(T-u)\,du=\Lambda\int_0^T D$, so
$P(\text{no diagnosis by }T)=e^{-m(T)}$ and $h=m'=\Lambda D$. The probability that the first
completion falls in $[T,T+dT]$ and comes from a founder in $[u,u+du]$ is
$\Lambda\,du\,D'(T-u)\,dT\,e^{-m(T)}$; dividing by $h(T)e^{-m(T)}dT$ gives
$f(u\mid T)=\Lambda D'(T-u)/h(T)=h'(T-u)/h(T)$. Integrating, $P(u\le x\mid T)=(h(T)-h(T-x))/h(T)$
and $E[u\mid T]=\int_0^T u\,h'(T-u)\,du/h(T)=\big[-u\,h(T-u)\big]_0^T/h(T)+\int_0^T h(T-u)\,du/h(T)=H(T)/h(T)$,
using $h(0)=\Lambda D(0)=0$. $\square$

Here a "founder" is a first-driver mutation whose clone survives drift, and $u$ is the time of
that mutation. In the simulation model (A) holds approximately: established founders arrive at
rate $U\,n_w(t)\,E[q(s-\bar x(t))]$, which is constant while the wild type $n_w\approx N$ and
$\bar x\approx0$, and falls once the person's mutant clones expand. (B) holds up to competition
between lineages of the same person; nothing in the model depends on age. The assumptions imply a
non-decreasing hazard $h=\Lambda D$; they are substantive restrictions, not consequences of
age-homogeneous rules alone.

*Numerical check* (`analysis/firstdriver_law.py`, 928 simulated cases; $h$ estimated from the
simulated diagnosis ages by a kernel estimate reflected at 85): law vs. simulation,
$P(u\le 20\mid T)$ averaged over cases 0.932 vs. 0.933; mean first-driver age 8.06 vs. 8.07
years. By quartile of diagnosis age the law matches to within 0.03 in probability and 1.6 years
in mean.

Applied to UK incidence (CRUK rates used in the paper's Fig. 5c), for a woman diagnosed at 72:
$P(u\le 20)=0.74$, $E[u]=14$-$16$ years (female rates; 0.82 and 12-13 years with male rates).

Consequences.
- The paper's Fig. 4b statement (first driver by age 20 in two thirds of cases) is about what any
  model satisfying (A) and (B) and calibrated to the incidence curve predicts. Agreement supports
  the model class; it does not discriminate the paper's specific mechanism from others in that
  class. The value at diagnosis age 72 is not the average over the cohort's diagnosis ages.
- The paper's own simulation, whose incidence is near zero before 40 and therefore steeper,
  predicts 93%. In simulated cases from that model the trajectory estimator under-reports early
  founding (`notes/recheck.md` §2b, §2b2); this is a statement about the tested designs, not a
  calibration of the paper's 68%.

## 4. Stage structure of the AML lineage

*Numerical* (928 simulated cases, `runs/repro/trunkstats.txt`): the 2nd driver arises while
clone 1 is small (median frequency 0.07, mean fitness 0.013). The 3rd and 4th arise when the
parent clone already dominates (median parent frequency 0.62 and 0.78, mean fitness 0.39 and
1.31, close to the parent's fitness). Median founding ages of drivers 1-4: 6.0, 39.1, 61.1,
74.1; diagnosis 79.1.

**Approximation A2 (pure branching).** Ignore competition ($\bar x=0$, no saturation), take each
established clone's growth as deterministic given its factor, $n(t)=We^{Ft}$, and the daughter
subtrees as independent given $W$. With A1 the
probability $G_k(a\mid F)$ that an established $k$-driver clone of fitness $F$ yields AML within
time $a$ satisfies, exactly under these assumptions,
$$G_k(a\mid F)=\frac{X_k(a)}{q(F)+X_k(a)},\qquad X_k(a)=U\!\int_0^a e^{Fv}\,E_s\big[q(F+s)\,G_{k+1}(a-v\mid F+s)\big]dv,$$
with the terminal proxy $G_4(a\mid F)=\exp(-q(F)N e^{-Fa})=P(n_4(a)>N)$, and
$h(T)=NU\,E_s[q(s)G_1(T\mid s)]$, because $E[e^{-WB}]=q/(q+B)$ for $W\sim\mathrm{Exp}(q)$. The proxy is
not the model's AML event: with wild type and ancestors retained, the half-population criterion is
$n_4>N+\sum_{j\ne4}n_j$ (for a fixed total it would be $n_4>N/2$). The 2nd-driver waiting time in
this regime satisfies $P(\tau_2>t)=\big(1+U Q_2 (e^{F_1t}-1)/(q(F_1)F_1)\big)^{-1}$ with
$Q_2=E_s[q(F_1+s)]$, i.e. $e^{F_1\tau_2}-1$ is log-logistic; it is independent of $h$ only to
leading order ($q\approx2Fh$). Result (`sim/theory_branch.cpp`): lifetime risk 4.98%, ten times the
simulation. The stage statistics above point to saturation of parent clones, but the tenfold gap
mixes several approximations and the terminal proxy and does not isolate that cause.

**Approximation A3 (single lineage, deterministic growth).** Keep drift only at establishment
(Proposition 2, A1), let all established clones follow the conditional-mean (mean-field) dynamics
of Proposition 1, and let only the newest clone of one lineage produce daughters. Result
(`sim/theory_reduced.cpp`, $3\times10^6$ lineages): lifetime risk 0.371% (simulation 0.464%),
mean founding ages of drivers 1-4: 8.2, 37.4, 58.2, 70.9 (simulation 8.1, 38.1, 59.9, 72.4).
Adding the person's other first-driver clones as inert competitors: 0.365% (mean founding ages
7.9, 37.1, 58.0, 70.9), a change of under 2%, within the Monte Carlo error of about 1.2% per run.
So the full model's incidence is reproduced to within about 20% by one lineage with
deterministic logistic growth after a random establishment factor. Adding the person's other
first-driver clones as inert competitors changes it by under 2%; other competition (their
daughters, off-trunk subclones) is not tested. Conjecture, not tested: the remaining 20% comes
from the restriction that only the newest clone of the lineage founds daughters.

In a two-type dominance approximation (A1 and Proposition 1), a clone with increment $s_k$ over a
dominant parent reaches half frequency about $(\log(Nq(s_k))-\log X)/s_k$ years after founding. An
earlier version derived stage durations (21, 12, 4.2 years) from this by inserting conditioned
increments into the at-birth arrival rate $UN\,E[q(s)]$; that substitution is not justified, and
those figures are withdrawn. The simulated medians are 21.1, 11.9 and 4.4 years.

## 5. Fitness of the drivers on the path to AML

*Numerical.* Conditioned on AML, the increments of drivers 1-4 have medians 0.176, 0.40, 1.00,
2.50: near the scales $a_k$ of the input distributions, about twice their unconditional means (0.089,
0.223, 0.557, 1.394). Clone fitness rises 0.18, 0.57, 1.54, 3.93 (ratios 3.2, 2.7, 2.6).

*Reason.* Conditioning on completing four stages within a lifetime tilts each increment's law by
a completion probability that rises steeply with $s$. For the light-tailed input
$\rho_a\propto e^{-(s/a)^3}$ the tilted law concentrates near the top of the input's bulk, at a
fixed multiple of $a$. The geometric escalation is therefore imposed by the input scales: even
before conditioning, expected cumulative fitness is proportional to $(r^k-1)/(r-1)$, whose
successive ratios approach $r$ (not exactly $r$, especially at the first step). Calibrating a model
that contains this assumption is legitimate, and quantitative agreement can test other aspects of
it, but the model's geometric ladder is an input rather than an independent prediction.

*Null model* (`runs/null/`): $r=1$, with $s_1=0.5$ chosen to match incidence (0.515% to 85,
1031 cases in 200,200 lives). Conditioned increments: medians 0.55, 0.50, 0.50, 0.55, i.e. the
common input scale at every stage. Clone fitness: 0.55, 1.05, 1.60, 2.15 (ratios 1.9, 1.5, 1.3),
an additive ladder. In this model, conditioning on AML does not create a geometric ladder; that does
not show it could never do so under other distributions or sampling schemes. Median founding ages of drivers 1-4: 4.6, 25.0, 44.0, 62.0 (published model:
6.0, 39.1, 61.1, 74.1).

## 6. Late evolution

The pattern "late" requires the first driver's lineage to stay below 0.2% of cells (about
$M=200$ cells; the classifier uses 0.2% of a changing population) until two years before diagnosis, and its 4-driver clone to exceed 50% at
diagnosis. Heuristically, with deterministic exponential growth, growing by a factor $\ge 250$ in
two years takes a growth rate above $\log 250/2=2.8$/yr, i.e. at least a 3-driver clone, so drivers
2 and 3 (and in practice 4) would have to be established while the lineage has at most $M$ cells.
This is not a necessary condition for the stochastic model: Poisson updates have unbounded upward
fluctuations, and later drivers can also arise after the lineage crosses the threshold and still
complete quickly.

**Heuristic H4 (illustrative restricted-path calculation; not a quantitative estimate).** A clone of fitness $F$ founds established daughters of
fitness $F'$ during the time it has at most $M$ cells at expected number
$\approx U\,q(F')\,M/F\approx 2hUM\,F'/F$. Chaining stages 2-4 and integrating over first
drivers ($q(F_1)\approx 2hF_1$ cancels the $1/F_1$):
$$\Lambda_{\rm late}\approx NU\,(2h)^4\,U^3M^3\,E[F_4]\approx 3\times10^{-11}\ \text{per year},$$
(with $NU=1$, $(2h)^4=1.6\times10^{-3}$, $(UM)^3=8\times10^{-9}$, $E[F_4]\approx2.3$), about
$2.5\times10^{-9}$ per lifetime. This product does not enforce the joint event that the whole lineage
stays below threshold and then completes within two years, contains no two-year completion
constraint, ignores paths through transient genotypes and later drivers arising after the
threshold crossing, and replaces stochastic occupancy by $M/F$. Its approximation error is
unknown, and it is not an estimate of the model's late rate. (An earlier version derived from it a
$4\times10^5$-fold gap and a required second channel of that size; withdrawn. Zero events in
200,200 simulated lives cannot validate a probability near $10^{-9}$.)

*Statistical statement* (no approximation): the simulation gives 0 late cases among 928 AML cases
(published: 0 of 943). The exact 95% upper bound on the model's late fraction is 0.40% of AML cases
(two-sided Clopper-Pearson; 0 of 200,200 lives: about $1.84\times10^{-5}$ per life). The observed
fraction is 10/47, two-sided 95% interval 10.7%-35.7%. Each interval fails on its relevant side
with probability at most 2.5%, so the endpoint ratio, about 27 (26.98), is a joint 95% lower
confidence bound on the ratio of late fractions (the Bonferroni construction at 97.5% each, giving
20, is valid but more conservative). Assumptions: independent cases; simulated and real cases
classified comparably. Comparability needs an observation model: the simulation classifies on the
complete first-driver lineage trajectory and its first-ever threshold crossing, the real data on
sparse samples and a finite panel, where an unobserved initiating driver can produce an apparent
late pattern without a separate biological mechanism. A sensitivity analysis (`notes/recheck.md` §2e: real sampling schedules, 0.1% VAF detection,
off-panel drivers, two artificial clone-selection rules) reduces the apparent discrepancy: the late
fraction ranges from 0% with every driver visible to values inside the observed interval at 30%
off-panel under one rule, so no robust fold-gap follows; cases with no visible trunk driver stay at
or below 1.6% against 5 of 47 real no-driver cases.

## 7. Positive predictive value of sweeps

The PPV of Fig. 5e averages, over everyone meeting the criterion at age $a$, the probability of AML
before the horizon. That probability depends on the person's whole clonal state (fitnesses and
other clones), not only on the number of drivers and the frequency. Less remaining time pushes the
pooled PPV down, but the composition of the eligible pool also changes with age, so a decline is
not forced mathematically; the table below establishes it in the model.

*Numerical* (`runs/horizon/`, the same 200,200 simulated lives followed to 100 instead of 85; the
lives are identical up to 85 because each life has its own random stream). PPV of "largest VAF
> 30% and second VAF > 10%":

| age | horizon 85 | horizon 100 |
|---|---|---|
| 40 | 50/56 = 89% | 56/56 = 100% |
| 50 | 217/303 = 72% | 282/303 = 93% |
| 60 | 553/1204 = 46% | 962/1204 = 80% |
| 70 | 733/3422 = 21% | 2117/3422 = 62% |

Without death from other causes, 62% of those meeting the criterion at 70 are diagnosed by 100,
40.4 points of them between 85 and 100: most progressors by 100 progress after 85. The model has no competing mortality, so the 85-year cut-off acts as a crude stand-in
for it. Weighting each simulated AML by the probability of surviving other causes from the
screening age $a$ to diagnosis, $S(T)/S(a)$ (the person is alive at screening; other-cause death
assumed independent of the clone dynamics; Gompertz mortality $\mu(x)=Ae^{Bx}$, $A$ set so that life expectancy at 70 is 16.5
years, roughly UK women; `analysis/ppv_mort.py`) gives PPVs of 81%, 64%, 44%, 28% at 40, 50, 60,
70 for $B=0.1$ (78-83%, 61-67%, 43-46%, 28% for $B=0.085$-$0.115$). These are close to the
published 85-year figures. AML after 100 is not counted. Its weight is at most $S(100)/S(a)$, which is 4.0-4.7% for
$B=0.1$ and 2.6-6.8% across $B=0.085$-$0.115$; multiplying by the fraction of eligible people not
diagnosed by 100 bounds the omitted tail below 2.6 points in all three scenarios (same
assumptions; checked independently by yang2). Under this
mortality scenario the PPV declines with age much as published; the cut-off is a reasonable proxy,
slightly too generous before 60 and too strict at 70.

## 8. Re-tuning at the stated scale

Question (raised by the Pro review): does the scale mismatch of Corollary 3 change the model's
outputs, or only its labels? Set $N=5\times10^5$, $h=0.1$ (the stated $N\tau=10^5$ in ref-39 units)
and re-tune one other parameter to the observed lifetime risk (`runs/retune/`,
`runs/retune/compare_*.txt`, 200,200 lives each).

| run | lifetime risk | patterns lin/branch/late | median driver ages 1-4 | conditioned clone fitness | PPV 40/50/60/70 |
|---|---|---|---|---|---|
| published ($Nh=10^4$, $U=10^{-5}$, $s_1=0.16$) | 0.464% | 69/31/0% | 6.0, 39.1, 61.1, 74.1 | 0.18, 0.57, 1.54, 3.93 | 89, 72, 46, 21% |
| stated scale, $U=3.8\times10^{-6}$ | 0.496% | 66/34/0% | 5.4, 39.8, 61.9, 74.2 | 0.19, 0.60, 1.59, 4.12 | 100, 94, 61, 26% |
| stated scale, $s_1=0.1165$ | 0.486% | 57/43/0% | 5.9, 34.2, 58.8, 72.4 | 0.13, 0.41, 1.15, 3.11 | 100, 98, 82, 39% |

PPV numerators, denominators and two-sided 95% Clopper-Pearson intervals, rounded to whole percent
(ages 40/50/60/70; exact lower ends for 24/24 and 19/19 are 85.75% and 82.35%):
published 50/56 [78-96%], 217/303 [66-77%], 553/1204 [43-49%], 733/3422 [20-23%]; re-tuned in $U$
24/24 [86-100%], 133/142 [88-97%], 480/783 [58-65%], 770/2984 [24-27%]; re-tuned in $s_1$ 19/19
[82-100%], 117/120 [93-99%], 384/470 [78-85%], 665/1725 [36-41%]. With 19-24 qualifiers at 40, the
100% point estimates are not population certainties.

Incidence shape (ratio to the CRUK male-female average at 60-64, 70-74, 80-84): 0.41, 0.93, 1.50
published; 0.41, 0.90, 1.67 re-tuned in $U$; 0.60, 1.08, 1.55 re-tuned in $s_1$. These show that the
re-tuned runs approximately preserve the published simulation's shape, not that either fits CRUK
well. Re-tuning $U$ approximately preserves the listed outputs: incidence shape, patterns, driver
timing and the conditioned fitness ladder. Re-tuning $s_1$ also fits incidence but lowers the
absolute conditioned fitness by 21-28% (adjacent ratios approximately similar), makes more cases
branched and moves the later drivers 2-5 years earlier. In
both, the sweep PPVs change, because at the larger scale fewer people reach the sweep criterion at a
given age (19-24 and 117-142 qualifying at 40 and 50 against 56 and 303) and those who do are more
advanced. So the mismatch is an inconsistency between stated and simulated parameters: the listed
outputs are approximately preserved at the stated scale, while the PPVs and (depending on which
parameter absorbs the change) the absolute fitness scale depend on it. Only two one-parameter
re-tunings were tried; they do not establish unique parameters or rule out joint tuning, and they
cover only the listed outputs, not every qualitative conclusion of the paper.

## Files

- `sim/amlsim.cpp` (port of `sim_fast.py`), `sim/wlaw.cpp`, `sim/theory_branch.cpp`,
  `sim/theory_renewal.cpp` (a cruder renewal version of A3, kept for the record: 1.07% risk),
  `sim/theory_reduced.cpp`.
- `analysis/firstdriver_law.py`, `analysis/trunkstats.py`, `analysis/ppv.py`.
- Run outputs under `runs/` (not committed; regenerate with the commands in `notes/runs.md`).
