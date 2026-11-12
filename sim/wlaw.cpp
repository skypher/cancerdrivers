// wlaw.cpp -- law of the establishment factor in the sim_fast update n <- Poisson(n e^{g h}).
// For a lineage started from one cell with constant relative fitness g, W = lim n_t e^{-g t}.
// Reports the survival fraction (vs q(g)), and for survivors the mean and s.d. of log(q W),
// to be compared with log of an Exp(1) variable: mean -0.5772, s.d. 1.2825.
// Build: clang++ -O3 -march=native -fopenmp=libomp -std=c++17 wlaw.cpp -o bin/wlaw
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <omp.h>
static double qsurv(double F, double h) {
  double m = std::exp(F * h); double q = 2.0 * (m - 1.0) / m;
  for (int it = 0; it < 200; it++) { double g = 1.0 - q - std::exp(-m * q), dg = -1.0 + m * std::exp(-m * q); double qn = q - g / dg; if (std::fabs(qn - q) < 1e-15) { q = qn; break; } q = qn; }
  return q;
}
int main(int argc, char **argv) {
  if (argc > 1 && (!std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help"))) {
    std::printf("wlaw [lineages per g]: law of W for g in {0.05,0.16,0.5,1,2.5,4}, h = 0.1\n"); return 0;
  }
  long long L = argc > 1 ? std::atoll(argv[1]) : 200000;
  const double h = 0.1, gs[6] = {0.05, 0.16, 0.5, 1.0, 2.5, 4.0};
  for (double g : gs) {
    double q = qsurv(g, h), m = std::exp(g * h);
    long long surv = 0; double s1 = 0, s2 = 0;
#pragma omp parallel for reduction(+ : surv, s1, s2) schedule(dynamic, 1000)
    for (long long i = 0; i < L; i++) {
      std::mt19937_64 rg(1234567ULL + 1000003ULL * i + (unsigned long long)(g * 1e6));
      double n = 1; int t = 0;
      while (n > 0 && n < 1e7) { std::poisson_distribution<long long> P(n * m); n = (double)P(rg); t++; }
      if (n > 0) { double W = n * std::exp(-g * h * t); double lw = std::log(q * W); surv++; s1 += lw; s2 += lw * lw; }
    }
    double mu = s1 / surv, sd = std::sqrt(s2 / surv - mu * mu);
    std::printf("g %.2f: survival %.5f (q %.5f, 2gh %.5f); log(qW): mean %+.4f (Exp: -0.5772), sd %.4f (Exp: 1.2825)\n",
                g, (double)surv / L, q, 2 * g * h, mu, sd);
    std::fflush(stdout);
  }
  return 0;
}
