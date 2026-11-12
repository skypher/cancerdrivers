// theory_branch.cpp -- age-incidence of the sim_fast model from the branching recursion.
//
// Pure branching approximation (all mutant clones small, mean fitness 0):
//   an established clone of fitness F has n(t) = W e^{F t} cells, W ~ Exp(rate q(F)),
//   q(F) the survival probability of the Poisson(e^{F h}) Galton-Watson chain (h = dt);
//   daughters arrive at rate U n(t), increment s from the stage's DFE table, establish w.p. q(F+s).
//   G_4(a|F) = P(4-mutant reaches half the population within a) = exp(-q(F) N e^{-F a})
//              (logistic growth against wild type: f/(1-f) = W e^{Fa}/N).
//   G_k(a|F) = X/(q+X),  X(a) = U * int_0^a e^{F v} H_{k+1}(a-v; F) dv,
//              H_{k+1}(u; F) = E_s[ q(F+s) G_{k+1}(u | F+s) ]           (E over W ~ Exp(q) exact)
//   hazard h(T) = N U E_s[ q(s) G_1(T|s) ],  cumulative risk 1 - exp(-int h).
// The DFE tables are those of sim_fast.py (18 atoms per stage).
//
// Build: clang++ -O3 -march=native -std=c++17 theory_branch.cpp -o bin/theory_branch
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>

static double gammp(double a, double x) {
  if (x <= 0) return 0.0;
  double gln = std::lgamma(a);
  if (x < a + 1.0) {
    double ap = a, sum = 1.0 / a, del = sum;
    for (int n = 0; n < 10000; n++) { ap += 1.0; del *= x / ap; sum += del; if (std::fabs(del) < std::fabs(sum) * 1e-16) break; }
    return sum * std::exp(-x + a * std::log(x) - gln);
  }
  double b = x + 1.0 - a, c = 1.0 / 1e-300, d = 1.0 / b, h = d;
  for (int i = 1; i < 10000; i++) {
    double an = -i * (i - a); b += 2.0;
    d = an * d + b; if (std::fabs(d) < 1e-300) d = 1e-300;
    c = b + an / c; if (std::fabs(c) < 1e-300) c = 1e-300;
    d = 1.0 / d; double del = d * c; h *= del; if (std::fabs(del - 1.0) < 1e-16) break;
  }
  return 1.0 - std::exp(-x + a * std::log(x) - gln) * h;
}
static std::vector<double> dfe_table(double p, double sb) {
  const int nq = 1000; std::vector<double> d(nq + 1, -1.0); int prev = 0; double slim = 100.0 * sb;
  for (int i = 0; i <= 1000; i++) {
    double ts = slim * i / 1000.0; double cdf = gammp(1.0 / p, std::pow(ts / sb, p));
    int cur = (int)std::nearbyint(cdf * nq);
    for (int q = prev; q < cur; q++) d[q] = ts;
    prev = cur; if (cur == nq) { d[cur] = ts; break; }
  }
  return d;
}
// atoms: value -> probability, with sim_fast's index rule idx = round(u*1000) (end points half weight)
static std::map<long, double> atoms(const std::vector<double> &d, double unit) {
  std::map<long, double> m;
  for (int q = 0; q <= 1000; q++) {
    double w = (q == 0 || q == 1000) ? 0.5 / 1000.0 : 1.0 / 1000.0;
    m[std::lround(d[q] / unit)] += w;
  }
  return m;
}
static double qsurv(double F, double h) {  // 1 - q = exp(-m q), m = e^{F h}
  double m = std::exp(F * h);
  if (m <= 1.0) return 0.0;
  double q = 2.0 * (m - 1.0) / m;
  for (int it = 0; it < 200; it++) {
    double g = 1.0 - q - std::exp(-m * q), dg = -1.0 + m * std::exp(-m * q);
    double qn = q - g / dg; if (qn <= 0) qn = q / 2; if (std::fabs(qn - q) < 1e-15) { q = qn; break; } q = qn;
  }
  return q;
}

int main(int argc, char **argv) {
  double N = 1e5, h = 0.1, U = 1e-5, s1 = 0.16, r = 2.5, p = 3, years = 85, da = 0.05;
  std::string out = "theory_branch.tsv";
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      std::printf("theory_branch: branching-recursion hazard for the sim_fast model.\n"
                  "options: --N --dt --U --s1 --r --p --years --da (grid step, yr) --out file\n");
      return 0;
    }
    auto v = [&]() { return std::atof(argv[++i]); };
    if (a == "--N") N = v(); else if (a == "--dt") h = v(); else if (a == "--U") U = v();
    else if (a == "--s1") s1 = v(); else if (a == "--r") r = v(); else if (a == "--p") p = v();
    else if (a == "--years") years = v(); else if (a == "--da") da = v();
    else if (a == "--out") out = argv[++i];
  }
  const double unit = 0.002 * s1 / 0.16;  // all table values are multiples of 0.1*s1*r^k; r=2.5 -> gcd 0.0125 s1
  std::vector<std::map<long, double>> A(4);
  long maxu = 0;
  for (int k = 0; k < 4; k++) {
    auto d = dfe_table(p, s1 * std::pow(r, k));
    A[k] = atoms(d, unit);
    maxu += A[k].rbegin()->first;
  }
  const int na = (int)std::lround(years / da) + 1;
  // G[k][F index][a index]; only F values reachable matter, store dense over [0, maxu]
  std::vector<std::vector<double>> Gnext(maxu + 1), Gcur;
  // stage 4
  for (long fu = 0; fu <= maxu; fu++) {
    double F = fu * unit, q = qsurv(F, h);
    Gnext[fu].assign(na, 0.0);
    if (q <= 0) continue;
    for (int j = 0; j < na; j++) { double a = j * da; Gnext[fu][j] = std::exp(-q * N * std::exp(-F * a)); }
  }
  for (int k = 3; k >= 1; k--) {  // clone with k drivers; daughters draw from A[k]
    Gcur.assign(maxu + 1, std::vector<double>());
    for (long fu = 0; fu <= maxu; fu++) {
      double F = fu * unit, q = qsurv(F, h);
      Gcur[fu].assign(na, 0.0);
      if (q <= 0) continue;
      std::vector<double> H(na, 0.0);
      bool any = false;
      for (auto &at : A[k]) {
        long fd = fu + at.first; if (fd > maxu || Gnext[fd].empty()) continue;
        double qd = qsurv(fd * unit, h);
        for (int j = 0; j < na; j++) H[j] += at.second * qd * Gnext[fd][j];
        any = true;
      }
      if (!any) continue;
      // X(a) = U e^{F a} int_0^a e^{-F u} H(u) du  (trapezoid, running)
      double I = 0;
      for (int j = 0; j < na; j++) {
        double a = j * da;
        if (j > 0) I += 0.5 * da * (std::exp(-F * (a - da)) * H[j - 1] + std::exp(-F * a) * H[j]);
        double X = U * std::exp(F * a) * I;
        Gcur[fu][j] = X / (q + X);
      }
    }
    Gnext.swap(Gcur);
  }
  // hazard
  std::vector<double> haz(na, 0.0);
  for (auto &at : A[0]) {
    long fu = at.first; double q = qsurv(fu * unit, h);
    for (int j = 0; j < na; j++) haz[j] += N * U * at.second * q * Gnext[fu][j];
  }
  FILE *f = std::fopen(out.c_str(), "w");
  std::fprintf(f, "age\thazard_per100k\tcum_risk\n");
  double H = 0;
  for (int j = 0; j < na; j++) {
    if (j > 0) H += 0.5 * da * (haz[j - 1] + haz[j]);
    std::fprintf(f, "%.2f\t%.5f\t%.6f\n", j * da, 1e5 * haz[j], 1 - std::exp(-H));
  }
  std::fclose(f);
  double H85 = 0; for (int j = 1; j < na; j++) H85 += 0.5 * da * (haz[j - 1] + haz[j]);
  std::printf("branching recursion: cumulative risk to %.0f = %.4f%%; hazard at 60/70/80 = %.2f/%.2f/%.2f per 100k/yr\n",
              years, 100 * (1 - std::exp(-H85)), 1e5 * haz[(int)(60 / da)], 1e5 * haz[(int)(70 / da)], 1e5 * haz[(int)(80 / da)]);
  return 0;
}
