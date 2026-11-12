// theory_renewal.cpp -- two-regime renewal approximation R for the sim_fast model.
//
// Single AML lineage; drift only at establishment; deterministic logistic growth afterwards.
//  1. first driver: F1 ~ p0(s) q(s) / Q1 (established first drivers), Lambda1 = N U Q1.
//  2. second driver in the branching regime: P(tau2 > t | W1) = exp(-U W1 Q2 (e^{F1 t}-1)/F1),
//     W1 ~ Exp(q(F1)); s2 ~ p1(s) q(F1+s).
//  3. clone 2 sweeps against wild type: f2(u) = 1/(1+exp(-F2 (u - t2 - sig2))), sig2 = ln(N/W2)/F2,
//     W2 ~ Exp(q(F2)).
//  4. third driver: intensity U N f2(u) sum_s p2(s) q(F2 + s - F2 f2(u)); s3 drawn at the arrival.
//  5. clone 3 sweeps against the clone-2 background with relative rate g3 = s3 + F2 (1 - f2(t3)):
//     f3(u) = 1/(1+exp(-g3 (u - t3 - sig3))), sig3 = ln(N/W3)/g3, W3 ~ Exp(q(g3)).
//  6. fourth driver: intensity U N f3(u) sum_s p3(s) q(s + g3 (1 - f3(u))).
//  7. AML when clone 4 reaches 1/2: t4 + ln(N/W4)/g4, g4 = s4 + g3 (1 - f3(t4)), W4 ~ Exp(q(g4)).
// Output: hazard h(T) = Lambda1 D(T), D the CDF of the completion time (Monte Carlo over lineages).
//
// Build: clang++ -O3 -march=native -fopenmp=libomp -std=c++17 theory_renewal.cpp -o bin/theory_renewal
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <omp.h>
#include <ctime>

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
static double qsurv(double F, double h) {
  double m = std::exp(F * h);
  if (m <= 1.0) return 0.0;
  double q = 2.0 * (m - 1.0) / m;
  for (int it = 0; it < 200; it++) {
    double g = 1.0 - q - std::exp(-m * q), dg = -1.0 + m * std::exp(-m * q);
    double qn = q - g / dg; if (qn <= 0) qn = q / 2; if (std::fabs(qn - q) < 1e-15) { q = qn; break; } q = qn;
  }
  return q;
}
struct DFE { std::vector<double> s, p; };
// tabulated q(x) for x in [0, QX], step QD, linear interpolation; q = 0 for x <= 0
static const double QX = 40.0, QD = 0.0005;
static std::vector<double> QT;
static inline double qtab(double x) {
  if (x <= 0) return 0.0;
  double z = x / QD; size_t i = (size_t)z;
  if (i + 1 >= QT.size()) return QT.back();
  double fr = z - i; return QT[i] * (1 - fr) + QT[i + 1] * fr;
}

int main(int argc, char **argv) {
  double N = 1e5, h = 0.1, U = 1e-5, s1 = 0.16, r = 2.5, p = 3, years = 85;
  long long M = 2000000; int threads = 8; std::string out = "renewal.tsv";
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      std::printf("theory_renewal: two-regime renewal approximation of the sim_fast incidence.\n"
                  "options: --N --dt --U --s1 --r --p --years --lineages M --threads T --out file\n");
      return 0;
    }
    auto v = [&]() { return std::atof(argv[++i]); };
    if (a == "--N") N = v(); else if (a == "--dt") h = v(); else if (a == "--U") U = v();
    else if (a == "--s1") s1 = v(); else if (a == "--r") r = v(); else if (a == "--p") p = v();
    else if (a == "--years") years = v(); else if (a == "--lineages") M = (long long)v();
    else if (a == "--threads") threads = (int)v(); else if (a == "--out") out = argv[++i];
  }
  std::vector<DFE> D(4);
  for (int k = 0; k < 4; k++) {
    auto d = dfe_table(p, s1 * std::pow(r, k));
    std::map<double, double> m;
    for (int q = 0; q <= 1000; q++) m[d[q]] += (q == 0 || q == 1000) ? 0.0005 : 0.001;
    for (auto &z : m) { D[k].s.push_back(z.first); D[k].p.push_back(z.second); }
  }
  for (double x = 0; x <= QX + QD; x += QD) QT.push_back(qsurv(x, h));
  // stage-1 law
  double Q1 = 0; for (size_t i = 0; i < D[0].s.size(); i++) Q1 += D[0].p[i] * qsurv(D[0].s[i], h);
  double Lambda1 = N * U * Q1;
  const double du = 0.05; const int nb = (int)std::lround(years / 0.5);
  std::vector<double> hist(nb, 0.0);
  std::vector<double> stage(6, 0.0); long long ndone = 0;
  omp_set_num_threads(threads);
#pragma omp parallel
  {
    std::vector<double> lh(nb, 0.0), lst(6, 0.0); long long lnd = 0;
    std::mt19937_64 g(9876 + 7919ULL * omp_get_thread_num());
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    auto pick = [&](const DFE &d, const std::vector<double> &w) {
      double tot = 0; for (double x : w) tot += x;
      double u = U01(g) * tot; for (size_t i = 0; i < w.size(); i++) { u -= w[i]; if (u <= 0) return (int)i; }
      return (int)w.size() - 1;
    };
    auto expo = [&](double rate) { return -std::log(1.0 - U01(g)) / rate; };
#pragma omp for schedule(dynamic, 1000)
    for (long long it = 0; it < M; it++) {
      if (it % (M / 10 > 0 ? M / 10 : 1) == 0) {
        time_t now = time(nullptr); char buf[32]; std::strftime(buf, sizeof buf, "%H:%M:%SZ", std::gmtime(&now));
#pragma omp critical
        { std::printf("[%s] lineage block starting at %lld/%lld\n", buf, it, M); std::fflush(stdout); }
      }
      // 1
      std::vector<double> w(D[0].s.size());
      for (size_t i = 0; i < w.size(); i++) w[i] = D[0].p[i] * qtab(D[0].s[i]);
      double F1 = D[0].s[pick(D[0], w)];
      double W1 = expo(qtab(F1));
      // 2
      std::vector<double> w2(D[1].s.size()); double Q2 = 0;
      for (size_t i = 0; i < w2.size(); i++) { w2[i] = D[1].p[i] * qtab(F1 + D[1].s[i]); Q2 += w2[i]; }
      double E = expo(1.0);  // U W1 Q2 (e^{F1 t}-1)/F1 = E
      double t2 = std::log1p(E * F1 / (U * W1 * Q2)) / F1;
      if (!(t2 < years)) continue;
      double F2 = F1 + D[1].s[pick(D[1], w2)];
      double W2 = expo(qtab(F2));
      double sig2 = std::log(N / W2) / F2;
      // 4: third driver
      double t3 = -1, cum = 0, E3 = expo(1.0);
      for (double u = t2; u < years; u += du) {
        double f2 = 1.0 / (1.0 + std::exp(-F2 * (u + 0.5 * du - t2 - sig2)));
        double lam = 0;
        for (size_t i = 0; i < D[2].s.size(); i++) lam += D[2].p[i] * qtab(F2 + D[2].s[i] - F2 * f2);
        cum += U * N * f2 * lam * du;
        if (cum >= E3) { t3 = u + du; break; }
      }
      if (t3 < 0) continue;
      double f2t3 = 1.0 / (1.0 + std::exp(-F2 * (t3 - t2 - sig2)));
      std::vector<double> w3(D[2].s.size());
      for (size_t i = 0; i < w3.size(); i++) w3[i] = D[2].p[i] * qtab(F2 + D[2].s[i] - F2 * f2t3);
      double s3 = D[2].s[pick(D[2], w3)];
      double g3 = s3 + F2 * (1 - f2t3);
      double W3 = expo(qtab(g3));
      double sig3 = std::log(N / W3) / g3;
      // 6: fourth driver
      double t4 = -1; cum = 0; double E4 = expo(1.0);
      for (double u = t3; u < years; u += du) {
        double f3 = 1.0 / (1.0 + std::exp(-g3 * (u + 0.5 * du - t3 - sig3)));
        double lam = 0;
        for (size_t i = 0; i < D[3].s.size(); i++) lam += D[3].p[i] * qtab(D[3].s[i] + g3 * (1 - f3));
        cum += U * N * f3 * lam * du;
        if (cum >= E4) { t4 = u + du; break; }
      }
      if (t4 < 0) continue;
      double f3t4 = 1.0 / (1.0 + std::exp(-g3 * (t4 - t3 - sig3)));
      std::vector<double> w4(D[3].s.size());
      for (size_t i = 0; i < w4.size(); i++) w4[i] = D[3].p[i] * qtab(D[3].s[i] + g3 * (1 - f3t4));
      double s4 = D[3].s[pick(D[3], w4)];
      double g4 = s4 + g3 * (1 - f3t4);
      double W4 = expo(qtab(g4));
      double tdx = t4 + std::log(N / W4) / g4;
      if (tdx < years) {
        int b = (int)(tdx / 0.5); if (b < nb) lh[b] += 1.0;
        lst[0] += t2; lst[1] += t3 - t2; lst[2] += t4 - t3; lst[3] += tdx - t4; lnd++;
      }
    }
#pragma omp critical
    { for (int b = 0; b < nb; b++) hist[b] += lh[b]; for (int k = 0; k < 6; k++) stage[k] += lst[k]; ndone += lnd; }
  }
  FILE *f = std::fopen(out.c_str(), "w");
  std::fprintf(f, "# Lambda1 = %.6g per yr, lineages %lld\nage\tD\thazard_per100k\n", Lambda1, M);
  double cdf = 0, H = 0;
  for (int b = 0; b < nb; b++) {
    cdf += hist[b] / M;
    double age = (b + 1) * 0.5;
    std::fprintf(f, "%.1f\t%.6g\t%.5f\n", age, cdf, 1e5 * Lambda1 * cdf);
  }
  std::fclose(f);
  // cumulative risk: int_0^85 Lambda1 D(T) dT
  cdf = 0; for (int b = 0; b < nb; b++) { cdf += hist[b] / M; H += Lambda1 * cdf * 0.5; }
  std::printf("renewal: Lambda1 = %.5f/yr; D(85) = %.5f; cumulative risk to %.0f = %.4f%%; "
              "mean stage durations (completed lineages): t2 %.1f, t3-t2 %.1f, t4-t3 %.1f, dx-t4 %.1f\n",
              Lambda1, cdf, years, 100 * (1 - std::exp(-H)), stage[0] / ndone, stage[1] / ndone, stage[2] / ndone, stage[3] / ndone);
  return 0;
}
