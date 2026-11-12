// theory_reduced.cpp -- reduced single-lineage model of the sim_fast dynamics.
//
// Ingredients (each switchable) used to explain the incidence of the full model:
//  * drift only at establishment: a clone founded when the mean fitness is xbar, with fitness F,
//    survives w.p. q(F - xbar) (Poisson Galton-Watson, generation dt) and then has W ~ Exp(q) cells
//    at its founding in the deterministic-growth representation;
//  * deterministic growth: unnormalised size w_j e^{F_j (u - t_j)}, frequencies normalised by
//    N + sum of unnormalised sizes (exact mean dynamics of the update n <- n e^{(F-xbar)dt});
//  * one AML lineage: only its newest clone produces daughters (rate U N f, increment from the
//    stage's DFE table, establishment w.p. q(F + s - xbar));
//  * optional inert background: first-driver clones of the same person (Poisson rate Lambda1,
//    uniform ages) that grow and compete but do not mutate further (--background 1).
// AML when the lineage's 4-driver clone exceeds half the population.  The founding age of the
// lineage is uniform on [0, years]; incidence(age) = Lambda1 * years * P(AML at age) per sample.
//
// Build: clang++ -O3 -march=native -fopenmp=libomp -std=c++17 theory_reduced.cpp -o bin/theory_reduced
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <omp.h>

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
static const double QX = 40.0, QD = 0.0005;
static std::vector<double> QT;
static inline double qtab(double x) {
  if (x <= 0) return 0.0;
  double z = x / QD; size_t i = (size_t)z;
  if (i + 1 >= QT.size()) return QT.back();
  double fr = z - i; return QT[i] * (1 - fr) + QT[i + 1] * fr;
}
struct DFE { std::vector<double> s, p; };
struct Cl { double F, w, t; int depth; };

int main(int argc, char **argv) {
  double N = 1e5, h = 0.1, U = 1e-5, s1 = 0.16, r = 2.5, p = 3, years = 85, du = 0.05;
  long long M = 1000000; int threads = 8, background = 1; std::string out = "reduced.tsv";
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      std::printf("theory_reduced: reduced single-lineage model of the sim_fast incidence.\n"
                  "options: --N --dt --U --s1 --r --p --years --du --samples M --threads T --background 0|1 --out file\n");
      return 0;
    }
    auto v = [&]() { return std::atof(argv[++i]); };
    if (a == "--N") N = v(); else if (a == "--dt") h = v(); else if (a == "--U") U = v();
    else if (a == "--s1") s1 = v(); else if (a == "--r") r = v(); else if (a == "--p") p = v();
    else if (a == "--years") years = v(); else if (a == "--du") du = v();
    else if (a == "--samples") M = (long long)v(); else if (a == "--threads") threads = (int)v();
    else if (a == "--background") background = (int)v(); else if (a == "--out") out = argv[++i];
  }
  for (double x = 0; x <= QX + QD; x += QD) QT.push_back(qsurv(x, h));
  std::vector<DFE> D(4);
  for (int k = 0; k < 4; k++) {
    auto d = dfe_table(p, s1 * std::pow(r, k));
    std::map<double, double> m;
    for (int q = 0; q <= 1000; q++) m[d[q]] += (q == 0 || q == 1000) ? 0.0005 : 0.001;
    for (auto &z : m) { D[k].s.push_back(z.first); D[k].p.push_back(z.second); }
  }
  double Q1 = 0; for (size_t i = 0; i < D[0].s.size(); i++) Q1 += D[0].p[i] * qtab(D[0].s[i]);
  double Lambda1 = N * U * Q1;
  const int nb = (int)std::lround(years);  // 1-year bins
  std::vector<double> hist(nb, 0.0);
  std::vector<double> occsum(4, 0.0); double ncase = 0, nlate = 0;
  long long prog = std::max<long long>(1, M / 10);
  time_t tstart = time(nullptr);
  omp_set_num_threads(threads);
#pragma omp parallel
  {
    std::vector<double> lh(nb, 0.0), locc(4, 0.0); double lc = 0, ll = 0;
    std::mt19937_64 g(424242 + 104729ULL * omp_get_thread_num());
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    auto expo = [&](double rate) { return -std::log(1.0 - U01(g)) / rate; };
    auto pickw = [&](const std::vector<double> &w) {
      double tot = 0; for (double x : w) tot += x;
      double u = U01(g) * tot; for (size_t i = 0; i < w.size(); i++) { u -= w[i]; if (u <= 0) return (int)i; }
      return (int)w.size() - 1;
    };
    std::vector<double> w0(D[0].s.size());
    for (size_t i = 0; i < w0.size(); i++) w0[i] = D[0].p[i] * qtab(D[0].s[i]);
#pragma omp for schedule(dynamic, 2000)
    for (long long it = 0; it < M; it++) {
      if (it % prog == 0) {
        time_t now = time(nullptr); char buf[32]; std::strftime(buf, sizeof buf, "%H:%M:%SZ", std::gmtime(&now));
#pragma omp critical
        { std::printf("[%s] sample block at %lld/%lld elapsed %lds\n", buf, it, M, (long)(now - tstart)); std::fflush(stdout); }
      }
      double u1 = U01(g) * years;
      double F1 = D[0].s[pickw(w0)];
      double W1 = expo(qtab(F1));
      // background clones: (age, F, W)
      std::vector<Cl> bg; std::vector<double> bgW, bgT;
      if (background) {
        std::poisson_distribution<int> P(Lambda1 * years);
        int K = P(g);
        for (int k = 0; k < K; k++) { double F = D[0].s[pickw(w0)]; bgT.push_back(U01(g) * years); bgW.push_back(expo(qtab(F))); bg.push_back({F, -1.0, 0.0, 0}); }
      }
      std::vector<Cl> tr;  // trunk clones
      double E = expo(1.0), cum = 0;
      double occ[4] = {u1, -1, -1, -1};
      double tdx = -1;
      for (double u = 0; u < years; u += du) {
        // births of background clones and of the trunk founder in this step
        double Z = N;
        for (auto &c : tr) Z += c.w * std::exp(c.F * (u - c.t));
        for (auto &c : bg) if (c.w > 0) Z += c.w * std::exp(c.F * (u - c.t));
        for (size_t k = 0; k < bg.size(); k++) if (bg[k].w < 0 && bgT[k] <= u) { bg[k].w = bgW[k] * Z / N; bg[k].t = u; }
        if (tr.empty() && u1 <= u) { tr.push_back({F1, W1 * Z / N, u, 1}); }
        if (tr.empty()) continue;
        // frequencies
        Z = N; double xs = 0;
        std::vector<double> nt(tr.size());
        for (size_t k = 0; k < tr.size(); k++) { nt[k] = tr[k].w * std::exp(tr[k].F * (u - tr[k].t)); Z += nt[k]; xs += tr[k].F * nt[k]; }
        for (auto &c : bg) if (c.w > 0) { double n = c.w * std::exp(c.F * (u - c.t)); Z += n; xs += c.F * n; }
        double xbar = xs / Z;
        int d = (int)tr.size();
        double ftip = nt[d - 1] / Z;
        if (d == 4) { if (ftip > 0.5) { tdx = u; break; } continue; }
        const DFE &Dk = D[d];
        double lam = 0;
        for (size_t i = 0; i < Dk.s.size(); i++) lam += Dk.p[i] * qtab(tr[d - 1].F + Dk.s[i] - xbar);
        cum += U * N * ftip * lam * du;
        if (cum >= E) {
          std::vector<double> w(Dk.s.size());
          for (size_t i = 0; i < w.size(); i++) w[i] = Dk.p[i] * qtab(tr[d - 1].F + Dk.s[i] - xbar);
          double Fn = tr[d - 1].F + Dk.s[pickw(w)];
          double Wn = expo(qtab(Fn - xbar));
          tr.push_back({Fn, Wn * Z / N, u, d + 1});
          occ[d] = u;
          cum = 0; E = expo(1.0);
        }
      }
      if (tdx >= 0) {
        int b = (int)tdx; if (b < nb) lh[b] += 1.0;
        for (int k = 0; k < 4; k++) locc[k] += occ[k];
        lc += 1;
        if (u1 <= 20) ll += 1;
      }
    }
#pragma omp critical
    { for (int b = 0; b < nb; b++) hist[b] += lh[b]; for (int k = 0; k < 4; k++) occsum[k] += locc[k]; ncase += lc; nlate += ll; }
  }
  FILE *f = std::fopen(out.c_str(), "w");
  std::fprintf(f, "# Lambda1 %.6g, samples %lld, background %d\nage\tincidence_per100k\n", Lambda1, M, background);
  double cumr = 0;
  for (int b = 0; b < nb; b++) {
    double inc = Lambda1 * years * hist[b] / M;  // per year (1-yr bins)
    cumr += inc;
    std::fprintf(f, "%d\t%.4f\n", b, 1e5 * inc);
  }
  std::fclose(f);
  std::printf("reduced (background %d): Lambda1 %.5f; cumulative risk to %.0f = %.4f%%; mean trunk ages %.1f %.1f %.1f %.1f; "
              "P(first driver <= 20 | AML) = %.3f; cases in sample %.0f\n",
              background, Lambda1, years, 100 * (1 - std::exp(-cumr)), occsum[0] / ncase, occsum[1] / ncase, occsum[2] / ncase,
              occsum[3] / ncase, nlate / ncase, ncase);
  return 0;
}
