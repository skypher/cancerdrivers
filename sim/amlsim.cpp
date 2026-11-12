// amlsim.cpp -- C++ port of the stochastic clone model in sim_fast.py
// (the-blundell-lab/preAML_evolutionary_dynamics, commit 5a53fe7), used for
// Fig. 5 of Watson et al., Nature 2026.  Same update rule, same DFE table,
// same AML criterion, same pattern classification; only the random stream
// differs (std::mt19937_64 seeded per simulation instead of numpy/random).
//
// Build: clang++ -O3 -march=native -fopenmp=libomp -std=c++17 amlsim.cpp -o bin/amlsim
// Run:   bin/amlsim --help
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <atomic>
#include <omp.h>

struct Params {
  long long sims = 1000;
  unsigned long long seed = 1000;
  int threads = 4;
  double N = 1e5;
  double dt = 0.1;
  double years = 85.0;
  double U = 1e-5;
  double s1 = 0.16;
  double r = 2.5;
  double p = 3.0;
  int maxmut = 4;
  double aml_frac = 0.5;
  std::string out = "amlsim";
  bool dump_cases = false;
  int dump_max = 1000000;
  long long progress_every = 0;  // 0: automatic
};

static void usage() {
  std::printf(
      "amlsim: C++ port of sim_fast.py (Watson et al. Nature 2026, Fig. 5 model).\n"
      "Options (defaults = published parameters):\n"
      "  --sims N        number of simulated lives (1000)\n"
      "  --seed S        base seed (1000); life i uses splitmix(S, i)\n"
      "  --threads T     OpenMP threads (4)\n"
      "  --N x           initial wild-type size (1e5)\n"
      "  --dt x          time step in years (0.1)\n"
      "  --years x       lifespan simulated (85)\n"
      "  --U x           driver rate per cell per year (1e-5)\n"
      "  --s1 x          DFE scale of first drivers (0.16)\n"
      "  --r x           DFE scale ratio per extra driver (2.5)\n"
      "  --p x           stretched-exponential exponent (3)\n"
      "  --out prefix    output prefix (amlsim)\n"
      "  --dump-cases    write every AML case's clones to <prefix>_cases.bin\n"
      "  --progress K    progress line every K lives (auto)\n"
      "Outputs: <prefix>_summary.txt, <prefix>_cases.tsv, <prefix>_ppv.tsv, optional _cases.bin\n");
}

static uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

// regularized lower incomplete gamma P(a,x) (series / continued fraction)
static double gammp(double a, double x) {
  if (x <= 0) return 0.0;
  double gln = std::lgamma(a);
  if (x < a + 1.0) {
    double ap = a, sum = 1.0 / a, del = sum;
    for (int n = 0; n < 10000; n++) {
      ap += 1.0; del *= x / ap; sum += del;
      if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
    }
    return sum * std::exp(-x + a * std::log(x) - gln);
  }
  double b = x + 1.0 - a, c = 1.0 / 1e-300, d = 1.0 / b, h = d;
  for (int i = 1; i < 10000; i++) {
    double an = -i * (i - a);
    b += 2.0;
    d = an * d + b; if (std::fabs(d) < 1e-300) d = 1e-300;
    c = b + an / c; if (std::fabs(c) < 1e-300) c = 1e-300;
    d = 1.0 / d;
    double del = d * c; h *= del;
    if (std::fabs(del - 1.0) < 1e-16) break;
  }
  return 1.0 - std::exp(-x + a * std::log(x) - gln) * h;
}

// Same table as mutation_fitness_generator(p, sb, offset=0): grid 0..100 sb in 1000 steps,
// d[q] = grid value at which round(1000*CDF) first exceeds q.
static std::vector<double> dfe_table(double p, double sb) {
  const int nq = 1000;
  std::vector<double> d(nq + 1, -1.0);
  int prev = 0;
  double slim = 100.0 * sb;
  for (int i = 0; i <= 1000; i++) {
    double ts = slim * i / 1000.0;
    double cdf = gammp(1.0 / p, std::pow(ts / sb, p)) / gammp(1.0 / p, std::pow(100.0, p));  // normalised over [0, 100 sb] as in sim_fast
    int cur = (int)std::nearbyint(cdf * nq);
    for (int q = prev; q < cur; q++) d[q] = ts;
    prev = cur;
    if (cur == nq) { d[cur] = ts; break; }
  }
  return d;
}

struct Clone {
  long long size;
  double fit;
  int nmut;
  int parent;
  int muts[4];
  double incs[4];
  double occ;
  bool has_children;
  bool purged;
  bool aml;
  float maxtraj;
  std::vector<float> traj;
};

struct CaseOut {
  long long sim;
  double dx;
  int pattern;  // 0 linear 1 branched 2 late
  double occ[4], inc[4], fit[4];
  double first_detect;
  std::vector<Clone> clones;  // only when dumping
};

static inline long long poisson(std::mt19937_64 &g, double mean) {
  if (mean <= 0) return 0;
  std::poisson_distribution<long long> P(mean);
  return P(g);
}

int main(int argc, char **argv) {
  Params P;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto nxt = [&](void) -> const char * {
      if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", a.c_str()); std::exit(2); }
      return argv[++i];
    };
    if (a == "-h" || a == "--help") { usage(); return 0; }
    else if (a == "--sims") P.sims = std::atoll(nxt());
    else if (a == "--seed") P.seed = std::strtoull(nxt(), nullptr, 10);
    else if (a == "--threads") P.threads = std::atoi(nxt());
    else if (a == "--N") P.N = std::atof(nxt());
    else if (a == "--dt") P.dt = std::atof(nxt());
    else if (a == "--years") P.years = std::atof(nxt());
    else if (a == "--U") P.U = std::atof(nxt());
    else if (a == "--s1") P.s1 = std::atof(nxt());
    else if (a == "--r") P.r = std::atof(nxt());
    else if (a == "--p") P.p = std::atof(nxt());
    else if (a == "--out") P.out = nxt();
    else if (a == "--dump-cases") P.dump_cases = true;
    else if (a == "--progress") P.progress_every = std::atoll(nxt());
    else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); usage(); return 2; }
  }
  const int T = (int)std::lround(P.years / P.dt);
  const double dt = P.dt;
  std::vector<std::vector<double>> DFE(4);
  for (int k = 0; k < 4; k++) DFE[k] = dfe_table(P.p, P.s1 * std::pow(P.r, k));
  for (int k = 0; k < 4; k++) {
    double m = 0; for (int q = 0; q <= 1000; q++) m += DFE[k][q];
    std::printf("DFE[%d] scale %.4f: min %.4f max %.4f table-mean %.4f\n", k, P.s1 * std::pow(P.r, k),
                DFE[k][0], DFE[k][1000], m / 1001.0);
  }
  const int AIDX[4] = {(int)std::lround(40 / dt), (int)std::lround(50 / dt), (int)std::lround(60 / dt), (int)std::lround(70 / dt)};
  // per-life PPV record: largest, second VAF at ages 40..70 (-1 if already diagnosed)
  std::vector<float> ppvrec((size_t)P.sims * 8, -1.0f);
  std::vector<double> dxage((size_t)P.sims, -1.0);
  std::vector<CaseOut> cases;
  long long prog = P.progress_every > 0 ? P.progress_every : std::max<long long>(1, P.sims / 20);
  std::atomic<long long> done(0), ncase(0);
  time_t t_start = time(nullptr);
  omp_set_num_threads(P.threads);
#pragma omp parallel
  {
    std::vector<CaseOut> mycases;
#pragma omp for schedule(dynamic, 64)
    for (long long sim = 0; sim < P.sims; sim++) {
      std::mt19937_64 g(splitmix64(P.seed * 0x100000001B3ULL + (uint64_t)sim));
      std::uniform_real_distribution<double> U01(0.0, 1.0);
      std::vector<Clone> C;
      C.reserve(256);
      {
        Clone root{};
        root.size = (long long)P.N; root.fit = 0; root.nmut = 0; root.parent = -1; root.occ = 0;
        root.has_children = false; root.purged = false; root.aml = false;
        root.traj.assign(T, 0.0f); root.traj[0] = (float)P.N; root.maxtraj = (float)P.N;
        C.push_back(std::move(root));
      }
      bool aml = false; int tdx = -1;
      for (int t = 0; t < T; t++) {
        // mutate (snapshot: only clones existing at the start of the step)
        size_t nb = C.size();
        for (size_t i = 0; i < nb; i++) {
          if (C[i].purged) continue;
          if (C[i].nmut < P.maxmut && C[i].size > 0) {
            long long nc = poisson(g, (double)C[i].size * P.U * dt);
            for (long long c = 0; c < nc; c++) {
              int k = C[i].nmut;
              int idx = (int)std::nearbyint(U01(g) * 1000.0);
              double s = DFE[k][idx];
              Clone ch{};
              ch.size = 1; ch.nmut = k + 1; ch.parent = (int)i; ch.occ = t * dt;
              for (int j = 0; j < k; j++) { ch.muts[j] = C[i].muts[j]; ch.incs[j] = C[i].incs[j]; }
              ch.muts[k] = (int)C.size(); ch.incs[k] = s;
              ch.fit = C[i].fit + s;
              ch.has_children = false; ch.purged = false; ch.aml = false;
              ch.traj.assign(T, 0.0f); ch.traj[t] = 1.0f; ch.maxtraj = 1.0f;
              C[i].has_children = true;
              C.push_back(std::move(ch));
            }
          }
        }
        // mean fitness
        double tot = 0, mf = 0;
        for (auto &c : C) if (!c.purged) tot += (double)c.size;
        for (auto &c : C) if (!c.purged) mf += c.fit * (double)c.size / tot;
        // select
        for (auto &c : C) {
          if (c.purged) continue;
          double ne = (double)c.size * std::exp((c.fit - mf) * dt);
          long long ns = ne > 0 ? poisson(g, ne) : 0;
          c.size = ns; c.traj[t] = (float)ns;
          if ((float)ns > c.maxtraj) c.maxtraj = (float)ns;
        }
        // purge
        for (auto &c : C) {
          if (c.purged) continue;
          if (c.size == 0 && !c.has_children && c.maxtraj < 10.0f) {
            c.purged = true; std::vector<float>().swap(c.traj);
          }
        }
        // PPV record at ages 40..70 (before the diagnosis check, as in fig5e which
        // reads the trajectory at tA; a life diagnosed at tA has pop > 0 there too)
        for (int a = 0; a < 4; a++) {
          if (t != AIDX[a]) continue;
          double pop = 0;
          for (auto &c : C) if (!c.purged) pop += (double)c.size;
          // variant cell counts
          std::vector<std::pair<int, double>> mut;
          for (auto &c : C) {
            if (c.purged || c.size <= 0) continue;
            for (int j = 0; j < c.nmut; j++) {
              int m = c.muts[j]; bool f = false;
              for (auto &pr : mut) if (pr.first == m) { pr.second += (double)c.size; f = true; break; }
              if (!f) mut.push_back({m, (double)c.size});
            }
          }
          double v1 = 0, v2 = 0;
          for (auto &pr : mut) {
            double v = 0.5 * pr.second / pop;
            if (v > v1) { v2 = v1; v1 = v; } else if (v > v2) v2 = v;
          }
          ppvrec[(size_t)sim * 8 + 2 * a] = (float)v1;
          ppvrec[(size_t)sim * 8 + 2 * a + 1] = (float)v2;
        }
        // AML diagnosis
        double tot2 = 0;
        for (auto &c : C) if (!c.purged) tot2 += (double)c.size;
        for (auto &c : C) {
          if (c.purged) continue;
          if (c.nmut == 4 && (double)c.size / tot2 > P.aml_frac) { aml = true; c.aml = true; }
        }
        if (aml) { tdx = t; break; }
      }
      if (aml) {
        dxage[sim] = tdx * dt;
        CaseOut co{};
        co.sim = sim; co.dx = tdx * dt;
        int amlc = -1;
        for (size_t i = 0; i < C.size(); i++) if (!C[i].purged && C[i].aml) { amlc = (int)i; break; }
        const Clone &A = C[amlc];
        // trunk clones: ids are the mutation ids
        for (int j = 0; j < 4; j++) {
          const Clone &K = C[A.muts[j]];
          co.occ[j] = K.occ; co.inc[j] = A.incs[j]; co.fit[j] = K.fit;
        }
        // classification exactly as figure5_pipeline.classify
        std::vector<double> pop(T, 0.1);
        for (auto &c : C) if (!c.purged) for (int t = 0; t < T; t++) pop[t] += c.traj[t];
        int m1 = A.muts[0];
        std::vector<double> cf(T, 0.0);
        for (auto &c : C) {
          if (c.purged) continue;
          bool has = false;
          for (int j = 0; j < c.nmut; j++) if (c.muts[j] == m1) has = true;
          if (has) for (int t = 0; t < T; t++) cf[t] += c.traj[t];
        }
        int fd = 0;
        for (int t = 0; t < T; t++) if (cf[t] / pop[t] >= 0.002) { fd = t; break; }
        co.first_detect = fd * dt;
        int pattern = 0;
        if (co.dx - co.first_detect <= 2.0) pattern = 2;
        else {
          for (size_t i = 0; i < C.size() && pattern == 0; i++) {
            const Clone &c = C[i];
            if (c.purged || c.nmut == 0) continue;
            bool intrunk = false;
            for (int j = 0; j < 4; j++) if ((int)i == A.muts[j]) intrunk = true;
            if (intrunk) continue;
            for (int t = 0; t < T; t++) if (c.traj[t] / pop[t] > 0.10) { pattern = 1; break; }
          }
        }
        co.pattern = pattern;
        if (P.dump_cases) {
          for (size_t i = 0; i < C.size(); i++) {
            if (C[i].purged) continue;
            if (C[i].maxtraj < 10.0f) continue;
            co.clones.push_back(C[i]);
          }
        }
        mycases.push_back(std::move(co));
        ncase++;
      }
      long long d = ++done;
      if (d % prog == 0) {
        time_t now = time(nullptr);
        char buf[32]; std::strftime(buf, sizeof buf, "%H:%M:%SZ", std::gmtime(&now));
        std::printf("[%s] lives %lld/%lld (%.1f%%) cases %lld elapsed %lds\n", buf, d, P.sims,
                    100.0 * d / P.sims, (long long)ncase, (long)(now - t_start));
        std::fflush(stdout);
      }
    }
#pragma omp critical
    for (auto &c : mycases) cases.push_back(std::move(c));
  }
  std::sort(cases.begin(), cases.end(), [](const CaseOut &a, const CaseOut &b) { return a.sim < b.sim; });
  // summary
  int cnt[3] = {0, 0, 0};
  std::vector<long long> inc(17, 0);
  for (auto &c : cases) {
    cnt[c.pattern]++;
    int b = (int)(c.dx / 5.0); if (b > 16) b = 16;
    inc[b]++;
  }
  std::string fs = P.out + "_summary.txt";
  FILE *f = std::fopen(fs.c_str(), "w");
  std::fprintf(f, "sims %lld\ncases %zu\nfraction %.6f\nlinear %d\nbranched %d\nlate %d\n", P.sims, cases.size(),
               (double)cases.size() / P.sims, cnt[0], cnt[1], cnt[2]);
  std::fprintf(f, "params N %g dt %g U %g s1 %g r %g p %g years %g seed %llu\n", P.N, P.dt, P.U, P.s1, P.r, P.p,
               P.years, P.seed);
  std::fprintf(f, "incidence_per100k_per_yr");
  for (int b = 0; b < 17; b++) std::fprintf(f, " %.3f", 1e5 * inc[b] / (5.0 * P.sims));
  std::fprintf(f, "\n");
  std::fclose(f);
  std::printf("cases %zu of %lld (%.4f%%): linear %d branched %d late %d\n", cases.size(), P.sims,
              100.0 * cases.size() / P.sims, cnt[0], cnt[1], cnt[2]);
  // per case table
  fs = P.out + "_cases.tsv";
  f = std::fopen(fs.c_str(), "w");
  std::fprintf(f, "sim\tdx\tpattern\tfirst_detect\tocc1\tocc2\tocc3\tocc4\tinc1\tinc2\tinc3\tinc4\tfit1\tfit2\tfit3\tfit4\n");
  for (auto &c : cases) {
    std::fprintf(f, "%lld\t%.1f\t%d\t%.1f", c.sim, c.dx, c.pattern, c.first_detect);
    for (int j = 0; j < 4; j++) std::fprintf(f, "\t%.1f", c.occ[j]);
    for (int j = 0; j < 4; j++) std::fprintf(f, "\t%.4f", c.inc[j]);
    for (int j = 0; j < 4; j++) std::fprintf(f, "\t%.4f", c.fit[j]);
    std::fprintf(f, "\n");
  }
  std::fclose(f);
  // PPV table: per life, dx age and VAF summaries
  fs = P.out + "_ppv.bin";
  f = std::fopen(fs.c_str(), "wb");
  long long ns = P.sims;
  std::fwrite(&ns, sizeof ns, 1, f);
  std::vector<float> dxf(P.sims);
  for (long long i = 0; i < P.sims; i++) dxf[i] = (float)dxage[i];
  std::fwrite(dxf.data(), sizeof(float), P.sims, f);
  std::fwrite(ppvrec.data(), sizeof(float), ppvrec.size(), f);
  std::fclose(f);
  if (P.dump_cases) {
    fs = P.out + "_cases.bin";
    f = std::fopen(fs.c_str(), "wb");
    int Ti = T;
    long long nc = (long long)cases.size();
    std::fwrite(&nc, sizeof nc, 1, f);
    std::fwrite(&Ti, sizeof Ti, 1, f);
    std::fwrite(&dt, sizeof dt, 1, f);
    for (auto &c : cases) {
      std::fwrite(&c.sim, sizeof c.sim, 1, f);
      std::fwrite(&c.dx, sizeof c.dx, 1, f);
      int ncl = (int)c.clones.size();
      std::fwrite(&ncl, sizeof ncl, 1, f);
      for (auto &cl : c.clones) {
        // id of a clone = its last mutation id (root: 0)
        int id = cl.nmut > 0 ? cl.muts[cl.nmut - 1] : 0;
        std::fwrite(&id, sizeof id, 1, f);
        std::fwrite(&cl.parent, sizeof cl.parent, 1, f);
        std::fwrite(&cl.nmut, sizeof cl.nmut, 1, f);
        std::fwrite(cl.muts, sizeof(int), 4, f);
        std::fwrite(cl.incs, sizeof(double), 4, f);
        std::fwrite(&cl.occ, sizeof cl.occ, 1, f);
        std::fwrite(&cl.fit, sizeof cl.fit, 1, f);
        char am = cl.aml ? 1 : 0;
        std::fwrite(&am, 1, 1, f);
        std::fwrite(cl.traj.data(), sizeof(float), T, f);
      }
    }
    std::fclose(f);
  }
  return 0;
}
