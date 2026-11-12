// refit.cpp -- maximum-likelihood re-fit of the clone model of Watson et al. (Nature 2026,
// Methods eqs. 1-2) to the public trajectories, with profile-likelihood intervals.
//
// Model (paper eqs. 1-2, code Code_for_inferring_acquisition_age_and_fitness_v16.py):
//   size_c(t) = (exp(s_c (t - t0_c)) - 1)/s_c for t > t0_c, else 0
//   freq_c(t) = size_c(t) / (1e5 + sum_c' size_c'(t))
//   VAF_v(t)  = sum_{c contains v} w(v,c) freq_c(t),  w = 1/2 for SNV/indel (1 for JAK2 V617F in a
//               clone with chr9p CN-LOH and for TET2 in a clone with chr4q CN-LOH), w = 1 for mCAs.
// Reported age: t_e = t0 + log(1+s)/s (deterministic size reaches 1).
// Likelihood: binomial on read counts; undetected timepoints censored, P(K <= 2 | depth, VAF);
// mCA cell fractions Gaussian with the phased-caller s.d.
// Constraints: s in [0.002, 15]/yr, t0 >= 0 (birth), t0(daughter) >= t0(parent).  No constraint
// that a daughter is fitter than its parent (the authors' code penalises that in some branches).
//
// Build: clang++ -O3 -march=native -fopenmp=libomp -std=c++17 refit.cpp -o bin/refit
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <omp.h>

struct Obs { int kind; double a, b; double c0 = 0, c1 = 0, c2 = 0; };  // kind 0 none, 1 binom (k=a, n=b), 2 censored (n=b), 3 gauss (f=a, sd=b)
struct Var { std::string name; char type; std::vector<Obs> obs; };
struct Person {
  std::string id; double dx;
  std::vector<double> t;
  std::vector<Var> vars;
  std::vector<std::vector<int>> clones;      // variant indices
  std::vector<std::string> cname;
  std::vector<int> parent;
  std::vector<std::vector<double>> w;        // w[c][v] weight of variant v in clone c
};

static std::vector<std::string> split(const std::string &s, char d) {
  std::vector<std::string> r; std::string x; std::stringstream ss(s);
  while (std::getline(ss, x, d)) r.push_back(x);
  return r;
}
static std::string trim(const std::string &s) {
  size_t a = s.find_first_not_of(" \t"), b = s.find_last_not_of(" \t\r\n");
  return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

static bool load(const std::string &path, Person &P) {
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  while (std::getline(f, line)) {
    std::stringstream ss(line); std::string tag; ss >> tag;
    if (tag == "PERSON") ss >> P.id >> P.dx;
    else if (tag == "TIMES") { double x; while (ss >> x) P.t.push_back(x); }
    else if (tag == "VAR") {
      Var v; std::string ty; ss >> v.name >> ty; v.type = ty[0];
      std::string tok;
      while (ss >> tok) {
        Obs o{0, 0, 0};
        if (tok == "NA") o.kind = 0;
        else if (tok[0] == 'U') { o.kind = 2; o.b = std::atof(tok.c_str() + 1); }
        else if (tok[0] == 'F') { auto p = split(tok.substr(1), ':'); o.kind = 3; o.a = std::atof(p[0].c_str()); o.b = std::atof(p[1].c_str()); }
        else { auto p = split(tok, '/'); o.kind = 1; o.a = std::atof(p[0].c_str()); o.b = std::atof(p[1].c_str()); }
        if (o.kind == 1) o.c0 = std::lgamma(o.b + 1) - std::lgamma(o.a + 1) - std::lgamma(o.b - o.a + 1);
        if (o.kind == 2) { o.c1 = std::log(o.b); o.c2 = std::log(o.b * (o.b - 1) / 2.0); }
        v.obs.push_back(o);
      }
      if (v.type == 'H') v.obs.assign(P.t.size(), Obs{0, 0, 0});
      P.vars.push_back(v);
    } else if (tag == "CLONE") {
      std::string rest; std::getline(ss, rest);
      auto parts = split(rest, '|');
      std::vector<int> cl; std::string nm;
      for (auto &q : parts) {
        std::string n = trim(q);
        for (size_t i = 0; i < P.vars.size(); i++) if (P.vars[i].name == n) cl.push_back((int)i);
        nm += (nm.empty() ? "" : ",") + n;
      }
      P.clones.push_back(cl); P.cname.push_back(nm);
    }
  }
  int C = (int)P.clones.size(), V = (int)P.vars.size();
  P.parent.assign(C, -1);
  for (int c = 0; c < C; c++) {
    for (int drop = 1; drop <= 2 && P.parent[c] < 0; drop++) {
      if ((int)P.clones[c].size() <= drop) break;
      std::vector<int> pre(P.clones[c].begin(), P.clones[c].end() - drop);
      for (int d = 0; d < C; d++) if (d != c && P.clones[d] == pre) { P.parent[c] = d; break; }
    }
  }
  P.w.assign(C, std::vector<double>(V, 0.0));
  for (int c = 0; c < C; c++) {
    bool has9p = false, has4q = false;
    for (int v : P.clones[c]) { if (P.vars[v].name == "chr9p_CNLOH") has9p = true; if (P.vars[v].name == "chr4q_CNLOH") has4q = true; }
    for (int v : P.clones[c]) {
      const std::string &n = P.vars[v].name;
      double wt = 0.5;
      if (P.vars[v].type == 'M') wt = 1.0;
      if (n == "JAK2_p.V617F" && has9p) wt = 1.0;
      if (n.rfind("TET2", 0) == 0 && has4q) wt = 1.0;
      P.w[c][v] = wt;
    }
  }
  return true;
}

static double lbinom(double k, double n, double p) {
  p = std::min(std::max(p, 1e-12), 1 - 1e-12);
  return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1) + k * std::log(p) + (n - k) * std::log1p(-p);
}
static double lcens(double n, double p) {  // log P(K <= 2)
  p = std::min(std::max(p, 1e-12), 1 - 1e-12);
  double s = 0;
  for (int k = 0; k <= 2; k++) s += std::exp(lbinom(k, n, p));
  return std::log(std::max(s, 1e-300));
}

// parameter vector x: for clone c, x[2c] = log s, x[2c+1] = t0
struct Model {
  const Person *P;
  int C, V, T;
  int fixc = -1, fixwhich = -1; double fixval = 0;  // profile: fix log s (0) or t_e (1) of clone fixc
  void preds(const std::vector<double> &x, std::vector<double> &vaf) const {
    vaf.assign((size_t)V * T, 0.0);
    std::vector<double> size(C);
    for (int j = 0; j < T; j++) {
      double t = P->t[j], tot = 1e5;
      for (int c = 0; c < C; c++) {
        double s = std::exp(x[2 * c]), t0 = x[2 * c + 1];
        double z = s * (t - t0);
        size[c] = z > 0 ? std::expm1(std::min(z, 700.0)) / s : 0.0;
        tot += size[c];
      }
      for (int c = 0; c < C; c++) {
        double fr = size[c] / tot;
        if (fr == 0) continue;
        for (int v : P->clones[c]) vaf[(size_t)v * T + j] += P->w[c][v] * fr;
      }
    }
  }
  void expand(std::vector<double> &x) const {  // apply profile constraint
    if (fixc < 0) return;
    if (fixwhich == 0) x[2 * fixc] = fixval;
    else { double s = std::exp(x[2 * fixc]); x[2 * fixc + 1] = fixval - std::log1p(s) / s; }
  }
  double nll(std::vector<double> x) const {
    expand(x);
    double pen = 0;
    for (int c = 0; c < C; c++) {
      double ls = x[2 * c], t0 = x[2 * c + 1];
      double lo = std::log(0.002), hi = std::log(15.0);
      if (ls < lo) pen += 1e4 * (lo - ls); if (ls > hi) pen += 1e4 * (ls - hi);
      if (t0 < 0) pen += 1e3 * (-t0);
      if (t0 > P->dx) pen += 1e3 * (t0 - P->dx);
      int p = P->parent[c];
      if (p >= 0) { double tp = x[2 * p + 1]; if (t0 < tp) pen += 1e3 * (tp - t0); }
    }
    std::vector<double> vaf; preds(x, vaf);
    double L = 0;
    for (int v = 0; v < V; v++) {
      const Var &vr = P->vars[v];
      for (int j = 0; j < T; j++) {
        const Obs &o = vr.obs[j]; double p = vaf[(size_t)v * T + j];
        if (o.kind == 1) {
          double pp = std::min(std::max(p, 1e-12), 1 - 1e-12);
          L -= o.c0 + o.a * std::log(pp) + (o.b - o.a) * std::log1p(-pp);
        } else if (o.kind == 2) {
          double pp = std::min(std::max(p, 1e-12), 1 - 1e-12);
          double l0 = o.b * std::log1p(-pp), lr = std::log(pp) - std::log1p(-pp);
          double t0v = l0, t1v = o.c1 + l0 + lr, t2v = o.c2 + l0 + 2 * lr;
          double mx = std::max(t0v, std::max(t1v, t2v));
          L -= mx + std::log(std::exp(t0v - mx) + std::exp(t1v - mx) + std::exp(t2v - mx));
        }
        else if (o.kind == 3) { double d = (o.a - std::min(p, 1.0)) / o.b; L += 0.5 * d * d; }
      }
    }
    return L + pen;
  }
};

// Nelder-Mead with adaptive coefficients (Gao & Han 2012)
static double nelder_mead(const Model &M, std::vector<double> &x, const std::vector<double> &step, int maxev) {
  int n = (int)x.size();
  double a = 1, g = 1 + 2.0 / n, rh = 0.75 - 1.0 / (2 * n), sg = 1 - 1.0 / n;
  std::vector<std::vector<double>> S(n + 1, x);
  std::vector<double> F(n + 1);
  for (int i = 0; i < n; i++) S[i + 1][i] += step[i];
  for (int i = 0; i <= n; i++) F[i] = M.nll(S[i]);
  int ev = n + 1;
  std::vector<int> idx(n + 1);
  while (ev < maxev) {
    for (int i = 0; i <= n; i++) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int i, int j) { return F[i] < F[j]; });
    if (std::fabs(F[idx[n]] - F[idx[0]]) < 1e-9 * (1 + std::fabs(F[idx[0]]))) {
      double sp = 0; for (int i = 0; i < n; i++) sp = std::max(sp, std::fabs(S[idx[n]][i] - S[idx[0]][i]));
      if (sp < 1e-7) break;
    }
    std::vector<double> xc(n, 0.0);
    for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) xc[i] += S[idx[k]][i] / n;
    auto lin = [&](double c) { std::vector<double> y(n); for (int i = 0; i < n; i++) y[i] = xc[i] + c * (xc[i] - S[idx[n]][i]); return y; };
    auto xr = lin(a); double fr = M.nll(xr); ev++;
    if (fr < F[idx[0]]) {
      auto xe = lin(a * g); double fe = M.nll(xe); ev++;
      if (fe < fr) { S[idx[n]] = xe; F[idx[n]] = fe; } else { S[idx[n]] = xr; F[idx[n]] = fr; }
    } else if (fr < F[idx[n - 1]]) { S[idx[n]] = xr; F[idx[n]] = fr; }
    else {
      bool outside = fr < F[idx[n]];
      auto xk = lin(outside ? a * rh : -rh); double fk = M.nll(xk); ev++;
      if (fk < (outside ? fr : F[idx[n]])) { S[idx[n]] = xk; F[idx[n]] = fk; }
      else {
        for (int k = 1; k <= n; k++) {
          for (int i = 0; i < n; i++) S[idx[k]][i] = S[idx[0]][i] + sg * (S[idx[k]][i] - S[idx[0]][i]);
          F[idx[k]] = M.nll(S[idx[k]]); ev++;
        }
      }
    }
  }
  int b = 0; for (int i = 1; i <= n; i++) if (F[i] < F[b]) b = i;
  x = S[b];
  return F[b];
}

// warm-started local refinement used along profiles
static double polish(const Model &M, std::vector<double> &x, int maxev) {
  int n = (int)x.size();
  std::vector<double> step(n);
  for (int i = 0; i < n; i++) step[i] = (i % 2 == 0) ? 0.15 : 1.5;
  double f = nelder_mead(M, x, step, maxev);
  for (int i = 0; i < n; i++) step[i] *= 0.3;
  double f2 = nelder_mead(M, x, step, maxev);
  return std::min(f, f2);
}

static double optimise(const Model &M, std::vector<double> &x, std::mt19937_64 &g, int restarts, int maxev) {
  int n = (int)x.size();
  std::vector<double> step(n);
  for (int i = 0; i < n; i++) step[i] = (i % 2 == 0) ? 0.5 : 5.0;
  std::vector<double> best = x; double fb = nelder_mead(M, best, step, maxev);
  for (int k = 0; k < 2; k++) {  // re-start at the optimum
    std::vector<double> y = best; double fy = nelder_mead(M, y, step, maxev);
    if (fy < fb) { fb = fy; best = y; }
  }
  std::normal_distribution<double> N01(0, 1);
  for (int r = 0; r < restarts; r++) {
    std::vector<double> y = best;
    for (int i = 0; i < n; i++) y[i] += N01(g) * ((i % 2 == 0) ? 0.7 : 8.0);
    for (int i = 1; i < n; i += 2) if (y[i] < 0) y[i] = -y[i];
    double fy = nelder_mead(M, y, step, maxev);
    for (int k = 0; k < 2; k++) { double f2 = nelder_mead(M, y, step, maxev); if (f2 >= fy - 1e-9) { fy = std::min(fy, f2); break; } fy = f2; }
    if (fy < fb) { fb = fy; best = y; }
  }
  x = best;
  return fb;
}

int main(int argc, char **argv) {
  if (argc < 3 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
    std::printf("refit: ML re-fit of the Watson et al. clone model with profile-likelihood intervals.\n"
                "usage: refit <out.tsv> <person1.txt> [person2.txt ...] [--threads T] [--restarts R]\n"
                "writes one row per clone: MLE, profile 95%% interval (binomial and quasi-binomial\n"
                "with the pooled dispersion), and flags for intervals that reach the search bounds.\n");
    return 0;
  }
  std::string out = argv[1], pubfile;
  std::vector<std::string> files; int threads = 4, restarts = 12; bool noprof = false;
  for (int i = 2; i < argc; i++) {
    if (!std::strcmp(argv[i], "--threads")) threads = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--restarts")) restarts = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--published")) pubfile = argv[++i];
    else if (!std::strcmp(argv[i], "--no-profiles")) noprof = true;
    else files.push_back(argv[i]);
  }
  // published estimates: lines "person<TAB>clone(comma list)<TAB>s<TAB>t0"
  std::vector<std::vector<std::string>> pub;
  if (!pubfile.empty()) {
    std::ifstream pf(pubfile); std::string ln;
    while (std::getline(pf, ln)) { auto q = split(ln, '\t'); if (q.size() >= 4) pub.push_back(q); }
  }
  int NP = (int)files.size();
  std::vector<Person> PS(NP);
  for (int i = 0; i < NP; i++) if (!load(files[i], PS[i])) { std::fprintf(stderr, "cannot read %s\n", files[i].c_str()); return 1; }
  // Stage 1: MLE per person; Pearson statistic for dispersion
  std::vector<std::vector<double>> X(NP);
  std::vector<double> F(NP), chi2(NP), nobs(NP), FPUB(NP, NAN);
  omp_set_num_threads(threads);
  time_t t0 = time(nullptr);
  int done = 0;
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < NP; i++) {
    const Person &P = PS[i];
    Model M{&P, (int)P.clones.size(), (int)P.vars.size(), (int)P.t.size()};
    std::mt19937_64 g(12345 + i);
    // initial guess: s from log-slope of the clone's defining variant, t0 back-extrapolated
    std::vector<double> x(2 * M.C);
    for (int c = 0; c < M.C; c++) {
      int v = P.clones[c].back();
      double s = 0.3, tfirst = -1, ffirst = 0;
      std::vector<std::pair<double, double>> pts;
      for (int j = 0; j < M.T; j++) {
        const Obs &o = P.vars[v].obs[j];
        double f = o.kind == 1 ? o.a / o.b : (o.kind == 3 ? o.a * 0.5 : -1);
        if (f > 0 && f < 0.4) pts.push_back({P.t[j], std::log(f)});
        if (f > 0 && tfirst < 0) { tfirst = P.t[j]; ffirst = f; }
      }
      if (pts.size() >= 2) {
        double mt = 0, my = 0; for (auto &q : pts) { mt += q.first; my += q.second; } mt /= pts.size(); my /= pts.size();
        double sxy = 0, sxx = 0; for (auto &q : pts) { sxy += (q.first - mt) * (q.second - my); sxx += (q.first - mt) * (q.first - mt); }
        if (sxx > 0) s = std::min(std::max(sxy / sxx, 0.02), 5.0);
      }
      if (tfirst < 0) { tfirst = P.t.back(); ffirst = 0.001; }
      double t0g = tfirst - std::log(std::max(2 * ffirst, 1e-4) * 1e5 * s) / s;
      x[2 * c] = std::log(s); x[2 * c + 1] = std::min(std::max(t0g, 0.5), P.dx - 0.5);
    }
    for (int c = 0; c < M.C; c++) { int p = P.parent[c]; if (p >= 0 && x[2 * c + 1] < x[2 * p + 1]) x[2 * c + 1] = x[2 * p + 1] + 0.5; }
    // published parameters, if all clones present
    std::vector<double> xp(2 * M.C, NAN); int found = 0;
    for (auto &q : pub) {
      if (q[0] != P.id) continue;
      std::string nm; for (auto &z : split(q[1], ',')) nm += (nm.empty() ? "" : ",") + trim(z);
      for (int c = 0; c < M.C; c++) if (P.cname[c] == nm) { xp[2 * c] = std::log(std::max(std::atof(q[2].c_str()), 1e-4)); xp[2 * c + 1] = std::atof(q[3].c_str()); found++; }
    }
    double fpub = NAN;
    if (found == M.C) {
      // public ages are dx rounded down minus exact months, so they lag the true ages by an
      // unknown delta in [0,1): evaluate the published parameters at the best common shift
      fpub = INFINITY; double bestd = 0;
      for (int k = 0; k <= 50; k++) {
        double d = k * 0.02; std::vector<double> y = xp;
        for (int c = 0; c < M.C; c++) y[2 * c + 1] -= d;
        double fy = M.nll(y);
        if (fy < fpub) { fpub = fy; bestd = d; }
      }
      for (int c = 0; c < M.C; c++) xp[2 * c + 1] -= bestd;
    }
    F[i] = optimise(M, x, g, restarts, 4000 * M.C);
    if (found == M.C) {
      std::vector<double> y = xp; double fy = optimise(M, y, g, restarts / 2, 4000 * M.C);
      if (fy < F[i]) { F[i] = fy; x = y; }
    }
    X[i] = x; FPUB[i] = fpub;
    // Pearson chi-square over detected SNV counts
    std::vector<double> vaf; M.preds(x, vaf);
    double c2 = 0, no = 0;
    for (int v = 0; v < M.V; v++) for (int j = 0; j < M.T; j++) {
      const Obs &o = P.vars[v].obs[j];
      if (o.kind != 1) continue;
      double p = std::min(std::max(vaf[(size_t)v * M.T + j], 1e-9), 1 - 1e-9);
      double e = o.b * p;
      c2 += (o.a - e) * (o.a - e) / (e * (1 - p)); no += 1;
    }
    chi2[i] = c2; nobs[i] = no;
#pragma omp critical
    {
      done++;
      time_t now = time(nullptr); char buf[32]; std::strftime(buf, sizeof buf, "%H:%M:%SZ", std::gmtime(&now));
      std::printf("[%s] MLE %d/%d %s clones %d nll %.2f (at published %.2f) chi2/obs %.1f/%d elapsed %lds\n", buf, done, NP, P.id.c_str(), M.C, F[i], FPUB[i], c2, (int)no, (long)(now - t0));
      std::fflush(stdout);
    }
  }
  double C2 = 0, NO = 0, NPAR = 0;
  for (int i = 0; i < NP; i++) { C2 += chi2[i]; NO += nobs[i]; NPAR += 2 * PS[i].clones.size(); }
  double phi = C2 / std::max(1.0, NO - NPAR);
  std::printf("pooled Pearson dispersion phi = %.2f (chi2 %.1f, obs %.0f, params %.0f)\n", phi, C2, NO, NPAR);
  std::fflush(stdout);
  // Stage 2: profiles
  struct Row { std::string person, clone; int nmut; double s, te, t0, slo, shi, slo_q, shi_q, telo, tehi, telo_q, tehi_q; int sflag, teflag; };
  std::vector<std::vector<Row>> rows(NP);
  done = 0;
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < NP; i++) {
    const Person &P = PS[i];
    Model M{&P, (int)P.clones.size(), (int)P.vars.size(), (int)P.t.size()};
    std::mt19937_64 g(777 + i);
    for (int c = 0; c < M.C; c++) {
      Row R; R.person = P.id; R.clone = P.cname[c]; R.nmut = (int)P.clones[c].size();
      double s = std::exp(X[i][2 * c]); R.s = s; R.t0 = X[i][2 * c + 1]; R.te = R.t0 + std::log1p(s) / s;
      R.slo = R.shi = R.slo_q = R.shi_q = R.telo = R.tehi = R.telo_q = R.tehi_q = NAN; R.sflag = R.teflag = 0;
      for (int which = 0; which < (noprof ? 0 : 2); which++) {
        // walk outward from the MLE with warm starts until the quasi threshold + 2 is passed
        double thrq = 1.92 * std::max(phi, 1.0);
        double v0 = which == 0 ? X[i][2 * c] : R.te;
        double stepv = which == 0 ? 0.08 : 0.5;
        double vmin = which == 0 ? std::log(0.002) : 0.0, vmax = which == 0 ? std::log(15.0) : P.dx;
        std::vector<std::pair<double, double>> pts; pts.push_back({v0, F[i]});
        int hitlo = 0, hithi = 0;
        for (int dir = -1; dir <= 1; dir += 2) {
          std::vector<double> x = X[i]; double v = v0;
          for (int k = 0; k < 4000; k++) {
            v += dir * stepv;
            if (v < vmin || v > vmax) { if (dir < 0) hitlo = 1; else hithi = 1; break; }
            Model Mp = M; Mp.fixc = c; Mp.fixwhich = which; Mp.fixval = v;
            double fv = polish(Mp, x, 1500 * M.C);
            pts.push_back({v, fv});
            if (fv - F[i] > thrq + 0.5) break;
          }
        }
        std::sort(pts.begin(), pts.end());
        double fmin = F[i]; for (auto &q : pts) fmin = std::min(fmin, q.second);
        auto ci = [&](double thr, double &lo, double &hi, int &flag) {
          int a = -1, b = -1;
          for (size_t k = 0; k < pts.size(); k++) if (pts[k].second - fmin <= thr) { if (a < 0) a = (int)k; b = (int)k; }
          flag = 0;
          if (a <= 0) { lo = pts[0].first; flag |= hitlo ? 1 : 0; }
          else { double x1 = pts[a - 1].first, y1 = pts[a - 1].second - fmin, x2 = pts[a].first, y2 = pts[a].second - fmin; lo = x1 + (thr - y1) * (x2 - x1) / (y2 - y1); }
          if (b >= (int)pts.size() - 1) { hi = pts.back().first; flag |= hithi ? 2 : 0; }
          else { double x1 = pts[b].first, y1 = pts[b].second - fmin, x2 = pts[b + 1].first, y2 = pts[b + 1].second - fmin; hi = x1 + (thr - y1) * (x2 - x1) / (y2 - y1); }
          if (which == 0) { lo = std::exp(lo); hi = std::exp(hi); }
        };
        int f1 = 0, f2 = 0;
        if (which == 0) { ci(1.92, R.slo, R.shi, f1); ci(thrq, R.slo_q, R.shi_q, f2); R.sflag = f2; }
        else { ci(1.92, R.telo, R.tehi, f1); ci(thrq, R.telo_q, R.tehi_q, f2); R.teflag = f2; }
      }
      rows[i].push_back(R);
    }
#pragma omp critical
    {
      done++;
      time_t now = time(nullptr); char buf[32]; std::strftime(buf, sizeof buf, "%H:%M:%SZ", std::gmtime(&now));
      std::printf("[%s] profiles %d/%d %s (%d clones) elapsed %lds\n", buf, done, NP, P.id.c_str(), M.C, (long)(now - t0));
      std::fflush(stdout);
    }
  }
  FILE *f = std::fopen(out.c_str(), "w");
  std::fprintf(f, "# pooled Pearson dispersion phi = %.3f\n", phi);
  std::fprintf(f, "person\tclone\tnmut\ts\tt0\tte\ts_lo\ts_hi\ts_lo_q\ts_hi_q\tte_lo\tte_hi\tte_lo_q\tte_hi_q\ts_bound\tte_bound\tnll\tnll_published\n");
  for (int i = 0; i < NP; i++) for (auto &R : rows[i])
    std::fprintf(f, "%s\t%s\t%d\t%.4f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.1f\t%.1f\t%.1f\t%.1f\t%d\t%d\t%.3f\t%.3f\n", R.person.c_str(), R.clone.c_str(), R.nmut,
                 R.s, R.t0, R.te, R.slo, R.shi, R.slo_q, R.shi_q, R.telo, R.tehi, R.telo_q, R.tehi_q, R.sflag, R.teflag, F[i], FPUB[i]);
  std::fclose(f);
  return 0;
}
