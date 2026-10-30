// Exact per-map reachability of 3D-mesh routing functions under link faults.
//
// For every fault map and every ordered pair (s, d), the probability that a
// packet injected at s reaches d is computed exactly:
//   xyz, zxy          dimension-order routing (deterministic)
//   o1t               O1TURN-3D, XYZ or ZYX chosen with probability 1/2
//   romm, val         two-phase XYZ through a uniform intermediate node in the
//                     minimal box (ROMM) or anywhere (Valiant)
//   nf_oracle/greedy  3D negative-first; oracle = some admissible healthy path
//                     exists; greedy = uniform choice among healthy admissible
//                     ports (local fault information only)
//   ham_oracle/greedy Hamiltonian label-monotone routing (may be non-minimal)
//   min_oracle/greedy minimal fully adaptive (FA family)
//   min_order         minimal, first healthy port in X, Y, Z order (ordered
//                     minimal adaptive routing)
//   connected         s and d connected in the faulty graph (upper bound for
//                     any fault-tolerant routing, including non-minimal)
// Planar links fail i.i.d. with probability p_h; vertical segments follow a
// Gaussian copula with exponential spatial covariance exp(-d/xi) and exact
// marginal p_v. Output: CSV, one row per (p_h, p_v, xi, policy) with the mean
// pair reachability, its 95% bootstrap interval over maps and the probability
// that at least one pair is unreachable.
#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

constexpr int KX = 8, KY = 8, KZ = 4, N = KX * KY * KZ;
enum { E = 0, W, NO, S, U, D, P };
using Row = std::bitset<N>;

int X(int r) { return r % KX; }
int Y(int r) { return (r / KX) % KY; }
int Z(int r) { return r / (KX * KY); }
int dist(int a, int b) { return std::abs(X(a) - X(b)) + std::abs(Y(a) - Y(b)) + std::abs(Z(a) - Z(b)); }
int nb(int r, int p) {
  const int x = X(r), y = Y(r), z = Z(r);
  switch (p) {
    case E: return x + 1 < KX ? r + 1 : -1;
    case W: return x > 0 ? r - 1 : -1;
    case NO: return y + 1 < KY ? r + KX : -1;
    case S: return y > 0 ? r - KX : -1;
    case U: return z + 1 < KZ ? r + KX * KY : -1;
    case D: return z > 0 ? r - KX * KY : -1;
  }
  return -1;
}

double inv_norm(double p) {
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                             6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                             -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                             3.754408661907416e+00};
  if (p <= 0) return -1e9;
  if (p < 0.02425) {
    const double q = std::sqrt(-2 * std::log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  if (p > 1 - 0.02425) return -inv_norm(1 - p);
  const double q = p - 0.5, r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
         (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}

std::vector<char> sample(std::mt19937_64& g, double ph, double pv, const std::vector<double>& L) {
  std::vector<char> ok(N * P, 1);
  std::uniform_real_distribution<double> u(0, 1);
  std::normal_distribution<double> nd(0, 1);
  for (int r = 0; r < N; ++r)
    for (int p : {E, NO}) {
      const int n = nb(r, p);
      if (n < 0) continue;
      if (u(g) < ph) {
        ok[r * P + p] = 0;
        ok[n * P + (p == E ? W : S)] = 0;
      }
    }
  const int C = KX * KY;
  const double thr = inv_norm(pv);
  for (int b = 0; b < KZ - 1; ++b) {
    std::vector<double> e(C), z(C, 0.0);
    for (auto& v : e) v = nd(g);
    for (int i = 0; i < C; ++i)
      for (int k = 0; k <= i; ++k) z[i] += L[i * C + k] * e[k];
    for (int i = 0; i < C; ++i)
      if (pv > 0 && z[i] < thr) {
        const int lo = i + b * C, hi = lo + C;
        ok[lo * P + U] = 0;
        ok[hi * P + D] = 0;
      }
  }
  return ok;
}

int dor_port(int r, int d, const int* order) {
  for (int k = 0; k < 3; ++k) {
    const int dim = order[k];
    const int delta = dim == 0 ? X(d) - X(r) : (dim == 1 ? Y(d) - Y(r) : Z(d) - Z(r));
    if (delta == 0) continue;
    if (dim == 0) return delta > 0 ? E : W;
    if (dim == 1) return delta > 0 ? NO : S;
    return delta > 0 ? U : D;
  }
  return -1;
}

void minimal_ports(int r, int d, int* ports, int& k) {
  k = 0;
  if (X(d) != X(r)) ports[k++] = X(d) > X(r) ? E : W;
  if (Y(d) != Y(r)) ports[k++] = Y(d) > Y(r) ? NO : S;
  if (Z(d) != Z(r)) ports[k++] = Z(d) > Z(r) ? U : D;
}

void nf_ports(int r, int d, int* ports, int& k) {
  k = 0;
  if (X(d) < X(r)) ports[k++] = W;
  if (Y(d) < Y(r)) ports[k++] = S;
  if (Z(d) < Z(r)) ports[k++] = D;
  if (k) return;
  if (X(d) > X(r)) ports[k++] = E;
  if (Y(d) > Y(r)) ports[k++] = NO;
  if (Z(d) > Z(r)) ports[k++] = U;
}

enum Policy {
  XYZ, ZXY, O1T, ROMM, VAL, NF_ORACLE, NF_GREEDY, HAM_ORACLE, HAM_GREEDY, MIN_ORACLE, MIN_GREEDY, MIN_ORDER,
  CONNECTED, NPOL
};
const char* kName[] = {"xyz", "zxy", "o1t", "romm", "val", "nf_oracle", "nf_greedy", "ham_oracle",
                       "ham_greedy", "min_oracle", "min_greedy", "min_order", "connected"};

struct Static {
  std::vector<std::vector<int>> by_dist;  // nodes sorted by distance to d
  std::vector<int> ham;                   // snake labels
  std::vector<std::vector<int>> by_label; // for each d: nodes in processing order for HAM
  std::vector<Row> box;                   // box[s*N+d]: nodes in the minimal box
};

Static build_static() {
  Static st;
  st.by_dist.assign(N, {});
  for (int d = 0; d < N; ++d) {
    std::vector<int> v(N);
    std::iota(v.begin(), v.end(), 0);
    std::stable_sort(v.begin(), v.end(), [&](int a, int b) { return dist(a, d) < dist(b, d); });
    st.by_dist[d] = v;
  }
  st.ham.assign(N, 0);
  const int plane = KX * KY;
  for (int r = 0; r < N; ++r) {
    int idx = (Y(r) % 2 == 0) ? Y(r) * KX + X(r) : Y(r) * KX + (KX - 1 - X(r));
    if (Z(r) % 2 == 1) idx = plane - 1 - idx;
    st.ham[r] = Z(r) * plane + idx;
  }
  st.box.assign(static_cast<size_t>(N) * N, Row());
  for (int s = 0; s < N; ++s)
    for (int d = 0; d < N; ++d) {
      Row& b = st.box[static_cast<size_t>(s) * N + d];
      for (int m = 0; m < N; ++m)
        if (X(m) >= std::min(X(s), X(d)) && X(m) <= std::max(X(s), X(d)) && Y(m) >= std::min(Y(s), Y(d)) &&
            Y(m) <= std::max(Y(s), Y(d)) && Z(m) >= std::min(Z(s), Z(d)) && Z(m) <= std::max(Z(s), Z(d)))
          b.set(m);
    }
  return st;
}

// Deterministic DOR reachability table ok[s] row over destinations, and the
// transposed table (rows over sources) for bitset intersections.
void dor_table(const std::vector<char>& ok, const Static& st, const int* order, std::vector<Row>& fwd,
               std::vector<Row>& rev) {
  fwd.assign(N, Row());
  rev.assign(N, Row());
  std::vector<char> good(N);
  for (int d = 0; d < N; ++d) {
    for (int u : st.by_dist[d]) {
      if (u == d) { good[u] = 1; continue; }
      const int p = dor_port(u, d, order);
      good[u] = ok[u * P + p] && good[nb(u, p)];
    }
    for (int u = 0; u < N; ++u)
      if (good[u]) { fwd[u].set(d); rev[d].set(u); }
  }
}

void evaluate(const std::vector<char>& ok, const Static& st, double* reach, int* broken) {
  static const int kXYZ[3] = {0, 1, 2}, kZXY[3] = {2, 0, 1}, kZYX[3] = {2, 1, 0};
  for (int pol = 0; pol < NPOL; ++pol) { reach[pol] = 0; broken[pol] = 0; }
  std::vector<Row> xyz_f, xyz_r, zxy_f, zxy_r, zyx_f, zyx_r;
  dor_table(ok, st, kXYZ, xyz_f, xyz_r);
  dor_table(ok, st, kZXY, zxy_f, zxy_r);
  dor_table(ok, st, kZYX, zyx_f, zyx_r);
  // Oblivious policies from the bit tables.
  for (int s = 0; s < N; ++s)
    for (int d = 0; d < N; ++d) {
      if (s == d) continue;
      const double vx = xyz_f[s][d], vz = zxy_f[s][d];
      const double vo = 0.5 * (xyz_f[s][d] + zyx_f[s][d]);
      const Row& box = st.box[static_cast<size_t>(s) * N + d];
      const Row through = xyz_f[s] & xyz_r[d];
      const double vr = static_cast<double>((through & box).count()) / box.count();
      const double vv = static_cast<double>(through.count()) / N;
      reach[XYZ] += vx; reach[ZXY] += vz; reach[O1T] += vo; reach[ROMM] += vr; reach[VAL] += vv;
      if (vx < 1) broken[XYZ] = 1;
      if (vz < 1) broken[ZXY] = 1;
      if (vo < 1) broken[O1T] = 1;
      if (vr < 1) broken[ROMM] = 1;
      if (vv < 1) broken[VAL] = 1;
    }
  // Topological connectivity (union-find on healthy links).
  std::vector<int> parent(N);
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&](int a) {
    while (parent[a] != a) a = parent[a] = parent[parent[a]];
    return a;
  };
  for (int r = 0; r < N; ++r)
    for (int p = 0; p < P; ++p)
      if (nb(r, p) >= 0 && ok[r * P + p]) parent[find(r)] = find(nb(r, p));
  std::vector<int> comp(N, 0);
  for (int r = 0; r < N; ++r) ++comp[find(r)];
  double con = 0;
  for (int r = 0; r < N; ++r) con += comp[find(r)] - 1;
  reach[CONNECTED] = con;
  broken[CONNECTED] = con < static_cast<double>(N) * (N - 1);
  // Adaptive and turn-model policies by dynamic programming.
  std::vector<double> prob(N);
  for (int d = 0; d < N; ++d) {
    for (int pol : {NF_ORACLE, NF_GREEDY, MIN_ORACLE, MIN_GREEDY, MIN_ORDER}) {
      for (int u : st.by_dist[d]) {
        if (u == d) { prob[u] = 1.0; continue; }
        int ports[3], k;
        if (pol == NF_ORACLE || pol == NF_GREEDY) nf_ports(u, d, ports, k);
        else minimal_ports(u, d, ports, k);
        double acc = 0.0, best = 0.0, first = -1.0;
        int healthy = 0;
        for (int i = 0; i < k; ++i) {
          if (!ok[u * P + ports[i]]) continue;
          const double v = prob[nb(u, ports[i])];
          ++healthy;
          acc += v;
          best = std::max(best, v);
          if (first < 0) first = v;
        }
        if (healthy == 0) prob[u] = 0.0;
        else if (pol == NF_ORACLE || pol == MIN_ORACLE) prob[u] = best > 0.5 ? 1.0 : 0.0;
        else if (pol == MIN_ORDER) prob[u] = first;
        else prob[u] = acc / healthy;
      }
      double s = 0.0;
      for (int u = 0; u < N; ++u)
        if (u != d) {
          s += prob[u];
          if (prob[u] < 1.0 - 1e-12) broken[pol] = 1;
        }
      reach[pol] += s;
    }
    // Hamiltonian routing: process nodes from the destination label outwards.
    for (int pol : {HAM_ORACLE, HAM_GREEDY}) {
      const int ld = st.ham[d];
      std::vector<int> order(N);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(),
                [&](int a, int b) { return std::abs(st.ham[a] - ld) < std::abs(st.ham[b] - ld); });
      for (int u : order) {
        if (u == d) { prob[u] = 1.0; continue; }
        const int lu = st.ham[u];
        std::vector<int> adm, red;
        for (int p = 0; p < P; ++p) {
          const int n = nb(u, p);
          if (n < 0 || !ok[u * P + p]) continue;
          const int ln = st.ham[n];
          const bool admissible = ld > lu ? (ln > lu && ln <= ld) : (ln < lu && ln >= ld);
          if (!admissible) continue;
          adm.push_back(n);
          if (dist(n, d) < dist(u, d)) red.push_back(n);
        }
        if (adm.empty()) { prob[u] = 0.0; continue; }
        if (pol == HAM_ORACLE) {
          double best = 0.0;
          for (int n : adm) best = std::max(best, prob[n]);
          prob[u] = best;
        } else if (!red.empty()) {
          double acc = 0.0;
          for (int n : red) acc += prob[n];
          prob[u] = acc / red.size();
        } else {
          int bestn = adm.front();
          for (int n : adm)
            if (std::abs(st.ham[n] - ld) < std::abs(st.ham[bestn] - ld)) bestn = n;
          prob[u] = prob[bestn];
        }
      }
      double s = 0.0;
      for (int u = 0; u < N; ++u)
        if (u != d) {
          s += prob[u];
          if (prob[u] < 1.0 - 1e-12) broken[pol] = 1;
        }
      reach[pol] += s;
    }
  }
  for (int pol = 0; pol < NPOL; ++pol) reach[pol] /= static_cast<double>(N) * (N - 1);
}

}  // namespace

int main(int argc, char** argv) {
  const int maps = argc > 1 ? std::atoi(argv[1]) : 500;
  std::printf("p_h,p_v,xi,policy,reach_mean,reach_lo,reach_hi,p_any_unreachable\n");
  struct Cfg { double ph, pv, xi; };
  std::vector<Cfg> cfgs;
  for (double ph : {0.0, 0.01, 0.02, 0.05, 0.10}) cfgs.push_back({ph, 0.0, 0.0});
  for (double pv : {0.05, 0.10, 0.20, 0.30})
    for (double xi : {0.0, 4.0}) cfgs.push_back({0.0, pv, xi});
  for (double xi : {0.0, 4.0}) cfgs.push_back({0.02, 0.10, xi});
  const Static st = build_static();
  const int C = KX * KY;
  unsigned long long seed = 20261024ULL;
  for (const Cfg& cf : cfgs) {
    std::vector<double> L(C * C, 0.0), S(C * C, 0.0);
    for (int i = 0; i < C; ++i)
      for (int j = 0; j < C; ++j) {
        const double dx = (i % KX) - (j % KX), dy = (i / KX) - (j / KX);
        const double dd = std::sqrt(dx * dx + dy * dy);
        S[i * C + j] = cf.xi > 0 ? std::exp(-dd / cf.xi) : (i == j ? 1.0 : 0.0);
      }
    for (int i = 0; i < C; ++i)
      for (int j = 0; j <= i; ++j) {
        double s = S[i * C + j];
        for (int k = 0; k < j; ++k) s -= L[i * C + k] * L[j * C + k];
        L[i * C + j] = (i == j) ? std::sqrt(std::max(s, 1e-12)) : s / L[j * C + j];
      }
    std::mt19937_64 g(seed++);
    std::vector<std::vector<double>> r(NPOL);
    std::vector<int> brk(NPOL, 0);
    for (int m = 0; m < maps; ++m) {
      const auto ok = sample(g, cf.ph, cf.pv, L);
      double reach[NPOL];
      int broken[NPOL];
      evaluate(ok, st, reach, broken);
      for (int p = 0; p < NPOL; ++p) { r[p].push_back(reach[p]); brk[p] += broken[p]; }
    }
    std::mt19937_64 bg(seed * 31ULL);
    for (int p = 0; p < NPOL; ++p) {
      double mean = 0;
      for (double v : r[p]) mean += v;
      mean /= maps;
      std::vector<double> boots;
      std::uniform_int_distribution<int> pick(0, maps - 1);
      for (int b = 0; b < 1000; ++b) {
        double s = 0;
        for (int i = 0; i < maps; ++i) s += r[p][pick(bg)];
        boots.push_back(s / maps);
      }
      std::sort(boots.begin(), boots.end());
      std::printf("%.3f,%.3f,%.1f,%s,%.6f,%.6f,%.6f,%.4f\n", cf.ph, cf.pv, cf.xi, kName[p], mean, boots[25],
                  boots[974], static_cast<double>(brk[p]) / maps);
    }
    std::fflush(stdout);
  }
  return 0;
}
