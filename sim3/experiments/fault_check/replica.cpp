// Replica of noc3d_sim::sample_faults (same RNG, same Cholesky) for the VPC stride-2 layout of an 8x8x4 mesh.
// For each (seed, pv, xi) it prints the realised fraction of failed elevator segments and the exact expected
// packet loss of Elevator-First (nearest elevator) under uniform random traffic, given that fault map.
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <vector>
static double inv(double p) {
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                             6.680131188771972e+01, -1.328068155288572e+01};
  static const double cc[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                              -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                             3.754408661907416e+00};
  const double plow = 0.02425;
  if (p < plow) { const double q = std::sqrt(-2 * std::log(p));
    return (((((cc[0] * q + cc[1]) * q + cc[2]) * q + cc[3]) * q + cc[4]) * q + cc[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  if (p > 1 - plow) return -inv(1 - p);
  const double q = p - 0.5, r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
         (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}
int main(int argc, char** argv) {
  const unsigned long long M = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 5;
  const int KX = 8, KY = 8, KZ = 4;
  std::vector<std::pair<int,int>> E;
  for (int x = 0; x < KX; ++x) for (int y = 0; y < KY; ++y) if (x % 2 == 1 && y % 2 == 1) E.push_back({x, y});
  const int ne = E.size();
  for (double pv : {0.1, 0.2}) for (double xi : {0.0, 4.0}) for (unsigned long long seed = 1; seed <= M; ++seed) {
    std::vector<double> L(ne * ne, 0.0), S(ne * ne, 0.0);
    for (int i = 0; i < ne; ++i) for (int j = 0; j < ne; ++j) {
      const double dx = E[i].first - E[j].first, dy = E[i].second - E[j].second, dd = std::sqrt(dx * dx + dy * dy);
      S[i * ne + j] = xi > 0 ? std::exp(-dd / xi) : (i == j ? 1.0 : 0.0); }
    for (int i = 0; i < ne; ++i) for (int j = 0; j <= i; ++j) { double s = S[i * ne + j];
      for (int k = 0; k < j; ++k) s -= L[i * ne + k] * L[j * ne + k];
      L[i * ne + j] = (i == j) ? std::sqrt(std::max(s, 1e-12)) : s / L[j * ne + j]; }
    std::mt19937_64 frng(seed * 7919ULL + 17ULL);
    std::normal_distribution<double> g(0.0, 1.0);
    const double thr = inv(pv);
    std::vector<char> ok(ne * (KZ - 1), 1);
    int failed = 0;
    for (int bnd = 0; bnd < KZ - 1; ++bnd) {
      std::vector<double> u(ne), zl(ne, 0.0);
      for (auto& v : u) v = g(frng);
      for (int i = 0; i < ne; ++i) for (int k = 0; k <= i; ++k) zl[i] += L[i * ne + k] * u[k];
      for (int e = 0; e < ne; ++e) { ok[e * (KZ - 1) + bnd] = zl[e] >= thr; failed += !ok[e * (KZ - 1) + bnd]; }
    }
    // nearest elevator (planar Manhattan distance, first in enumeration order on ties) -- as nearest_elevator()
    double lost = 0, tot = 0;
    for (int s = 0; s < KX * KY * KZ; ++s) {
      const int sx = s % KX, sy = (s / KX) % KY, sz = s / (KX * KY);
      int best = -1, bd = 1 << 30;
      for (int e = 0; e < ne; ++e) { int dd = std::abs(sx - E[e].first) + std::abs(sy - E[e].second);
        if (dd < bd) { bd = dd; best = e; } }
      for (int dz = 0; dz < KZ; ++dz) {
        const int m = (dz == sz) ? KX * KY - 1 : KX * KY;   // destinations on tier dz (excluding s)
        bool span = true;
        for (int b2 = std::min(sz, dz); b2 < std::max(sz, dz); ++b2) span = span && ok[best * (KZ - 1) + b2];
        tot += m; if (!span) lost += m;
      }
    }
    std::printf("%.1f %.1f %llu %.4f %.4f\n", pv, xi, seed, failed / double(ne * (KZ - 1)), lost / tot);
  }
}
