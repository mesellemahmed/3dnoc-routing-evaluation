// Cycle-level wormhole/virtual-channel simulator for 3D mesh networks-on-chip.
//
// Version 3. One binary provides (i) the simulator and (ii) an exact
// channel-dependency-graph (CDG) checker that evaluates the routing relation
// used by the simulator (option --check 1).
//
// Topologies
//   fc  : fully connected 3D mesh, every router has vertical links.
//   vpc : vertically partially connected mesh; elevator columns at
//         (x % 2 == 1, y % 2 == 1) [stride2] or (x % 4 == 1, y % 4 == 1)
//         [stride4]; vertical segments may fail (Gaussian copula, Eq. in text).
//
// Router model: input-queued, V virtual channels (VCs) per port, credit-based
// flow control, atomic VC reallocation (an output VC is granted only when the
// downstream buffer is empty and the previous packet's tail has left),
// rotating-priority VC allocation, separable input-first round-robin switch
// allocation, per-hop delay H (--hop-delay) and credit delay (--credit-delay).
//
// Deadlock handling: none, only detection. A run is reported deadlocked when
// no flit moves for --watchdog cycles while flits are in flight, or when the
// network fails to drain after injection stops (drain test; it also detects
// partial deadlocks that a global watchdog misses).
//
// Routing classes (--routing); "escape" = Duato escape network on VC0 (XYZ)
// with adaptive VCs 1..V-1 and no return from the escape network.
//   xyz, zxy        dimension-order routing
//   nf              3D negative-first turn model
//   dyad            NF, deterministic until the preferred port is congested
//   ham             Hamiltonian-path (snake labelling) routing, minimal first
//   o1t             O1TURN-3D: XYZ on VC class 0 or ZYX on class 1, chosen at
//                   injection
//   romm, val       two-phase oblivious: random intermediate in the minimal
//                   box (ROMM) or anywhere (Valiant), XYZ in each phase
//   fa              minimal fully adaptive, max-credit selection, escape
//   rca             minimal adaptive, regional congestion awareness (RCA-1D)
//   qr              minimal adaptive, Q-routing selection, escape
//   ta              minimal adaptive, thermal-aware selection, escape
//   ef, efr, aesg,  elevator-based VPC routing with direction virtual networks
//   adele           (Elevator-First): nearest elevator, nearest healthy
//                   elevator, global-load elevator selection, AdEle selection
//   fl              First-Last (positive-first / positive-last virtual
//                   networks), strict VC separation between VN0 and VN2
//
// Outputs (prefix --out):
//   .trace  16-byte records {int32 t, int32 p, int16 r, int16 next, uint8 xi,
//           uint8 c, uint16 port}; xi in {0 inject, 1 forward, 2 receive,
//           3 drop}; c = 1 for vertical links (ERA schema).
//   .buf    float32[W][N][7] mean input-buffer occupancy ratio per window.
//   .pkt    int32[P][11] {src, dst, t_create, t_head_out, t_eject, hops,
//           vhops, dropped, elevator, len, attacker}.
//   .act    int64[N][4] per-router counts in the measurement window
//           {crossbar flits, planar out flits, vertical out flits,
//           buffer writes}.
//   .topo   "kx ky kz n_windows window" then one line per directed link.
//   .stats  JSON summary.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <random>
#include <string>
#include <vector>

namespace {

enum Port { kLocal = 0, kEast, kWest, kNorth, kSouth, kUp, kDown, kPorts };
constexpr int kOpposite[kPorts] = {kLocal, kWest, kEast, kSouth, kNorth, kDown, kUp};
enum Event : uint8_t { kInject = 0, kForward = 1, kReceive = 2, kDrop = 3 };

enum Alg {
  kXYZ, kZXY, kNF, kDyAD, kHam, kO1T, kROMM, kVal,
  kFA, kRCA, kQR, kTA,
  kEF, kEFR, kAESG, kAdEle, kFL, kNumAlg
};

struct AlgSpec {
  const char* name;
  bool vpc;
  bool escape;
  int min_vcs;
};

constexpr AlgSpec kSpec[kNumAlg] = {
    {"xyz", false, false, 1},   {"zxy", false, false, 1},  {"nf", false, false, 1},
    {"dyad", false, false, 1},  {"ham", false, false, 1},  {"o1t", false, false, 2},
    {"romm", false, false, 2},  {"val", false, false, 2},  {"fa", false, true, 2},
    {"rca", false, true, 2},    {"qr", false, true, 2},    {"ta", false, true, 2},
    {"ef", true, false, 2},     {"efr", true, false, 2},   {"aesg", true, false, 2},
    {"adele", true, false, 2},  {"fl", true, false, 2},
};

struct Config {
  int kx = 8, ky = 8, kz = 4;
  std::string topo = "fc";
  std::string routing = "xyz";
  std::string traffic = "random";
  double pir = 0.006;
  uint64_t seed = 1;
  int cycles = 18000;
  int warmup = 3000;
  int window = 1000;
  int depth = 8;
  int vcs = 2;
  int pkt_flits = 8;
  int hop_delay = 1;
  int credit_delay = 0;
  double pv = 0.0;
  double xi = 0.0;
  double hotspot_frac = 0.2;
  std::string elevators = "stride2";
  std::string out = "run";
  bool trace = true;
  int watchdog = 2000;
  bool drain = true;
  int drain_max = 200000;
  bool check = false;
  bool dump_stuck = false;
  bool atomic = true;
  bool inject_escape = false;  // experiment: BookSim-like injection, half of the packets stay on VC0 from the source
  bool rand_sel = false;  // experiment: equal-priority (random) choice among adaptive candidates, as in BookSim
  bool alloc_cache = true;
  // Class parameters, fixed a priori.
  double aes_beta = 0.5;
  double q_eta = 0.5;
  double rca_weight = 0.5;
  int ta_period = 1000;
  int adele_k = 4;
  double adele_a = 0.2;
  double adele_xi = 0.05;
  double adele_thr = 0.1;
  double dyad_frac = 0.25;
  // Traffic extensions.
  double burst_factor = 8.0;
  double burst_on_mean = 40.0;
  double pareto_shape = 1.4;
  int attackers = 0;
  double attack_pir = 0.0;
  // Compact thermal model used by the thermal-aware class (SI units).
  double g_lat = 6.0e-3, g_vert = 0.14, g_sink = 0.0625;
  double p_core = 0.25, p_leak = 0.01, e_flit = 30e-12, freq = 1.0e9;
};

struct Flit {
  int32_t pkt;
  int32_t ready;
  bool head;
  bool tail;
};

struct Packet {
  int32_t src = 0, dst = 0, t_create = 0;
  int32_t t_head_out = -1, t_eject = -1, hops = 0, vhops = 0;
  int32_t elevator = -1, inter = -1;
  int16_t len = 8;
  int8_t dir = 0;
  uint8_t order = 0;
  bool dropped = false;
  bool attacker = false;
};

struct TraceRec {
  int32_t t;
  int32_t p;
  int16_t r;
  int16_t next;
  uint8_t xi;
  uint8_t c;
  uint16_t port;
};
static_assert(sizeof(TraceRec) == 16, "trace record must be 16 bytes");

class Ring {
 public:
  void init(int cap) {
    buf_.assign(cap, Flit{});
    head_ = size_ = 0;
  }
  bool empty() const { return size_ == 0; }
  int size() const { return size_; }
  const Flit& front() const { return buf_[head_]; }
  void push(const Flit& f) {
    buf_[(head_ + size_) % static_cast<int>(buf_.size())] = f;
    ++size_;
  }
  void pop() {
    head_ = (head_ + 1) % static_cast<int>(buf_.size());
    --size_;
  }

 private:
  std::vector<Flit> buf_;
  int head_ = 0, size_ = 0;
};

struct InputVC {
  Ring q;
  int out_port = -1;
  int out_vc = -1;
  int wait_from = -1;
  // Allocation-retry cache: a head whose VC allocation failed is retried only
  // after one of its candidate output ports may have freed a VC.
  bool failed = false;
  uint8_t n_cand_ports = 0;
  uint8_t cand_ports[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  uint64_t fail_sig = 0;
};

struct OutputVC {
  int owner_port = -1;
  int owner_vc = -1;
  int credits = 0;
};

struct Cand {
  int port;
  int vc_mask;
};

struct CreditMsg {
  int r, port, vc;
};

class Simulator {
 public:
  explicit Simulator(const Config& cfg)
      : c_(cfg), rng_(cfg.seed * 0x9E3779B97F4A7C15ULL + 1ULL), trng_(cfg.seed) {
    alg_ = kNumAlg;
    for (int a = 0; a < kNumAlg; ++a)
      if (c_.routing == kSpec[a].name) alg_ = static_cast<Alg>(a);
    if (alg_ == kNumAlg) die("unknown routing class");
    vpc_ = c_.topo == "vpc";
    if (!vpc_ && c_.topo != "fc") die("unknown topology");
    if (kSpec[alg_].vpc != vpc_) die("routing class incompatible with topology");
    if (c_.vcs < kSpec[alg_].min_vcs) die("routing class needs more virtual channels");
    n_ = c_.kx * c_.ky * c_.kz;
    xs_.resize(n_);
    ys_.resize(n_);
    zs_.resize(n_);
    for (int r = 0; r < n_; ++r) {
      xs_[r] = r % c_.kx;
      ys_[r] = (r / c_.kx) % c_.ky;
      zs_[r] = r / (c_.kx * c_.ky);
    }
    V_ = c_.vcs;
    all_ = (1 << V_) - 1;
    c0_ = (1 << std::max(1, V_ / 2)) - 1;
    c1_ = all_ & ~c0_;
    if (c1_ == 0) c1_ = all_;
    esc_ = 1;
    adp_ = all_ & ~1;
    if (adp_ == 0) adp_ = all_;
    build_topology();
    build_hamiltonian();
    in_.assign(static_cast<size_t>(n_) * kPorts * V_, InputVC{});
    outvc_.assign(static_cast<size_t>(n_) * kPorts * V_, OutputVC{});
    for (auto& v : in_) v.q.init(c_.depth);
    for (int r = 0; r < n_; ++r)
      for (int p = 0; p < kPorts; ++p)
        for (int v = 0; v < V_; ++v) ov(r, p, v).credits = (p == kLocal) ? (1 << 30) : c_.depth;
    srcq_.assign(n_, {});
    inj_flit_.assign(n_, 0);
    inj_vc_.assign(n_, -1);
    rr_in_.assign(static_cast<size_t>(n_) * kPorts, 0);
    rr_out_.assign(static_cast<size_t>(n_) * kPorts, 0);
    rot_.assign(n_, 0);
    n_windows_ = (c_.cycles - c_.warmup + c_.window - 1) / c_.window;
    buf_acc_.assign(static_cast<size_t>(n_windows_) * n_ * kPorts, 0.0);
    elev_load_.assign(elev_xy_.size(), 0);
    act_.assign(static_cast<size_t>(n_) * 4, 0);
    act_period_.assign(n_, 0);
    credit_pipe_.assign(c_.credit_delay + 1, {});
    pver_.assign(static_cast<size_t>(n_) * kPorts, 0);
    occ_.assign(n_, 0);
    // Classes whose candidate set is a pure function of the head's position
    // and header (no credit thresholds, tables, timers or random choices).
    cacheable_ = c_.alloc_cache && alg_ != kDyAD;
    if (alg_ == kRCA) {
      agg_.assign(static_cast<size_t>(n_) * kPorts, 0.0);
      agg_next_ = agg_;
    }
    if (alg_ == kQR) init_q();
    if (alg_ == kTA) init_thermal();
    if (alg_ == kAdEle) init_adele();
    init_traffic();
  }

  void run() {
    if (c_.trace) trace_ = std::fopen((c_.out + ".trace").c_str(), "wb");
    int t = 0;
    for (; t < c_.cycles; ++t) {
      step(t, false);
      if (in_flight_flits_ > 0 && t - last_move_ > c_.watchdog) {
        deadlock_cycle_ = t;
        break;
      }
    }
    if (trace_) std::fclose(trace_);
    trace_ = nullptr;
    cycles_run_ = t;
    if (deadlock_cycle_ < 0 && c_.drain) {
      const int start = std::max(t, c_.cycles);
      int u = start;
      last_move_ = std::max(last_move_, start);
      while ((in_flight_flits_ > 0 || partial_injections() > 0) && u - last_move_ <= c_.watchdog &&
             u - start < c_.drain_max) {
        step(u, true);
        ++u;
      }
      drain_cycles_ = u - start;
      stuck_flits_ = in_flight_flits_ + partial_injections();
      drained_ = stuck_flits_ == 0;
    } else {
      stuck_flits_ = in_flight_flits_;
      drained_ = in_flight_flits_ == 0;
    }
    write_outputs();
    if (c_.dump_stuck && !drained_) dump_stuck();
  }

  // Forensics of a run that failed to drain: one JSON line per non-empty
  // input VC with its queue content, its allocation and, for an unallocated
  // head, the candidate output VCs with their owners and credits.
  void dump_stuck() {
    std::FILE* f = std::fopen((c_.out + ".stuck").c_str(), "w");
    for (int r = 0; r < n_; ++r)
      for (int ip = 0; ip < kPorts; ++ip)
        for (int v = 0; v < V_; ++v) {
          InputVC& in = iv(r, ip, v);
          if (in.q.empty()) continue;
          std::fprintf(f, "{\"r\": %d, \"ip\": %d, \"v\": %d, \"out_port\": %d, \"out_vc\": %d, \"pkts\": [",
                       r, ip, v, in.out_port, in.out_vc);
          int last = -1;
          bool first = true;
          Ring copy = in.q;
          while (!copy.empty()) {
            const Flit fl = copy.front();
            copy.pop();
            if (fl.pkt != last) {
              std::fprintf(f, "%s[%d, %d, %d]", first ? "" : ", ", fl.pkt, fl.head ? 1 : 0, pkts_[fl.pkt].dst);
              first = false;
              last = fl.pkt;
            }
          }
          std::fprintf(f, "], \"cands\": [");
          const Flit& fr = in.q.front();
          if (in.out_port < 0 && fr.head) {
            candidates(r, ip, v, pkts_[fr.pkt], cands_);
            bool fc = true;
            for (const Cand& cd : cands_)
              for (int w = 0; w < V_; ++w) {
                if (!((cd.vc_mask >> w) & 1)) continue;
                const OutputVC& o = ov(r, cd.port, w);
                std::fprintf(f, "%s[%d, %d, %d, %d, %d]", fc ? "" : ", ", cd.port, w, o.owner_port, o.owner_vc,
                             cd.port == kLocal ? -1 : o.credits);
                fc = false;
              }
          }
          std::fprintf(f, "], \"credits\": %d}\n",
                       (in.out_port > 0) ? ov(r, in.out_port, in.out_vc).credits : -1);
        }
    std::fclose(f);
  }

  // ----- channel-dependency checker ----------------------------------------
  void check_cdg() {
    const int C = n_ * kPorts * V_;
    std::vector<uint64_t> bits((static_cast<size_t>(C) * C + 63) / 64, 0);
    std::vector<int> stamp(C, -1);
    std::vector<int> queue;
    std::vector<Cand> cands;
    long long states = 0, dead_ends = 0, missing_links = 0, nonminimal = 0, no_escape = 0;
    long long esc_to_adaptive = 0;
    int ctx_id = 0;
    auto chan = [&](int r, int p, int v) { return (r * kPorts + p) * V_ + v; };
    auto add_edge = [&](int a, int b) {
      const size_t k = static_cast<size_t>(a) * C + b;
      bits[k >> 6] |= 1ULL << (k & 63);
    };
    auto explore = [&](const Packet& pk, const std::vector<int>& sources) {
      ++ctx_id;
      queue.clear();
      for (int s : sources) {
        relation(s, kLocal, 0, pk, cands);
        for (const Cand& cd : cands) {
          if (cd.port == kLocal) continue;
          if (neighbor(s, cd.port) < 0) { ++missing_links; continue; }
          for (int v = 0; v < V_; ++v)
            if ((cd.vc_mask >> v) & 1) {
              const int ch = chan(s, cd.port, v);
              if (stamp[ch] != ctx_id) { stamp[ch] = ctx_id; queue.push_back(ch); }
            }
        }
      }
      for (size_t qi = 0; qi < queue.size(); ++qi) {
        const int ch = queue[qi];
        ++states;
        const int r = ch / (kPorts * V_), p = (ch / V_) % kPorts, v = ch % V_;
        const int nr = neighbor(r, p), ip = kOpposite[p];
        relation(nr, ip, v, pk, cands);
        if (cands.empty()) { ++dead_ends; continue; }
        const bool adaptive_state = kSpec[alg_].escape && v != 0;
        bool has_escape = false;
        const int dist = manhattan(nr, pk.dst);
        for (const Cand& cd : cands) {
          if (cd.port == kLocal) { has_escape = true; continue; }
          const int nn = neighbor(nr, cd.port);
          if (nn < 0) { ++missing_links; continue; }
          if (manhattan(nn, pk.dst) > dist && pk.inter < 0) ++nonminimal;
          if (cd.vc_mask & esc_) has_escape = true;
          for (int w = 0; w < V_; ++w)
            if ((cd.vc_mask >> w) & 1) {
              const int ch2 = chan(nr, cd.port, w);
              add_edge(ch, ch2);
              if (kSpec[alg_].escape && v == 0 && w != 0) ++esc_to_adaptive;
              if (stamp[ch2] != ctx_id) { stamp[ch2] = ctx_id; queue.push_back(ch2); }
            }
        }
        if (adaptive_state && !has_escape) ++no_escape;
      }
    };
    std::vector<int> sources;
    Packet pk;
    for (int d = 0; d < n_; ++d) {
      pk = Packet{};
      pk.dst = d;
      if (alg_ == kO1T) {
        for (int order = 0; order < 2; ++order) {
          pk.order = static_cast<uint8_t>(order);
          sources.clear();
          for (int s = 0; s < n_; ++s) if (s != d) sources.push_back(s);
          explore(pk, sources);
        }
      } else if (alg_ == kROMM || alg_ == kVal) {
        for (int m = 0; m < n_; ++m) {
          pk.inter = m;
          sources.clear();
          for (int s = 0; s < n_; ++s) if (s != d) sources.push_back(s);
          explore(pk, sources);
        }
      } else if (vpc_) {
        const int ne = static_cast<int>(elev_xy_.size());
        pk.dir = 0;
        pk.elevator = -1;
        sources.clear();
        for (int s = 0; s < n_; ++s) if (s != d && Z(s) == Z(d)) sources.push_back(s);
        explore(pk, sources);
        for (int dir = -1; dir <= 1; dir += 2)
          for (int e = 0; e < ne; ++e) {
            pk.dir = static_cast<int8_t>(dir);
            pk.elevator = e;
            sources.clear();
            for (int s = 0; s < n_; ++s)
              if ((dir > 0 && Z(s) < Z(d)) || (dir < 0 && Z(s) > Z(d))) sources.push_back(s);
            explore(pk, sources);
          }
      } else {
        sources.clear();
        for (int s = 0; s < n_; ++s) if (s != d) sources.push_back(s);
        explore(pk, sources);
      }
    }
    // Adjacency lists and cycle detection (iterative DFS, colours).
    std::vector<std::vector<int>> adj(C);
    long long edges = 0;
    for (int a = 0; a < C; ++a)
      for (int b = 0; b < C; ++b) {
        const size_t k = static_cast<size_t>(a) * C + b;
        if (bits[k >> 6] >> (k & 63) & 1ULL) { adj[a].push_back(b); ++edges; }
      }
    auto find_cycle = [&](const std::vector<char>& include, std::vector<int>& cycle) {
      std::vector<char> colour(C, 0);
      std::vector<int> parent(C, -1);
      for (int s = 0; s < C; ++s) {
        if (!include[s] || colour[s]) continue;
        std::vector<std::pair<int, size_t>> st{{s, 0}};
        colour[s] = 1;
        while (!st.empty()) {
          auto& [u, i] = st.back();
          if (i < adj[u].size()) {
            const int w = adj[u][i++];
            if (!include[w]) continue;
            if (colour[w] == 0) {
              colour[w] = 1;
              parent[w] = u;
              st.push_back({w, 0});
            } else if (colour[w] == 1) {
              cycle.clear();
              for (int x = u; x != w; x = parent[x]) cycle.push_back(x);
              cycle.push_back(w);
              std::reverse(cycle.begin(), cycle.end());
              return true;
            }
          } else {
            colour[u] = 2;
            st.pop_back();
          }
        }
      }
      return false;
    };
    std::vector<char> all_ch(C, 1), esc_ch(C, 0);
    for (int ch = 0; ch < C; ++ch) esc_ch[ch] = (ch % V_) == 0;
    // Shortest dependency cycle (girth) by breadth-first search from every
    // channel; used as the reported witness.
    auto shortest_cycle = [&](const std::vector<char>& include, std::vector<int>& best) {
      best.clear();
      std::vector<int> dist(C), par(C);
      for (int s = 0; s < C; ++s) {
        if (!include[s] || adj[s].empty()) continue;
        std::fill(dist.begin(), dist.end(), -1);
        std::vector<int> q{s};
        dist[s] = 0;
        bool found = false;
        for (size_t qi = 0; qi < q.size() && !found; ++qi) {
          const int u = q[qi];
          if (!best.empty() && dist[u] + 1 >= static_cast<int>(best.size())) break;
          for (int w : adj[u]) {
            if (!include[w]) continue;
            if (w == s) {
              std::vector<int> cyc;
              for (int x = u; x != s; x = par[x]) cyc.push_back(x);
              cyc.push_back(s);
              std::reverse(cyc.begin(), cyc.end());
              if (best.empty() || cyc.size() < best.size()) best = cyc;
              found = true;
              break;
            }
            if (dist[w] < 0) {
              dist[w] = dist[u] + 1;
              par[w] = u;
              q.push_back(w);
            }
          }
        }
      }
      return !best.empty();
    };
    std::vector<int> cycle, esc_cycle;
    bool cyclic = find_cycle(all_ch, cycle);
    if (cyclic) shortest_cycle(all_ch, cycle);
    bool esc_cyclic = false;
    if (kSpec[alg_].escape) esc_cyclic = find_cycle(esc_ch, esc_cycle);
    std::string verdict;
    if (!cyclic) verdict = "acyclic CDG (deadlock-free)";
    else if (kSpec[alg_].escape && !esc_cyclic && esc_to_adaptive == 0 && no_escape == 0)
      verdict = "cyclic CDG, acyclic connected escape subnetwork without return (deadlock-free, Duato)";
    else verdict = "cyclic CDG without escape subnetwork (deadlock possible)";
    std::printf("{\"routing\": \"%s\", \"topo\": \"%s\", \"channels_with_deps\": %d, \"edges\": %lld, "
                "\"states\": %lld, \"cyclic\": %s, \"escape\": %s, \"escape_cyclic\": %s, "
                "\"escape_to_adaptive\": %lld, \"adaptive_states_without_escape\": %lld, "
                "\"dead_ends\": %lld, \"missing_links\": %lld, \"nonminimal_moves\": %lld, "
                "\"verdict\": \"%s\", \"cycle\": [",
                kSpec[alg_].name, c_.topo.c_str(), count_used(adj), edges, states,
                cyclic ? "true" : "false", kSpec[alg_].escape ? "true" : "false",
                esc_cyclic ? "true" : "false", esc_to_adaptive, no_escape, dead_ends, missing_links,
                nonminimal, verdict.c_str());
    for (size_t i = 0; i < cycle.size() && i < 64; ++i) {
      const int ch = cycle[i];
      const int r = ch / (kPorts * V_), p = (ch / V_) % kPorts, v = ch % V_;
      std::printf("%s[%d, %d, %d, %d, %d]", i ? ", " : "", X(r), Y(r), Z(r), p, v);
    }
    std::printf("], \"cycle_length\": %zu}\n", cycle.size());
  }

 private:
  [[noreturn]] static void die(const char* msg) {
    std::fprintf(stderr, "error: %s\n", msg);
    std::exit(2);
  }

  static int count_used(const std::vector<std::vector<int>>& adj) {
    int k = 0;
    for (const auto& a : adj) k += !a.empty();
    return k;
  }

  // ----- geometry ---------------------------------------------------------
  int id(int x, int y, int z) const { return x + c_.kx * (y + c_.ky * z); }
  int X(int r) const { return xs_[r]; }
  int Y(int r) const { return ys_[r]; }
  int Z(int r) const { return zs_[r]; }
  int manhattan(int a, int b) const {
    return std::abs(X(a) - X(b)) + std::abs(Y(a) - Y(b)) + std::abs(Z(a) - Z(b));
  }
  int planar_dist(int x0, int y0, int x1, int y1) const { return std::abs(x0 - x1) + std::abs(y0 - y1); }

  InputVC& iv(int r, int p, int v) { return in_[(static_cast<size_t>(r) * kPorts + p) * V_ + v]; }
  OutputVC& ov(int r, int p, int v) { return outvc_[(static_cast<size_t>(r) * kPorts + p) * V_ + v]; }
  const OutputVC& ov(int r, int p, int v) const {
    return outvc_[(static_cast<size_t>(r) * kPorts + p) * V_ + v];
  }
  int neighbor(int r, int p) const { return nbr_[static_cast<size_t>(r) * kPorts + p]; }

  void build_topology() {
    for (int x = 0; x < c_.kx; ++x)
      for (int y = 0; y < c_.ky; ++y) {
        bool keep;
        if (!vpc_ || c_.elevators == "full") keep = true;
        else if (c_.elevators == "stride2") keep = (x % 2 == 1) && (y % 2 == 1);
        else if (c_.elevators == "stride4") keep = (x % 4 == 1) && (y % 4 == 1);
        else die("unknown elevator layout");
        if (keep) elev_xy_.push_back({x, y});
      }
    const int ne = static_cast<int>(elev_xy_.size());
    seg_ok_.assign(static_cast<size_t>(ne) * std::max(1, c_.kz - 1), 1);
    if (c_.pv > 0.0) sample_faults(ne);
    nbr_.assign(static_cast<size_t>(n_) * kPorts, -1);
    elev_at_.assign(c_.kx * c_.ky, -1);
    for (int e = 0; e < ne; ++e) elev_at_[elev_xy_[e][0] + c_.kx * elev_xy_[e][1]] = e;
    for (int r = 0; r < n_; ++r) {
      const int x = X(r), y = Y(r), z = Z(r);
      if (x + 1 < c_.kx) nbr_[r * kPorts + kEast] = id(x + 1, y, z);
      if (x > 0) nbr_[r * kPorts + kWest] = id(x - 1, y, z);
      if (y + 1 < c_.ky) nbr_[r * kPorts + kNorth] = id(x, y + 1, z);
      if (y > 0) nbr_[r * kPorts + kSouth] = id(x, y - 1, z);
      const int e = elev_at_[x + c_.kx * y];
      if (e >= 0) {
        if (z + 1 < c_.kz && seg_ok_[e * (c_.kz - 1) + z]) nbr_[r * kPorts + kUp] = id(x, y, z + 1);
        if (z > 0 && seg_ok_[e * (c_.kz - 1) + z - 1]) nbr_[r * kPorts + kDown] = id(x, y, z - 1);
      }
    }
  }

  // Snake labelling: layer z traversed boustrophedon, reversed on odd layers,
  // so that consecutive labels are always adjacent routers.
  void build_hamiltonian() {
    ham_.assign(n_, 0);
    const int plane = c_.kx * c_.ky;
    for (int r = 0; r < n_; ++r) {
      const int x = X(r), y = Y(r), z = Z(r);
      int idx = (y % 2 == 0) ? y * c_.kx + x : y * c_.kx + (c_.kx - 1 - x);
      if (z % 2 == 1) idx = plane - 1 - idx;
      ham_[r] = z * plane + idx;
    }
  }

  // Gaussian copula with exponential spatial covariance and exact marginal pv.
  void sample_faults(int ne) {
    std::vector<double> L(static_cast<size_t>(ne) * ne, 0.0), S(static_cast<size_t>(ne) * ne, 0.0);
    for (int i = 0; i < ne; ++i)
      for (int j = 0; j < ne; ++j) {
        const double dx = elev_xy_[i][0] - elev_xy_[j][0], dy = elev_xy_[i][1] - elev_xy_[j][1];
        const double d = std::sqrt(dx * dx + dy * dy);
        S[i * ne + j] = (c_.xi > 0.0) ? std::exp(-d / c_.xi) : (i == j ? 1.0 : 0.0);
      }
    for (int i = 0; i < ne; ++i)
      for (int j = 0; j <= i; ++j) {
        double s = S[i * ne + j];
        for (int k = 0; k < j; ++k) s -= L[i * ne + k] * L[j * ne + k];
        L[i * ne + j] = (i == j) ? std::sqrt(std::max(s, 1e-12)) : s / L[j * ne + j];
      }
    std::mt19937_64 frng(c_.seed * 7919ULL + 17ULL);
    std::normal_distribution<double> g(0.0, 1.0);
    const double thr = inverse_normal_cdf(c_.pv);
    for (int b = 0; b < c_.kz - 1; ++b) {
      std::vector<double> u(ne), zl(ne, 0.0);
      for (auto& v : u) v = g(frng);
      for (int i = 0; i < ne; ++i)
        for (int k = 0; k <= i; ++k) zl[i] += L[i * ne + k] * u[k];
      for (int e = 0; e < ne; ++e) seg_ok_[e * (c_.kz - 1) + b] = zl[e] >= thr;
    }
  }

  static double inverse_normal_cdf(double p) {
    // Acklam's rational approximation, |relative error| < 1.2e-9.
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
    static const double cc[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                                -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                               3.754408661907416e+00};
    const double plow = 0.02425;
    if (p < plow) {
      const double q = std::sqrt(-2 * std::log(p));
      return (((((cc[0] * q + cc[1]) * q + cc[2]) * q + cc[3]) * q + cc[4]) * q + cc[5]) /
             ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    if (p > 1 - plow) return -inverse_normal_cdf(1 - p);
    const double q = p - 0.5, r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }

  bool column_spans(int e, int z0, int z1) const {
    const int lo = std::min(z0, z1), hi = std::max(z0, z1);
    for (int b = lo; b < hi; ++b)
      if (!seg_ok_[e * (c_.kz - 1) + b]) return false;
    return true;
  }

  // ----- traffic ----------------------------------------------------------
  void init_traffic() {
    hotspot_ = id(c_.kx / 2 - 1, c_.ky / 2 - 1, std::max(0, c_.kz / 2 - 1));
    attacker_.assign(n_, 0);
    if (c_.attackers > 0) {
      std::mt19937_64 arng(c_.seed * 101ULL + 7ULL);
      std::vector<int> pool;
      for (int r = 0; r < n_; ++r) if (r != hotspot_) pool.push_back(r);
      std::shuffle(pool.begin(), pool.end(), arng);
      for (int i = 0; i < c_.attackers && i < static_cast<int>(pool.size()); ++i) attacker_[pool[i]] = 1;
    }
    if (c_.traffic == "selfsim") {
      p_on_ = std::min(1.0, c_.burst_factor * c_.pir);
      const double duty = c_.pir / p_on_;
      const double mean_on = c_.burst_on_mean;
      const double mean_off = mean_on * (1.0 - duty) / duty;
      xm_on_ = mean_on * (c_.pareto_shape - 1.0) / c_.pareto_shape;
      xm_off_ = mean_off * (c_.pareto_shape - 1.0) / c_.pareto_shape;
      on_.assign(n_, 0);
      left_.assign(n_, 0);
      for (int s = 0; s < n_; ++s) {
        on_[s] = uni_(trng_) < duty;
        left_[s] = pareto(on_[s] ? xm_on_ : xm_off_);
      }
    }
  }

  int pareto(double xm) {
    const double u = std::max(1e-12, 1.0 - uni_(trng_));
    return std::max(1, static_cast<int>(std::lround(xm / std::pow(u, 1.0 / c_.pareto_shape))));
  }

  bool onoff(int s) {
    if (left_[s] <= 0) {
      on_[s] = !on_[s];
      left_[s] = pareto(on_[s] ? xm_on_ : xm_off_);
    }
    --left_[s];
    return on_[s] && uni_(trng_) < p_on_;
  }

  int destination(int s) {
    const int bits = static_cast<int>(std::lround(std::log2(n_)));
    std::uniform_int_distribution<int> any(0, n_ - 1);
    auto uniform = [&]() {
      int d;
      do { d = any(trng_); } while (d == s);
      return d;
    };
    const std::string& tp = c_.traffic;
    if (tp == "random" || tp == "selfsim" || tp == "dos") return uniform();
    if (tp == "hotspot") {
      if (uni_(trng_) < c_.hotspot_frac && s != hotspot_) return hotspot_;
      return uniform();
    }
    if (tp == "transpose1") return id(c_.kx - 1 - Y(s), c_.ky - 1 - X(s), c_.kz - 1 - Z(s));
    if ((1 << bits) != n_) die("bit-permutation traffic needs a power-of-two node count");
    if (tp == "shuffle") return ((s << 1) | (s >> (bits - 1))) & (n_ - 1);
    if (tp == "bitreversal") {
      int d = 0;
      for (int i = 0; i < bits; ++i)
        if (s & (1 << i)) d |= 1 << (bits - 1 - i);
      return d;
    }
    if (tp == "butterfly") {
      const int msb = (s >> (bits - 1)) & 1, lsb = s & 1;
      const int d = s & ~(1 | (1 << (bits - 1)));
      return d | (lsb << (bits - 1)) | msb;
    }
    die("unknown traffic pattern");
  }

  // ----- elevator selection (VPC) -----------------------------------------
  int nearest_elevator(int s, bool healthy, int zd, bool prefer_negative) const {
    const int ne = static_cast<int>(elev_xy_.size());
    int best = -1, best_d = 1 << 30;
    bool best_neg = false;
    for (int e = 0; e < ne; ++e) {
      if (healthy && !column_spans(e, Z(s), zd)) continue;
      const int ex = elev_xy_[e][0], ey = elev_xy_[e][1];
      const int dd = planar_dist(X(s), Y(s), ex, ey);
      const bool neg = ex <= X(s) && ey <= Y(s);
      if (dd < best_d || (prefer_negative && dd == best_d && neg && !best_neg)) {
        best_d = dd;
        best = e;
        best_neg = neg;
      }
    }
    return best;
  }

  void init_adele() {
    const int ne = static_cast<int>(elev_xy_.size());
    const int k = std::min(c_.adele_k, ne);
    adele_subset_.assign(n_, {});
    adele_cost_.assign(n_, std::vector<double>(k, 0.0));
    adele_ptr_.assign(n_, 0);
    for (int s = 0; s < n_; ++s) {
      std::vector<std::pair<int, int>> order;
      for (int e = 0; e < ne; ++e)
        order.push_back({planar_dist(X(s), Y(s), elev_xy_[e][0], elev_xy_[e][1]), e});
      std::stable_sort(order.begin(), order.end());
      for (int i = 0; i < k; ++i) adele_subset_[s].push_back(order[i].second);
    }
  }

  int adele_choose(int s, int d) {
    const auto& sub = adele_subset_[s];
    auto& cost = adele_cost_[s];
    const int k = static_cast<int>(sub.size());
    bool low = true;
    for (double c : cost) low = low && c < c_.adele_thr;
    if (low) {
      int best = -1, best_d = 1 << 30;
      for (int i = 0; i < k; ++i) {
        const int ex = elev_xy_[sub[i]][0], ey = elev_xy_[sub[i]][1];
        const int dd = planar_dist(X(s), Y(s), ex, ey) + planar_dist(ex, ey, X(d), Y(d));
        if (dd < best_d) { best_d = dd; best = i; }
      }
      return sub[best];
    }
    double total = 0.0;
    for (double c : cost) total += c;
    for (int pass = 0; pass < 16; ++pass)
      for (int j = 0; j < k; ++j) {
        const int i = (adele_ptr_[s] + j) % k;
        const double rel = total > 0.0 ? cost[i] / total : 1.0 / k;
        double ps;
        if (rel >= 2.0 / k) ps = 1.0 - c_.adele_xi;
        else if (rel >= 1.0 / k) ps = k * (rel - 1.0 / k) * (1.0 - c_.adele_xi);
        else ps = 0.0;
        if (uni_(rng_) >= ps) {
          adele_ptr_[s] = (i + 1) % k;
          return sub[i];
        }
      }
    const int i = static_cast<int>(std::min_element(cost.begin(), cost.end()) - cost.begin());
    return sub[i];
  }

  void adele_feedback(const Packet& pk, int t_tail) {
    const auto& sub = adele_subset_[pk.src];
    for (size_t i = 0; i < sub.size(); ++i)
      if (sub[i] == pk.elevator) {
        const double te = std::max(0.0, (t_tail - pk.t_head_out + 1.0 - pk.len) / pk.len);
        double& c = adele_cost_[pk.src][i];
        c = c_.adele_a * te + (1.0 - c_.adele_a) * c;
        return;
      }
  }

  int choose_elevator(int s, int d) {
    const int ne = static_cast<int>(elev_xy_.size());
    switch (alg_) {
      case kEF: {
        const int e = nearest_elevator(s, false, Z(d), false);
        return column_spans(e, Z(s), Z(d)) ? e : -1;
      }
      case kEFR: return nearest_elevator(s, true, Z(d), false);
      case kFL: {
        const int e = nearest_elevator(s, false, Z(d), true);
        return column_spans(e, Z(s), Z(d)) ? e : -1;
      }
      case kAdEle: {
        const int e = adele_choose(s, d);
        return column_spans(e, Z(s), Z(d)) ? e : -1;
      }
      case kAESG: {
        int best = -1;
        double best_cost = 1e30;
        for (int e = 0; e < ne; ++e) {
          if (!column_spans(e, Z(s), Z(d))) continue;
          const int ex = elev_xy_[e][0], ey = elev_xy_[e][1];
          const double cost = planar_dist(X(s), Y(s), ex, ey) + planar_dist(ex, ey, X(d), Y(d)) +
                              c_.aes_beta * elev_load_[e];
          if (cost < best_cost - 1e-12) { best_cost = cost; best = e; }
        }
        return best;
      }
      default: die("not a VPC class");
    }
  }

  void create_packet(int s, int d, int t, bool attacker) {
    Packet pk;
    pk.src = s;
    pk.dst = d;
    pk.t_create = t;
    pk.len = static_cast<int16_t>(c_.pkt_flits);
    pk.attacker = attacker;
    if (alg_ == kO1T) pk.order = uni_(rng_) < 0.5 ? 0 : 1;
    if (c_.inject_escape && kSpec[alg_].escape) pk.order = uni_(rng_) < 0.5 ? 0 : 1;
    if (alg_ == kVal) pk.inter = std::uniform_int_distribution<int>(0, n_ - 1)(rng_);
    if (alg_ == kROMM) {
      auto pick = [&](int a, int b) {
        return std::uniform_int_distribution<int>(std::min(a, b), std::max(a, b))(rng_);
      };
      pk.inter = id(pick(X(s), X(d)), pick(Y(s), Y(d)), pick(Z(s), Z(d)));
    }
    const int pid = static_cast<int>(pkts_.size());
    if (vpc_) {
      pk.dir = static_cast<int8_t>(Z(d) > Z(s) ? 1 : (Z(d) < Z(s) ? -1 : 0));
      if (pk.dir != 0) {
        pk.elevator = choose_elevator(s, d);
        if (pk.elevator < 0) {
          pk.dropped = true;
          pkts_.push_back(pk);
          log(t, pid, s, -1, kDrop, 0, 0);
          return;
        }
        ++elev_load_[pk.elevator];
      }
    }
    pkts_.push_back(pk);
    srcq_[s].push_back(pid);
  }

  void generate(int t) {
    for (int s = 0; s < n_; ++s) {
      if (attacker_[s]) {
        if (uni_(trng_) < c_.attack_pir) create_packet(s, hotspot_, t, true);
        continue;
      }
      const bool fire = (c_.traffic == "selfsim") ? onoff(s) : uni_(trng_) < c_.pir;
      if (!fire) continue;
      const int d = destination(s);
      if (d == s) continue;
      create_packet(s, d, t, false);
    }
  }

  // ----- routing relation -------------------------------------------------
  static int dor(int dx, int dy, int dz, const int* order) {
    for (int k = 0; k < 3; ++k) {
      const int dim = order[k];
      const int delta = dim == 0 ? dx : (dim == 1 ? dy : dz);
      if (delta == 0) continue;
      if (dim == 0) return delta > 0 ? kEast : kWest;
      if (dim == 1) return delta > 0 ? kNorth : kSouth;
      return delta > 0 ? kUp : kDown;
    }
    return kLocal;
  }
  int dor_to(int r, int d, const int* order) const {
    return dor(X(d) - X(r), Y(d) - Y(r), Z(d) - Z(r), order);
  }
  static constexpr int kOrdXYZ[3] = {0, 1, 2};
  static constexpr int kOrdZXY[3] = {2, 0, 1};
  static constexpr int kOrdZYX[3] = {2, 1, 0};

  void minimal_ports(int r, int d, std::vector<Cand>& out, int mask) const {
    if (X(d) != X(r)) out.push_back({X(d) > X(r) ? kEast : kWest, mask});
    if (Y(d) != Y(r)) out.push_back({Y(d) > Y(r) ? kNorth : kSouth, mask});
    if (Z(d) != Z(r)) out.push_back({Z(d) > Z(r) ? kUp : kDown, mask});
  }

  void nf_ports(int r, int d, std::vector<Cand>& out, int mask) const {
    const size_t before = out.size();
    if (X(d) < X(r)) out.push_back({kWest, mask});
    if (Y(d) < Y(r)) out.push_back({kSouth, mask});
    if (Z(d) < Z(r)) out.push_back({kDown, mask});
    if (out.size() > before) return;
    if (X(d) > X(r)) out.push_back({kEast, mask});
    if (Y(d) > Y(r)) out.push_back({kNorth, mask});
    if (Z(d) > Z(r)) out.push_back({kUp, mask});
  }

  void ham_ports(int r, int d, std::vector<Cand>& out) const {
    const int lr = ham_[r], ld = ham_[d];
    for (int p = 1; p < kPorts; ++p) {
      const int nb = neighbor(r, p);
      if (nb < 0) continue;
      const int ln = ham_[nb];
      const bool ok = ld > lr ? (ln > lr && ln <= ld) : (ln < lr && ln >= ld);
      if (ok) out.push_back({p, all_});
    }
  }

  void ef_relation(int r, int ip, int v, const Packet& pk, std::vector<Cand>& out) const {
    const int d = pk.dst;
    int mask;
    if (pk.dir > 0) mask = c0_;
    else if (pk.dir < 0) mask = c1_;
    else mask = (ip == kLocal) ? all_ : (((c0_ >> v) & 1) ? c0_ : c1_);
#ifdef NEGATIVE_TEST_EF_MIXED
    // Known-unsafe variant used only to test the checker: same-layer packets
    // may change virtual network at every hop.
    if (pk.dir == 0) mask = all_;
#endif
    if (Z(r) == Z(d)) {
      out.push_back({dor_to(r, d, kOrdXYZ), mask});
      return;
    }
    const int ex = elev_xy_[pk.elevator][0], ey = elev_xy_[pk.elevator][1];
    int p;
    if (X(r) != ex) p = ex > X(r) ? kEast : kWest;
    else if (Y(r) != ey) p = ey > Y(r) ? kNorth : kSouth;
    else p = Z(d) > Z(r) ? kUp : kDown;
    out.push_back({p, mask});
  }

  void fl_relation(int r, int ip, int v, const Packet& pk, std::vector<Cand>& out) const {
    const int d = pk.dst;
    int vn = -1;
    if (ip == kWest || ip == kSouth) vn = ((c0_ >> v) & 1) ? 0 : 2;
    else if (ip != kLocal) vn = 1;
    const bool other_tier = Z(r) != Z(d);
    const int tx = other_tier ? elev_xy_[pk.elevator][0] : X(d);
    const int ty = other_tier ? elev_xy_[pk.elevator][1] : Y(d);
    const bool need_e = tx > X(r), need_n = ty > Y(r), need_w = tx < X(r), need_s = ty < Y(r);
    if (other_tier) {
      if (!need_e && !need_n && !need_w && !need_s) {
        out.push_back({Z(d) > Z(r) ? kUp : kDown, all_});
        return;
      }
      if (vn <= 0 && (need_e || need_n)) {
        if (need_e) out.push_back({kEast, c0_});
        if (need_n) out.push_back({kNorth, c0_});
        return;
      }
      if (need_w) out.push_back({kWest, all_});
      if (need_s) out.push_back({kSouth, all_});
      return;
    }
    if ((need_w || need_s) && vn <= 1) {
      if (need_w) out.push_back({kWest, all_});
      if (need_s) out.push_back({kSouth, all_});
      return;
    }
    if (need_e) out.push_back({kEast, c1_});
    if (need_n) out.push_back({kNorth, c1_});
  }

  // Admissible (port, VC mask) pairs for the head of pk at router r, input
  // port ip and VC v. Shared by the simulator and the CDG checker.
  void relation(int r, int ip, int v, const Packet& pk, std::vector<Cand>& out) const {
    out.clear();
    const int d = pk.dst;
    if (r == d) {
      out.push_back({kLocal, all_});
      return;
    }
    switch (alg_) {
      case kXYZ: out.push_back({dor_to(r, d, kOrdXYZ), all_}); return;
      case kZXY: out.push_back({dor_to(r, d, kOrdZXY), all_}); return;
      case kNF:
      case kDyAD: nf_ports(r, d, out, all_); return;
      case kHam: ham_ports(r, d, out); return;
      case kO1T: {
        const int order = (ip == kLocal) ? pk.order : ((c1_ >> v) & 1);
        if (order == 0) out.push_back({dor_to(r, d, kOrdXYZ), c0_});
        else out.push_back({dor_to(r, d, kOrdZYX), c1_});
        return;
      }
      case kROMM:
      case kVal: {
        const bool phase1 = (r == pk.inter) || (ip != kLocal && ((c1_ >> v) & 1));
        if (phase1) out.push_back({dor_to(r, d, kOrdXYZ), c1_});
        else out.push_back({dor_to(r, pk.inter, kOrdXYZ), c0_});
        return;
      }
      case kFA:
      case kRCA:
      case kQR:
      case kTA: {
        const bool escaped = (ip != kLocal && v == 0) || (c_.inject_escape && ip == kLocal && pk.order == 1);
        if (!escaped) minimal_ports(r, d, out, adp_);
        out.push_back({dor_to(r, d, kOrdXYZ), esc_});
        return;
      }
      case kEF:
      case kEFR:
      case kAESG:
      case kAdEle: ef_relation(r, ip, v, pk, out); return;
      case kFL: fl_relation(r, ip, v, pk, out); return;
      default: die("routing relation not defined");
    }
  }

  int port_credits(int r, int p) const {
    int s = 0;
    for (int v = 0; v < V_; ++v) s += std::min(ov(r, p, v).credits, c_.depth);
    return s;
  }

  // Candidate set used by the simulator: a subset of relation() that encodes
  // class-specific restrictions of the selection function.
  void candidates(int r, int ip, int v, const Packet& pk, std::vector<Cand>& out) {
    relation(r, ip, v, pk, out);
    if (out.empty() || out[0].port == kLocal) return;
    switch (alg_) {
      case kDyAD: {
        const int det = out.front().port;
        if (port_credits(r, det) >= c_.dyad_frac * V_ * c_.depth) out.erase(out.begin() + 1, out.end());
        return;
      }
      case kHam: {
        const int dist = manhattan(r, pk.dst);
        tmp_.clear();
        for (const Cand& cd : out)
          if (manhattan(neighbor(r, cd.port), pk.dst) < dist) tmp_.push_back(cd);
        if (!tmp_.empty()) {
          out = tmp_;
          return;
        }
        Cand best = out.front();
        int best_gap = 1 << 30;
        for (const Cand& cd : out) {
          const int gap = std::abs(ham_[pk.dst] - ham_[neighbor(r, cd.port)]);
          if (gap < best_gap) { best_gap = gap; best = cd; }
        }
        out.assign(1, best);
        return;
      }
      default: return;
    }
  }

  // Selection key of a free output VC, compared lexicographically.
  struct Key {
    int tier;
    double k1;
    double k2;
    bool operator>(const Key& o) const {
      if (tier != o.tier) return tier > o.tier;
      if (k1 != o.k1) return k1 > o.k1;
      return k2 > o.k2;
    }
  };

  Key select_key(int r, int port, int vc, const Packet& pk) {
    const double credits = port_credits(r, port);
    if (!kSpec[alg_].escape || vc == 0 || port == kLocal) return {0, credits, 0.0};
    if (c_.rand_sel && alg_ == kFA) return {1, uni_(rng_), 0.0};
    switch (alg_) {
      case kRCA: return {1, -agg_[static_cast<size_t>(r) * kPorts + port], credits};
      case kQR: return {1, -q_[qidx(r, pk.dst, port)], credits};
      case kTA: return {1, -temp_[neighbor(r, port)], credits};
      default: return {1, credits, 0.0};
    }
  }

  uint64_t port_signature(int r, const InputVC& in) const {
    uint64_t sig = 0;
    for (int i = 0; i < in.n_cand_ports; ++i) sig += pver_[static_cast<size_t>(r) * kPorts + in.cand_ports[i]];
    return sig;
  }

  bool output_free(int r, int p, int v) const {
    const OutputVC& o = ov(r, p, v);
    if (o.owner_port >= 0) return false;
    return p == kLocal || !c_.atomic || o.credits >= c_.depth;
  }

  bool allocate(int r, int ip, int v) {
    InputVC& in = iv(r, ip, v);
    const Packet& pk = pkts_[in.q.front().pkt];
    candidates(r, ip, v, pk, cands_);
    int best_p = -1, best_v = -1;
    Key best{-1, -1e300, -1e300};
    for (const Cand& cd : cands_) {
      if (cd.port != kLocal && neighbor(r, cd.port) < 0) continue;
      for (int ovc = 0; ovc < V_; ++ovc) {
        if (!((cd.vc_mask >> ovc) & 1)) continue;
        if (!output_free(r, cd.port, ovc)) continue;
        const Key k = select_key(r, cd.port, ovc, pk);
        if (k > best) {
          best = k;
          best_p = cd.port;
          best_v = ovc;
        }
      }
    }
    if (best_p < 0) {
      if (cacheable_) {
        in.n_cand_ports = 0;
        for (const Cand& cd : cands_)
          if (in.n_cand_ports < 8) in.cand_ports[in.n_cand_ports++] = static_cast<uint8_t>(cd.port);
        in.fail_sig = port_signature(r, in);
        in.failed = true;
      }
      return false;
    }
    in.failed = false;
    in.out_port = best_p;
    in.out_vc = best_v;
    ov(r, best_p, best_v).owner_port = ip;
    ov(r, best_p, best_v).owner_vc = v;
    return true;
  }

  // ----- class-specific state ----------------------------------------------
  size_t qidx(int r, int d, int p) const {
    return (static_cast<size_t>(r) * n_ + d) * kPorts + p;
  }

  void init_q() {
    q_.assign(static_cast<size_t>(n_) * n_ * kPorts, 0.0f);
    for (int r = 0; r < n_; ++r)
      for (int d = 0; d < n_; ++d)
        for (int p = 0; p < kPorts; ++p)
          q_[qidx(r, d, p)] = static_cast<float>(c_.hop_delay * manhattan(r, d));
  }

  float q_min(int r, int d) const {
    if (r == d) return 0.0f;
    float best = 1e30f;
    if (X(d) != X(r)) best = std::min(best, q_[qidx(r, d, X(d) > X(r) ? kEast : kWest)]);
    if (Y(d) != Y(r)) best = std::min(best, q_[qidx(r, d, Y(d) > Y(r) ? kNorth : kSouth)]);
    if (Z(d) != Z(r)) best = std::min(best, q_[qidx(r, d, Z(d) > Z(r) ? kUp : kDown)]);
    return best;
  }

  void init_thermal() {
    // Steady-state compact model G T = P with lateral, inter-layer and sink
    // conductances; R = G^{-1} by Cholesky factorisation.
    std::vector<double> G(static_cast<size_t>(n_) * n_, 0.0);
    for (int r = 0; r < n_; ++r) {
      for (int p = 1; p < kPorts; ++p) {
        const int x = X(r), y = Y(r), z = Z(r);
        int nb = -1;
        double g = 0.0;
        if (p == kEast && x + 1 < c_.kx) { nb = id(x + 1, y, z); g = c_.g_lat; }
        if (p == kWest && x > 0) { nb = id(x - 1, y, z); g = c_.g_lat; }
        if (p == kNorth && y + 1 < c_.ky) { nb = id(x, y + 1, z); g = c_.g_lat; }
        if (p == kSouth && y > 0) { nb = id(x, y - 1, z); g = c_.g_lat; }
        if (p == kUp && z + 1 < c_.kz) { nb = id(x, y, z + 1); g = c_.g_vert; }
        if (p == kDown && z > 0) { nb = id(x, y, z - 1); g = c_.g_vert; }
        if (nb < 0) continue;
        G[static_cast<size_t>(r) * n_ + r] += g;
        G[static_cast<size_t>(r) * n_ + nb] -= g;
      }
      if (Z(r) == 0) G[static_cast<size_t>(r) * n_ + r] += c_.g_sink;
    }
    std::vector<double> L(static_cast<size_t>(n_) * n_, 0.0);
    for (int i = 0; i < n_; ++i)
      for (int j = 0; j <= i; ++j) {
        double s = G[static_cast<size_t>(i) * n_ + j];
        for (int k = 0; k < j; ++k) s -= L[static_cast<size_t>(i) * n_ + k] * L[static_cast<size_t>(j) * n_ + k];
        L[static_cast<size_t>(i) * n_ + j] = (i == j) ? std::sqrt(s) : s / L[static_cast<size_t>(j) * n_ + j];
      }
    rth_.assign(static_cast<size_t>(n_) * n_, 0.0);
    std::vector<double> col(n_), y(n_);
    for (int c = 0; c < n_; ++c) {
      for (int i = 0; i < n_; ++i) {
        double s = (i == c) ? 1.0 : 0.0;
        for (int k = 0; k < i; ++k) s -= L[static_cast<size_t>(i) * n_ + k] * y[k];
        y[i] = s / L[static_cast<size_t>(i) * n_ + i];
      }
      for (int i = n_ - 1; i >= 0; --i) {
        double s = y[i];
        for (int k = i + 1; k < n_; ++k) s -= L[static_cast<size_t>(k) * n_ + i] * col[k];
        col[i] = s / L[static_cast<size_t>(i) * n_ + i];
      }
      for (int i = 0; i < n_; ++i) rth_[static_cast<size_t>(i) * n_ + c] = col[i];
    }
    temp_.assign(n_, 0.0);
    update_temperature(1);
  }

  void update_temperature(int period) {
    std::vector<double> power(n_);
    for (int r = 0; r < n_; ++r)
      power[r] = c_.p_core + c_.p_leak + c_.e_flit * c_.freq * act_period_[r] / std::max(1, period);
    for (int i = 0; i < n_; ++i) {
      double s = 0.0;
      for (int j = 0; j < n_; ++j) s += rth_[static_cast<size_t>(i) * n_ + j] * power[j];
      temp_[i] = s;
    }
    std::fill(act_period_.begin(), act_period_.end(), 0);
  }

  void update_rca() {
    const double cap = static_cast<double>(V_) * c_.depth;
    for (int r = 0; r < n_; ++r)
      for (int p = 1; p < kPorts; ++p) {
        const int nb = neighbor(r, p);
        if (nb < 0) { agg_next_[static_cast<size_t>(r) * kPorts + p] = 1.0; continue; }
        const double local = 1.0 - port_credits(r, p) / cap;
        const double remote = neighbor(nb, p) >= 0 ? agg_[static_cast<size_t>(nb) * kPorts + p] : 0.0;
        agg_next_[static_cast<size_t>(r) * kPorts + p] = c_.rca_weight * local + (1.0 - c_.rca_weight) * remote;
      }
    agg_.swap(agg_next_);
  }

  // ----- simulation step ----------------------------------------------------
  bool in_window(int t) const { return t >= c_.warmup && t < c_.cycles; }

  long long partial_injections() const {
    long long k = 0;
    for (int s = 0; s < n_; ++s)
      if (inj_flit_[s] > 0) k += pkts_[srcq_[s].front()].len - inj_flit_[s];
    return k;
  }

  void inject(int t, bool draining) {
    for (int s = 0; s < n_; ++s) {
      if (srcq_[s].empty()) continue;
      if (draining && inj_flit_[s] == 0) continue;
      const int pid = srcq_[s].front();
      if (inj_vc_[s] < 0) {
        for (int v = 0; v < V_; ++v) {
          InputVC& in = iv(s, kLocal, v);
          if (in.q.empty() && in.out_port < 0) { inj_vc_[s] = v; break; }
        }
        if (inj_vc_[s] < 0) continue;
      }
      InputVC& in = iv(s, kLocal, inj_vc_[s]);
      if (in.q.size() >= c_.depth) continue;
      const int k = inj_flit_[s];
      const int len = pkts_[pid].len;
      in.q.push(Flit{pid, t, k == 0, k == len - 1});
      ++occ_[s];
      ++in_flight_flits_;
      last_move_ = t;
      if (in_window(t)) ++act_[static_cast<size_t>(s) * 4 + 3];
      log(t, pid, s, -1, kInject, 0, kLocal);
      if (++inj_flit_[s] == len) {
        inj_flit_[s] = 0;
        inj_vc_[s] = -1;
        srcq_[s].pop_front();
      }
    }
  }

  struct Move {
    int r, ip, v;
  };

  void step(int t, bool draining) {
    now_ = t;
    if (alg_ == kRCA) update_rca();
    if (alg_ == kTA && t > 0 && t % c_.ta_period == 0) update_temperature(c_.ta_period);
    if (!draining) generate(t);
    inject(t, draining);
    moves_.clear();
    const int nvc = kPorts * V_;
    for (int r = 0; r < n_; ++r) {
      const int start = rot_[r];
      rot_[r] = (rot_[r] + 1) % nvc;
      if (occ_[r] == 0) continue;
      for (int k = 0; k < nvc; ++k) {
        int idx = start + k;
        if (idx >= nvc) idx -= nvc;
        const int ip = idx / V_, v = idx - ip * V_;
        InputVC& in = iv(r, ip, v);
        if (in.q.empty() || in.out_port >= 0) continue;
        const Flit& f = in.q.front();
        if (!f.head || f.ready > t) continue;
        if (in.wait_from < 0) in.wait_from = t;
        if (cacheable_ && in.failed && port_signature(r, in) == in.fail_sig) continue;
        if (allocate(r, ip, v)) {
          if (in_window(t)) {
            const int w = t - in.wait_from;
            max_va_wait_ = std::max(max_va_wait_, w);
            if (w > 1000) ++va_wait_gt1000_;
          }
          in.wait_from = -1;
        }
      }
      std::array<int, kPorts> req_vc;
      std::array<uint32_t, kPorts> req_mask{};
      req_vc.fill(-1);
      for (int ip = 0; ip < kPorts; ++ip) {
        const int ptr = rr_in_[r * kPorts + ip];
        for (int k = 0; k < V_; ++k) {
          int v = ptr + k;
          if (v >= V_) v -= V_;
          InputVC& in = iv(r, ip, v);
          if (in.q.empty() || in.out_port < 0) continue;
          if (in.q.front().ready > t) continue;
          if (in.out_port != kLocal && ov(r, in.out_port, in.out_vc).credits <= 0) continue;
          req_vc[ip] = v;
          req_mask[in.out_port] |= 1u << ip;
          break;
        }
      }
      for (int op = 0; op < kPorts; ++op) {
        const uint32_t m = req_mask[op];
        if (!m) continue;
        int& ptr = rr_out_[r * kPorts + op];
        // First requesting input at or after the round-robin pointer.
        const uint32_t rot = ((m >> ptr) | (m << (kPorts - ptr))) & ((1u << kPorts) - 1u);
        int ip = ptr + __builtin_ctz(rot);
        if (ip >= kPorts) ip -= kPorts;
        moves_.push_back({r, ip, req_vc[ip]});
        int nv = req_vc[ip] + 1;
        if (nv >= V_) nv -= V_;
        rr_in_[r * kPorts + ip] = nv;
        ptr = ip + 1 == kPorts ? 0 : ip + 1;
      }
    }
    auto& due = credit_pipe_[t % (c_.credit_delay + 1)];
    for (const Move& m : moves_) {
      InputVC& in = iv(m.r, m.ip, m.v);
      const Flit f = in.q.front();
      in.q.pop();
      --occ_[m.r];
      if (m.ip != kLocal) credit_pipe_[(t + c_.credit_delay) % (c_.credit_delay + 1)].push_back(
          {neighbor(m.r, m.ip), kOpposite[m.ip], m.v});
      const int op = in.out_port, ovc = in.out_vc;
      Packet& pk = pkts_[f.pkt];
      last_move_ = t;
      const bool win = in_window(t);
      if (win) ++act_[static_cast<size_t>(m.r) * 4 + 0];
      ++act_period_[m.r];
      if (m.ip == kLocal && f.head) pk.t_head_out = t;
      if (m.ip == kLocal && f.tail && alg_ == kAdEle && pk.elevator >= 0) adele_feedback(pk, t);
      if (op == kLocal) {
        --in_flight_flits_;
        log(t, f.pkt, m.r, -1, kReceive, 0, kLocal);
        if (f.tail) {
          pk.t_eject = t;
          if (pk.elevator >= 0) --elev_load_[pk.elevator];
        }
      } else {
        const int nr = neighbor(m.r, op);
        const bool vertical = op == kUp || op == kDown;
        ov(m.r, op, ovc).credits -= 1;
        iv(nr, kOpposite[op], ovc).q.push(Flit{f.pkt, t + c_.hop_delay, f.head, f.tail});
        ++occ_[nr];
        log(t, f.pkt, m.r, nr, kForward, vertical ? 1 : 0, op);
        if (win) {
          ++act_[static_cast<size_t>(m.r) * 4 + (vertical ? 2 : 1)];
          ++act_[static_cast<size_t>(nr) * 4 + 3];
          if (vertical) ++flit_hops_v_; else ++flit_hops_h_;
        }
        if (alg_ == kQR && f.head) {
          const float target = static_cast<float>(t - f.ready + c_.hop_delay) + q_min(nr, pk.dst);
          float& q = q_[qidx(m.r, pk.dst, op)];
          q += static_cast<float>(c_.q_eta) * (target - q);
          if (win) ++q_updates_;
        }
        if (f.head) {
          ++pk.hops;
          if (vertical) ++pk.vhops;
        }
      }
      if (f.tail) {
        ov(m.r, op, ovc).owner_port = -1;
        ov(m.r, op, ovc).owner_vc = -1;
        ++pver_[static_cast<size_t>(m.r) * kPorts + op];
        in.out_port = -1;
        in.out_vc = -1;
      }
    }
    for (const CreditMsg& cr : due) {
      OutputVC& o = ov(cr.r, cr.port, cr.vc);
      o.credits += 1;
      if (o.credits >= c_.depth) ++pver_[static_cast<size_t>(cr.r) * kPorts + cr.port];
    }
    due.clear();
    if (in_window(t)) {
      const int w = (t - c_.warmup) / c_.window;
      double* acc = &buf_acc_[static_cast<size_t>(w) * n_ * kPorts];
      for (int r = 0; r < n_; ++r)
        for (int p = 0; p < kPorts; ++p) {
          int occ = 0;
          for (int v = 0; v < V_; ++v) occ += iv(r, p, v).q.size();
          acc[r * kPorts + p] += occ;
        }
    }
  }

  void log(int t, int p, int r, int next, Event xi, int c, int port) {
    if (!trace_) return;
    const TraceRec rec{t, p, static_cast<int16_t>(r), static_cast<int16_t>(next), xi, static_cast<uint8_t>(c),
                       static_cast<uint16_t>(port)};
    std::fwrite(&rec, sizeof(rec), 1, trace_);
  }

  void write_outputs() {
    std::FILE* fb = std::fopen((c_.out + ".buf").c_str(), "wb");
    const double cap = static_cast<double>(V_) * c_.depth;
    for (int w = 0; w < n_windows_; ++w) {
      const int span = std::min(c_.window, c_.cycles - c_.warmup - w * c_.window);
      for (int i = 0; i < n_ * kPorts; ++i) {
        const float v = static_cast<float>(buf_acc_[static_cast<size_t>(w) * n_ * kPorts + i] / (span * cap));
        std::fwrite(&v, sizeof(v), 1, fb);
      }
    }
    std::fclose(fb);
    std::FILE* fp = std::fopen((c_.out + ".pkt").c_str(), "wb");
    for (const Packet& pk : pkts_) {
      const int32_t row[11] = {pk.src, pk.dst, pk.t_create, pk.t_head_out, pk.t_eject, pk.hops, pk.vhops,
                               pk.dropped ? 1 : 0, pk.elevator, pk.len, pk.attacker ? 1 : 0};
      std::fwrite(row, sizeof(row), 1, fp);
    }
    std::fclose(fp);
    std::FILE* fa = std::fopen((c_.out + ".act").c_str(), "wb");
    std::fwrite(act_.data(), sizeof(long long), act_.size(), fa);
    std::fclose(fa);
    std::FILE* fs = std::fopen((c_.out + ".stats").c_str(), "w");
    std::fprintf(fs,
                 "{\"deadlock_cycle\": %d, \"drained\": %s, \"drain_cycles\": %d, \"stuck_flits\": %lld, "
                 "\"cycles_run\": %d, \"flit_hops_h\": %lld, \"flit_hops_v\": %lld, "
                 "\"q_updates\": %lld, \"max_va_wait\": %d, "
                 "\"va_wait_gt1000\": %lld, \"n_packets\": %zu, \"in_flight_end\": %lld}\n",
                 deadlock_cycle_, drained_ ? "true" : "false", drain_cycles_, stuck_flits_, cycles_run_,
                 flit_hops_h_, flit_hops_v_, q_updates_, max_va_wait_,
                 va_wait_gt1000_, pkts_.size(), in_flight_flits_);
    std::fclose(fs);
    std::FILE* ft = std::fopen((c_.out + ".topo").c_str(), "w");
    std::fprintf(ft, "%d %d %d %d %d\n", c_.kx, c_.ky, c_.kz, n_windows_, c_.window);
    for (int r = 0; r < n_; ++r)
      for (int p = 1; p < kPorts; ++p)
        if (neighbor(r, p) >= 0) std::fprintf(ft, "%d %d %d\n", r, neighbor(r, p), (p == kUp || p == kDown) ? 1 : 0);
    std::fclose(ft);
  }

  Config c_;
  Alg alg_;
  std::mt19937_64 rng_;   // routing decisions
  std::mt19937_64 trng_;  // traffic generation (identical across classes)
  std::uniform_real_distribution<double> uni_{0.0, 1.0};
  bool vpc_ = false;
  int n_ = 0, V_ = 2, n_windows_ = 0;
  int all_ = 3, c0_ = 1, c1_ = 2, esc_ = 1, adp_ = 2;
  std::vector<std::array<int, 2>> elev_xy_;
  std::vector<int> elev_at_;
  std::vector<uint8_t> seg_ok_;
  std::vector<int> nbr_, ham_, xs_, ys_, zs_;
  std::vector<InputVC> in_;
  std::vector<OutputVC> outvc_;
  std::vector<std::deque<int>> srcq_;
  std::vector<int> inj_flit_, inj_vc_, rr_in_, rr_out_, rot_, elev_load_;
  std::vector<Packet> pkts_;
  std::vector<double> buf_acc_;
  std::vector<Cand> cands_, tmp_;
  std::vector<Move> moves_;
  std::vector<std::vector<CreditMsg>> credit_pipe_;
  std::vector<uint64_t> pver_;
  std::vector<int> occ_;
  bool cacheable_ = false;
  std::FILE* trace_ = nullptr;
  std::vector<double> agg_, agg_next_, rth_, temp_;
  std::vector<float> q_;
  std::vector<long long> act_, act_period_;
  std::vector<std::vector<int>> adele_subset_;
  std::vector<std::vector<double>> adele_cost_;
  std::vector<int> adele_ptr_;
  std::vector<uint8_t> attacker_, on_;
  std::vector<int> left_;
  double p_on_ = 0.0, xm_on_ = 1.0, xm_off_ = 1.0;
  int hotspot_ = 0;
  int now_ = 0, last_move_ = 0, deadlock_cycle_ = -1, cycles_run_ = 0, drain_cycles_ = 0, max_va_wait_ = 0;
  bool drained_ = false;
  long long in_flight_flits_ = 0, flit_hops_h_ = 0, flit_hops_v_ = 0, q_updates_ = 0;
  long long stuck_flits_ = 0, va_wait_gt1000_ = 0;
};

}  // namespace

int main(int argc, char** argv) {
  Config c;
  for (int i = 1; i + 1 < argc; i += 2) {
    const std::string k = argv[i], v = argv[i + 1];
    if (k == "--topo") c.topo = v;
    else if (k == "--routing") c.routing = v;
    else if (k == "--traffic") c.traffic = v;
    else if (k == "--pir") c.pir = std::atof(v.c_str());
    else if (k == "--seed") c.seed = std::strtoull(v.c_str(), nullptr, 10);
    else if (k == "--cycles") c.cycles = std::atoi(v.c_str());
    else if (k == "--warmup") c.warmup = std::atoi(v.c_str());
    else if (k == "--window") c.window = std::atoi(v.c_str());
    else if (k == "--pv") c.pv = std::atof(v.c_str());
    else if (k == "--xi") c.xi = std::atof(v.c_str());
    else if (k == "--elevators") c.elevators = v;
    else if (k == "--out") c.out = v;
    else if (k == "--trace") c.trace = v != "0";
    else if (k == "--kx") c.kx = std::atoi(v.c_str());
    else if (k == "--ky") c.ky = std::atoi(v.c_str());
    else if (k == "--kz") c.kz = std::atoi(v.c_str());
    else if (k == "--flits") c.pkt_flits = std::atoi(v.c_str());
    else if (k == "--depth") c.depth = std::atoi(v.c_str());
    else if (k == "--vcs") c.vcs = std::atoi(v.c_str());
    else if (k == "--hop-delay") c.hop_delay = std::atoi(v.c_str());
    else if (k == "--credit-delay") c.credit_delay = std::atoi(v.c_str());
    else if (k == "--watchdog") c.watchdog = std::atoi(v.c_str());
    else if (k == "--drain") c.drain = v != "0";
    else if (k == "--drain-max") c.drain_max = std::atoi(v.c_str());
    else if (k == "--check") c.check = v != "0";
    else if (k == "--atomic") c.atomic = v != "0";
    else if (k == "--inject-escape") c.inject_escape = v != "0";
    else if (k == "--rand-sel") c.rand_sel = v != "0";
    else if (k == "--dump-stuck") c.dump_stuck = v != "0";
    else if (k == "--alloc-cache") c.alloc_cache = v != "0";
    else if (k == "--hotspot-frac") c.hotspot_frac = std::atof(v.c_str());
    else if (k == "--attackers") c.attackers = std::atoi(v.c_str());
    else if (k == "--attack-pir") c.attack_pir = std::atof(v.c_str());
    else if (k == "--burst-factor") c.burst_factor = std::atof(v.c_str());
    else if (k == "--aes-beta") c.aes_beta = std::atof(v.c_str());
    else if (k == "--q-eta") c.q_eta = std::atof(v.c_str());
    else if (k == "--adele-k") c.adele_k = std::atoi(v.c_str());
    else {
      std::fprintf(stderr, "unknown option %s\n", k.c_str());
      return 2;
    }
  }
  Simulator sim(c);
  if (c.check) sim.check_cdg();
  else sim.run();
  return 0;
}
