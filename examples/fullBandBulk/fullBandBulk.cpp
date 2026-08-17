/*! fullBandBulk - full-band ensemble MC for bulk transport (tasks 2.2-2.4).
 *
 * Reads a .matpkg material package (bands on a tetrahedral BZ mesh + Phase-A
 * scattering tables) and runs an ensemble of electrons under a homogeneous
 * field: event-driven free flight with self-scattering (constant Gamma0),
 * k-update dk/dt = -eF/hbar, table-driven scattering with DOS-weighted
 * final-state selection on constant-energy surfaces.
 *
 * Output: time-averaged drift velocity, mobility and mean energy per field.
 */
#include <FullBand/emcFullBandScattering.hpp>
#include <FullBand/emcNumericBandStructure.hpp>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <thread>
#include <vector>

int main(int argc, char **argv) {
  using NBS = emcNumericBandStructure<double>;
  using FBS = emcFullBandScattering<double>;
  constexpr double HBAR = 1.054571817e-34, QE = 1.602176634e-19,
                   KB = 8.617333262e-5; // eV/K
  const std::string pkg =
      argc > 1 ? argv[1] : "/home/filipov/Software/FullBandMC/matforge/Si.matpkg";
  const double T = 300.0;               // K
  // argv[5] = particles per replica, argv[6] = independent replicas.
  // Replicas give an honest standard error: single-run mu noise is +-5-6%
  // at the default 2000 x 20 ps, which hides every few-percent effect.
  const int nPart = argc > 5 ? std::stoi(argv[5]) : 2000;
  const int nRep = argc > 6 ? std::stoi(argv[6]) : 1;
  const double tTotal = 20e-12, tTransient = 5e-12;

  NBS bs(pkg);
  // QUADE=0 disables the quadratic energy interpolation (A/B testing)
  if (const char *q = std::getenv("QUADE"))
    bs.setQuadraticEnergy(std::atoi(q) != 0);
  const double cbm = bs.getBandMinimum(0);
  const int mode = argc > 3 ? std::stoi(argv[3]) : 0;  // 0 auto, 1 phaseA, 2 no-g2
  FBS scat(pkg, bs, 0, T, cbm + 1.0, 2e-3, mode);
  const double gamma0 = scat.getGamma0();
  std::printf("# %s: CBM %.4f eV, Gamma0 = %.3e 1/s, %zu mechanism(s), %s\n",
              pkg.c_str(), cbm, gamma0, scat.getNrMechanisms(),
              scat.hasG2Tables() ? "k-RESOLVED + g2 FINAL STATES (Stage 2)"
              : (scat.isKResolved() ? "k-RESOLVED (phaseB)" : "energy tables (phaseA)"));
  std::printf("# ensemble: %d particles x %d replica(s) x %.0f ps\n",
              nPart, nRep, tTotal * 1e12);
  std::printf("# %8s %14s %14s %10s %14s %12s %10s\n", "F[V/cm]", "vd[m/s]",
              "mu[cm2/Vs]", "SE(mu)", "<E>-CBM[eV]", "selfscat[%]",
              "<cos_th>", "fail[%]", "accel");

  std::vector<double> fields{500.0, 1000.0, 2000.0, 5000.0};
  if (argc > 2) {                      // comma-separated F list [V/cm]
    fields.clear();
    std::string s(argv[2]);
    size_t p = 0;
    while (p < s.size()) {
      size_t c = s.find(',', p);
      if (c == std::string::npos) c = s.size();
      fields.push_back(std::stod(s.substr(p, c - p)));
      p = c + 1;
    }
  }
  // field direction (argv[4]): "100" (default), "110", "111"
  double dir[3] = {1.0, 0.0, 0.0};
  if (argc > 4) {
    const std::string d(argv[4]);
    for (int i = 0; i < 3 && i < (int)d.size(); i++)
      dir[i] = d[i] - '0';
    const double n = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    for (int i = 0; i < 3; i++) dir[i] /= n;
    std::printf("# field along <%s>\n", d.c_str());
  }

  // one independent replica -> (mu, <E>, selfscat%)
  struct Res { double mu, energy, self, cosTheta, failFrac, accel; };
  auto runReplica = [&](int rep, double F) {
    std::mt19937_64 rng(1234ull + 7919ull * static_cast<unsigned>(rep));
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    // thermal initialization: E ~ Boltzmann above CBM, k via iso-surface
    // sampler (also exercises the final-state machinery)
    struct P {
      NBS::Vec3 k;
      std::int64_t hint = -1;
    };
    std::vector<P> ps(nPart);
    for (auto &p : ps) {
      double E;
      do
        E = cbm - KB * T * std::log(1 - U01(rng));
      while (!scat.sampleFinalState(E, rng, p.k, p.hint));
    }

    double sumVdt = 0, sumEdt = 0, sumT = 0;
    long nSelf = 0, nReal = 0;
    double cosAcc = 0;      // <cos(theta)> momentum-randomization diagnostic
    long cosN = 0;
    long nFail = 0;         // real events whose final state could not be sampled
    // effective acceleration between real collisions: must equal eF/m*.
    // The ballistic test only verifies this over long flights (many cells);
    // here the field moves k by ~2% of a cell between collisions, which is
    // the regime that actually sets the mobility.
    double dvAcc = 0, dtAcc = 0;
    for (auto &p : ps) {
      double t = 0;
      double tLast = 0, vxLast = 0;
      bool haveLast = false;
      // trapezoidal drift accumulation: sampling only the pre-flight velocity
      // loses the in-flight acceleration, a bias of order Gamma/Gamma0 (it
      // destroyed a constant-rate test package where Gamma0 = 1.2 Gamma, and
      // costs ~1% where self-scattering dominates). The end-of-flight velocity
      // is the next flight's start velocity, so the trapezoid is free.
      auto vCur = bs.getVelocity(p.k, 0, p.hint);
      while (t < tTotal) {
        const double tau = -std::log(1 - U01(rng)) / gamma0;
        const double dt = std::min(tau, tTotal - t);
        // free flight: dk = -e F dt / hbar (electron)
        for (int i = 0; i < 3; i++)
          p.k[i] -= QE * F * dir[i] * dt / HBAR;
        const auto vEnd = bs.getVelocity(p.k, 0, p.hint);
        if (t + dt > tTransient) {
          const double w = std::min(dt, t + dt - tTransient);
          const double vMid = 0.5 * ((vCur[0] + vEnd[0]) * dir[0] +
                                     (vCur[1] + vEnd[1]) * dir[1] +
                                     (vCur[2] + vEnd[2]) * dir[2]);
          sumVdt += vMid * w;
          sumEdt += (bs.getEnergy(p.k, 0, p.hint) - cbm) * w;
          sumT += w;
        }
        vCur = vEnd;
        t += dt;
        if (dt < tau)
          break; // reached tTotal mid-flight
        const double E = bs.getEnergy(p.k, 0, p.hint);
        if (U01(rng) * gamma0 < scat.getTotalRate(p.k, E, p.hint)) {
          nReal++;
          const std::size_t m = scat.selectMechanism(p.k, E, p.hint, rng);
          const double Ef = E + scat.getDeltaE(m);
          NBS::Vec3 kNew;
          std::int64_t hNew = p.hint;
          const auto vPre = bs.getVelocity(p.k, 0, p.hint);
          if (scat.sampleFinalStateG2(m, p.k, Ef, rng, kNew, hNew)) {
            // momentum-randomization diagnostic: <cos(theta)> between the
            // velocity before and after a real scattering event. Isotropic
            // final states give 0; a positive value means the event only
            // partially destroys momentum, so the transport (momentum)
            // relaxation time is tau/(1-<cos>) - the quantity that sets
            // mobility, not the total scattering time.
            const auto vPost = bs.getVelocity(kNew, 0, hNew);
            const double n1 = std::sqrt(vPre[0]*vPre[0] + vPre[1]*vPre[1] +
                                        vPre[2]*vPre[2]);
            const double n2 = std::sqrt(vPost[0]*vPost[0] + vPost[1]*vPost[1] +
                                        vPost[2]*vPost[2]);
            if (n1 > 0 && n2 > 0) {
              cosAcc += (vPre[0]*vPost[0] + vPre[1]*vPost[1] +
                         vPre[2]*vPost[2]) / (n1 * n2);
              cosN++;
            }
            if (haveLast && t > tLast) {
              const double vxPre = vPre[0]*dir[0] + vPre[1]*dir[1] +
                                   vPre[2]*dir[2];
              dvAcc += vxPre - vxLast;
              dtAcc += t - tLast;
            }
            tLast = t;
            vxLast = vPost[0]*dir[0] + vPost[1]*dir[1] + vPost[2]*dir[2];
            haveLast = true;
            p.k = kNew;
            p.hint = hNew;
            vCur = vPost;      // discontinuous change at a real event
          } else {
            nFail++;   // no final state: momentum NOT randomized this event
          }
        } else {
          nSelf++;
        }
      }
    }
    const double vd = sumVdt / sumT;
    return Res{-vd / F * 1e4,        // electrons drift against F; cm^2/Vs
               sumEdt / sumT,
               100.0 * nSelf / std::max(1L, nSelf + nReal),
               cosN ? cosAcc / cosN : 0.0,
               100.0 * nFail / std::max(1L, nReal),
               dtAcc > 0 ? dvAcc / dtAcc : 0.0};
  };

  const unsigned nThreads = std::min<unsigned>(
      static_cast<unsigned>(nRep),
      std::max(1u, std::thread::hardware_concurrency()));

  for (const double Fcm : fields) {
    const double F = Fcm * 100.0; // V/m, along dir
    std::vector<Res> res(nRep);
    std::atomic<int> next{0};
    auto worker = [&]() {
      int r;
      while ((r = next.fetch_add(1)) < nRep)
        res[r] = runReplica(r, F);
    };
    std::vector<std::thread> pool;
    for (unsigned i = 0; i < nThreads; i++)
      pool.emplace_back(worker);
    for (auto &th : pool)
      th.join();

    double muM = 0, eM = 0, sM = 0, cM = 0, fM = 0, aM = 0;
    for (const auto &r : res) {
      muM += r.mu; eM += r.energy; sM += r.self; cM += r.cosTheta;
      fM += r.failFrac; aM += r.accel;
    }
    muM /= nRep; eM /= nRep; sM /= nRep; cM /= nRep; fM /= nRep; aM /= nRep;
    double se = 0;
    if (nRep > 1) {
      double var = 0;
      for (const auto &r : res) var += (r.mu - muM) * (r.mu - muM);
      se = std::sqrt(var / (nRep - 1) / nRep);   // standard error of the mean
    }
    std::printf("  %8.0f %14.4e %14.1f %10.1f %14.4f %12.1f %10.3f %9.1f"
                " %11.3e\n",
                Fcm, -muM * F * 1e-4, muM, se, eM, sM, cM, fM, aM);
  }
  return 0;
}
