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

#include <cstdio>
#include <random>
#include <vector>

int main(int argc, char **argv) {
  using NBS = emcNumericBandStructure<double>;
  using FBS = emcFullBandScattering<double>;
  constexpr double HBAR = 1.054571817e-34, QE = 1.602176634e-19,
                   KB = 8.617333262e-5; // eV/K
  const std::string pkg =
      argc > 1 ? argv[1] : "/home/filipov/Software/FullBandMC/matforge/Si.matpkg";
  const double T = 300.0;               // K
  const int nPart = 2000;
  const double tTotal = 20e-12, tTransient = 5e-12;

  NBS bs(pkg);
  const double cbm = bs.getBandMinimum(0);
  FBS scat(pkg, bs, 0, T, cbm + 1.0);
  const double gamma0 = scat.getGamma0();
  std::printf("# %s: CBM %.4f eV, Gamma0 = %.3e 1/s, %zu mechanism(s)\n",
              pkg.c_str(), cbm, gamma0, scat.getNrMechanisms());
  std::printf("# %8s %14s %14s %14s %12s\n", "F[V/cm]", "vd_x[m/s]",
              "mu[cm2/Vs]", "<E>-CBM[eV]", "selfscat[%]");

  for (const double Fcm : {500.0, 1000.0, 2000.0, 5000.0}) {
    const double F = Fcm * 100.0; // V/m, along x
    std::mt19937_64 rng(1234);
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
    for (auto &p : ps) {
      double t = 0;
      while (t < tTotal) {
        const double tau = -std::log(1 - U01(rng)) / gamma0;
        const double dt = std::min(tau, tTotal - t);
        // measure with pre-flight velocity (piecewise-constant per tet)
        const auto v = bs.getVelocity(p.k, 0, p.hint);
        // free flight: dk = -e F dt / hbar (electron)
        p.k[0] -= QE * F * dt / HBAR;
        if (t + dt > tTransient) {
          const double w = std::min(dt, t + dt - tTransient);
          sumVdt += v[0] * w;
          sumEdt += (bs.getEnergy(p.k, 0, p.hint) - cbm) * w;
          sumT += w;
        }
        t += dt;
        if (dt < tau)
          break; // reached tTotal mid-flight
        const double E = bs.getEnergy(p.k, 0, p.hint);
        if (U01(rng) * gamma0 < scat.getTotalRate(E)) {
          nReal++;
          const std::size_t m = scat.selectMechanism(E, rng);
          const double Ef = E + scat.getDeltaE(m);
          NBS::Vec3 kNew;
          std::int64_t hNew = p.hint;
          if (scat.sampleFinalState(Ef, rng, kNew, hNew)) {
            p.k = kNew;
            p.hint = hNew;
          } // else: energy left the table range - keep state (counts as self)
        } else {
          nSelf++;
        }
      }
    }
    const double vd = sumVdt / sumT;
    const double mu = -vd / F * 1e4; // electrons drift against F; cm^2/Vs
    std::printf("  %8.0f %14.4e %14.1f %14.4f %12.1f\n", Fcm, vd, mu,
                sumEdt / sumT, 100.0 * nSelf / std::max(1L, nSelf + nReal));
  }
  return 0;
}
