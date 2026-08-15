#ifndef EMC_FULL_BAND_SCATTERING_HPP
#define EMC_FULL_BAND_SCATTERING_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

#include <FullBand/emcNumericBandStructure.hpp>

/*! \brief Table-driven full-band scattering from a .matpkg package
 * (spec v0.1, /scattering/phaseA).
 *
 * Loads energy-resolved rate tables Gamma_nu(T, band, E) and provides:
 *  - total scatter rate Gamma(E) (linear interpolation, clamped - the spec
 *    forbids extrapolation),
 *  - self-scattering constant Gamma0 = safety * max(Gamma),
 *  - mechanism selection proportional to partial rates,
 *  - DOS-weighted final-state selection: k' is sampled on the constant-energy
 *    surface E' of the tetrahedral mesh (per-tet weights ~ tet DOS, uniform
 *    sampling on the iso-surface polygon inside the chosen tetrahedron).
 *
 * Phase-A tables are isotropic in the coupling; anisotropy enters only through
 * the band structure itself (spec 8.1). Phase-B (k-resolved) replaces the
 * final-state sampler, not this interface.
 */
template <class T> class emcFullBandScattering {
public:
  using Vec3 = std::array<T, 3>;

  struct Mechanism {
    std::string name;
    T deltaE;              // energy change per event [eV]
    std::vector<T> grid;   // energy grid [eV]
    std::vector<T> rates;  // [1/s] at selected temperature (single band, v0.1)
  };

  emcFullBandScattering(const std::string &packageFile,
                        const emcNumericBandStructure<T> &bandStructure,
                        std::size_t band, T temperature, T maxEnergy,
                        T binWidth = static_cast<T>(2e-3))
      : bs(bandStructure), bandIdx(band) {
    loadMechanisms(packageFile, temperature);
    buildEnergyBins(maxEnergy, binWidth);
  }

  /// total scatter rate at energy [eV] -> [1/s]; clamped at table ends
  T getTotalRate(T energy) const {
    T sum = 0;
    for (const auto &m : mechanisms)
      sum += interp(m, energy);
    return sum;
  }

  /// self-scattering constant (>= max total rate on the tables)
  T getGamma0(T safety = static_cast<T>(1.2)) const {
    T mx = 0;
    for (std::size_t i = 0; i < mechanisms[0].grid.size(); i++) {
      T s = 0;
      for (const auto &m : mechanisms)
        s += m.rates[std::min(i, m.rates.size() - 1)];
      mx = std::max(mx, s);
    }
    return safety * mx;
  }

  /// selects mechanism index proportional to partial rates at given energy
  template <class RNG> std::size_t selectMechanism(T energy, RNG &rng) const {
    std::uniform_real_distribution<T> U(0, getTotalRate(energy));
    T r = U(rng), acc = 0;
    for (std::size_t i = 0; i < mechanisms.size(); i++) {
      acc += interp(mechanisms[i], energy);
      if (r <= acc)
        return i;
    }
    return mechanisms.size() - 1;
  }

  T getDeltaE(std::size_t mech) const { return mechanisms[mech].deltaE; }
  std::size_t getNrMechanisms() const { return mechanisms.size(); }
  const std::string &getMechanismName(std::size_t i) const {
    return mechanisms[i].name;
  }

  /// samples k' (Cartesian [1/m]) on the constant-energy surface at energy
  /// [eV]; returns false if the energy is outside the binned band range
  template <class RNG>
  bool sampleFinalState(T energy, RNG &rng, Vec3 &kOut,
                        std::int64_t &tetHint) const {
    const std::int64_t bin = binOf(energy);
    if (bin < 0 || bins[bin].tets.empty())
      return false;
    const auto &B = bins[bin];
    // pick tet ~ weight
    std::uniform_real_distribution<T> U(0, B.totalWeight);
    T r = U(rng), acc = 0;
    std::size_t pick = B.tets.size() - 1;
    for (std::size_t i = 0; i < B.tets.size(); i++) {
      acc += B.weights[i];
      if (r <= acc) {
        pick = i;
        break;
      }
    }
    const std::int64_t t = B.tets[pick];
    // iso-surface polygon of level `energy` inside tet t (in frac coords)
    const auto &tet = bs.getTetrahedra()[t];
    const auto &pts = bs.getPoints();
    const auto &E = bs.getBandEnergies(bandIdx);
    std::array<Vec3, 4> v;
    std::array<T, 4> e;
    for (int i = 0; i < 4; i++) {
      v[i] = pts[tet[i]];
      e[i] = E[tet[i]];
    }
    Vec3 poly[4];
    int np = 0;
    static const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                    {1, 2}, {1, 3}, {2, 3}};
    for (const auto &ed : edges) {
      const T ea = e[ed[0]], eb = e[ed[1]];
      if ((ea - energy) * (eb - energy) < 0) {
        const T f = (energy - ea) / (eb - ea);
        for (int c = 0; c < 3; c++)
          poly[np][c] = v[ed[0]][c] + f * (v[ed[1]][c] - v[ed[0]][c]);
        if (++np == 4)
          break;
      }
    }
    if (np < 3) { // degenerate (energy grazes a vertex): fall back to centroid
      Vec3 c{};
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
          c[j] += v[i][j] / 4;
      kOut = bs.fracToCart(c);
      tetHint = t;
      return true;
    }
    if (np == 4) // order quad vertices to avoid bow-tie (swap if needed)
      orderQuad(poly);
    // triangulate fan (1 or 2 triangles), pick by area, sample uniformly
    const int ntri = np - 2;
    T areas[2] = {0, 0};
    for (int i = 0; i < ntri; i++)
      areas[i] = triArea(poly[0], poly[i + 1], poly[i + 2]);
    std::uniform_real_distribution<T> U01(0, 1);
    const int tri =
        (ntri == 2 && U01(rng) * (areas[0] + areas[1]) > areas[0]) ? 1 : 0;
    T a = U01(rng), b2 = U01(rng);
    if (a + b2 > 1) {
      a = 1 - a;
      b2 = 1 - b2;
    }
    Vec3 f{};
    for (int c = 0; c < 3; c++)
      f[c] = poly[0][c] + a * (poly[tri + 1][c] - poly[0][c]) +
             b2 * (poly[tri + 2][c] - poly[0][c]);
    kOut = bs.fracToCart(f);
    tetHint = t;
    return true;
  }

private:
  const emcNumericBandStructure<T> &bs;
  std::size_t bandIdx;
  std::vector<Mechanism> mechanisms;

  struct Bin {
    std::vector<std::int64_t> tets;
    std::vector<T> weights;
    T totalWeight = 0;
  };
  T binLo = 0, binW = 0;
  std::vector<Bin> bins;

  std::int64_t binOf(T energy) const {
    const std::int64_t b = static_cast<std::int64_t>((energy - binLo) / binW);
    return (b >= 0 && b < (std::int64_t)bins.size()) ? b : -1;
  }

  static T interp(const Mechanism &m, T energy) {
    if (energy <= m.grid.front())
      return m.rates.front(); // clamp (spec: no extrapolation)
    if (energy >= m.grid.back())
      return m.rates.back();
    const auto it = std::upper_bound(m.grid.begin(), m.grid.end(), energy) - 1;
    const std::size_t i = it - m.grid.begin();
    const T f = (energy - m.grid[i]) / (m.grid[i + 1] - m.grid[i]);
    return m.rates[i] * (1 - f) + m.rates[i + 1] * f;
  }

  static T triArea(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
    const T u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const T w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const T x = u[1] * w[2] - u[2] * w[1], y = u[2] * w[0] - u[0] * w[2],
            z = u[0] * w[1] - u[1] * w[0];
    return std::sqrt(x * x + y * y + z * z) / 2;
  }

  /// order 4 coplanar points into a simple quad (angle sort around centroid)
  static void orderQuad(Vec3 *p) {
    Vec3 c{};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
        c[j] += p[i][j] / 4;
    // plane basis from first two edges
    Vec3 u{p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
    Vec3 n{}; // normal via cross with another edge
    const Vec3 w{p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
    n = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
         u[0] * w[1] - u[1] * w[0]};
    Vec3 v2{n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2],
            n[0] * u[1] - n[1] * u[0]};
    std::array<std::pair<T, int>, 4> ang;
    for (int i = 0; i < 4; i++) {
      const Vec3 d{p[i][0] - c[0], p[i][1] - c[1], p[i][2] - c[2]};
      const T x = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
      const T y = d[0] * v2[0] + d[1] * v2[1] + d[2] * v2[2];
      ang[i] = {std::atan2(y, x), i};
    }
    std::sort(ang.begin(), ang.end());
    const Vec3 tmp[4] = {p[ang[0].second], p[ang[1].second], p[ang[2].second],
                         p[ang[3].second]};
    for (int i = 0; i < 4; i++)
      p[i] = tmp[i];
  }

  void loadMechanisms(const std::string &file, T temperature) {
    hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0)
      throw std::runtime_error("matpkg: cannot open " + file);
    hid_t grp = H5Gopen2(f, "/scattering/phaseA", H5P_DEFAULT);
    if (grp < 0)
      throw std::runtime_error("matpkg: no /scattering/phaseA group");
    H5G_info_t info;
    H5Gget_info(grp, &info);
    for (hsize_t i = 0; i < info.nlinks; i++) {
      char name[256];
      H5Lget_name_by_idx(grp, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name,
                         sizeof(name), H5P_DEFAULT);
      Mechanism m;
      m.name = name;
      hid_t mg = H5Gopen2(grp, name, H5P_DEFAULT);
      hid_t at = H5Aopen(mg, "delta_E", H5P_DEFAULT);
      double dE = 0;
      H5Aread(at, H5T_NATIVE_DOUBLE, &dE);
      H5Aclose(at);
      m.deltaE = static_cast<T>(dE);
      // temperature index: nearest entry on the temperature grid
      std::vector<double> tg = read1D(mg, "temperature_grid");
      std::size_t ti = 0;
      for (std::size_t j = 1; j < tg.size(); j++)
        if (std::abs(tg[j] - temperature) < std::abs(tg[ti] - temperature))
          ti = j;
      std::vector<double> eg = read1D(mg, "energy_grid");
      m.grid.assign(eg.begin(), eg.end());
      // rates shape (NT, Nb, NE) - v0.1 prototype consumes band slice 0
      hid_t ds = H5Dopen2(mg, "rates", H5P_DEFAULT);
      hid_t sp = H5Dget_space(ds);
      hsize_t dims[3];
      H5Sget_simple_extent_dims(sp, dims, nullptr);
      std::vector<double> buf(dims[0] * dims[1] * dims[2]);
      H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              buf.data());
      m.rates.resize(dims[2]);
      for (hsize_t e = 0; e < dims[2]; e++)
        m.rates[e] = static_cast<T>(buf[(ti * dims[1] + 0) * dims[2] + e]);
      H5Sclose(sp);
      H5Dclose(ds);
      H5Gclose(mg);
      mechanisms.push_back(std::move(m));
    }
    H5Gclose(grp);
    H5Fclose(f);
    if (mechanisms.empty())
      throw std::runtime_error("matpkg: no phaseA mechanisms found");
  }

  static std::vector<double> read1D(hid_t loc, const char *name) {
    hid_t ds = H5Dopen2(loc, name, H5P_DEFAULT);
    hid_t sp = H5Dget_space(ds);
    hsize_t n;
    H5Sget_simple_extent_dims(sp, &n, nullptr);
    std::vector<double> out(n);
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Sclose(sp);
    H5Dclose(ds);
    return out;
  }

  /// per-tet DOS-proportional weights per energy bin:
  /// weight = vol(tet) * overlap(bin, [emin,emax]) / (emax - emin)
  void buildEnergyBins(T maxEnergy, T width) {
    const auto &tets = bs.getTetrahedra();
    const auto &pts = bs.getPoints();
    const auto &E = bs.getBandEnergies(bandIdx);
    binLo = bs.getBandMinimum(bandIdx);
    binW = width;
    const std::size_t nb =
        static_cast<std::size_t>((maxEnergy - binLo) / width) + 1;
    bins.assign(nb, Bin{});
    for (std::int64_t t = 0; t < (std::int64_t)tets.size(); t++) {
      T lo = E[tets[t][0]], hi = lo;
      for (int i = 1; i < 4; i++) {
        lo = std::min(lo, E[tets[t][i]]);
        hi = std::max(hi, E[tets[t][i]]);
      }
      if (hi <= binLo || lo >= maxEnergy || hi - lo < 1e-12)
        continue;
      // frac-space volume as DOS proxy (constant Jacobian to cartesian)
      const auto &a = pts[tets[t][0]];
      T e1[3], e2[3], e3[3];
      for (int c = 0; c < 3; c++) {
        e1[c] = pts[tets[t][1]][c] - a[c];
        e2[c] = pts[tets[t][2]][c] - a[c];
        e3[c] = pts[tets[t][3]][c] - a[c];
      }
      const T vol = std::abs(
          e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) -
          e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) +
          e1[2] * (e2[0] * e3[1] - e2[1] * e3[0])) / 6;
      const std::int64_t b0 = std::max<std::int64_t>(0, binOf(lo));
      const std::int64_t b1 =
          std::min<std::int64_t>(bins.size() - 1, binOf(std::min(hi, maxEnergy)) < 0
                                                      ? bins.size() - 1
                                                      : binOf(std::min(hi, maxEnergy)));
      for (std::int64_t b = b0; b <= b1; b++) {
        const T bl = binLo + b * binW, bh = bl + binW;
        const T ov = std::min(hi, bh) - std::max(lo, bl);
        if (ov <= 0)
          continue;
        const T w = vol * ov / (hi - lo);
        bins[b].tets.push_back(t);
        bins[b].weights.push_back(w);
        bins[b].totalWeight += w;
      }
    }
  }
};

#endif // EMC_FULL_BAND_SCATTERING_HPP
