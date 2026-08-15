# Hot-Carrier Dynamics in Metal-Halide Perovskites

Simulates photo-excited hot-carrier cooling and device operation in a bulk
metal-halide perovskite (MHP). It extends the single-carrier
hot-phonon-bottleneck model of Faber, Filipovic, Koster, *J. Phys. Chem.
Lett.* **2024**, 15, 12601 into a bipolar, three-temperature device
simulation with a wavevector-resolved phonon population.

**Executable:** `build/examples/hotCarrierMHP/hotCarrierMHP`

## What it models

- **Hot-phonon bottleneck (HPB)** — carriers emit LO phonons (Fröhlich);
  the LO occupation rises above its equilibrium value and stimulated
  re-absorption slows further cooling. `emcPhononBath` resolves the
  occupation N_q per wavevector bin and updates each bin every step as

  ```
  Nq[i] += (nEmitted[i] - nAbsorbed[i]) / Dph  -  (dt / tauLO) * (Nq[i] - Nstar)
  ```

  where `Dph = q² dq Vsim / (2π²)` is the number of phonon modes in the bin
  and `tauLO` the Klemens lifetime. The Fröhlich rates are rebuilt every
  step from the current N_q. With `--use_qres 1` (default) both the rate
  and the final-state angle are drawn from the coupling-weighted occupation
  over the accessible wavevector window rather than the
  density-of-states mean.
- **Free-carrier screening** — with `--use_screening 1` (default) the
  Fröhlich vertex is screened, 1/q² → 1/(q² + q_s²), with the Debye
  wavevector q_s refreshed each step from the live density and carrier
  temperature of both species.
- **Acoustic cascade** — a three-temperature carrier/LO/acoustic model in
  which the LO decay target follows a finite-lifetime acoustic reservoir
  (`--tau_ac`) that the LO modes back-heat. An optional **Ridley** channel
  (LO → LA + TO) competes with Klemens decay.
- **Non-parabolic (Kane) bands** — per-species two-band Kane parameters
  α = (1/E_gap)(1 − m*/m₀)² are applied by default; override with
  `--alpha`, `--alpha_e`, `--alpha_h` (0 recovers parabolic bands).
- **Bipolar transport** — electrons and holes share the phonon bath, with
  intra-species carrier–carrier and inter-species electron–hole scattering.
- **Recombination** — radiative (Bnp) and Auger (CCCH/CHHS); Auger returns
  the gap energy to a surviving carrier (Auger reheating).
- **Band filling** — Lugli–Ferry Pauli blocking on a k-space occupancy grid.
- **Energy-selective contact (ESC)** — extraction through a window
  [E_ex ± ΔE/2], from which the working open-circuit voltage and PCE are
  computed.

Optionally, `--use_multimode 1` resolves the two IR-active MAPbI₃ LO
branches instead of a single effective mode.

## Default model and the averaged limit

The defaults run the resolved model: 40 phonon bins of width 5×10⁷ m⁻¹,
screened and q-resolved coupling, Kane bands, and a 150 nm box. The
wavevector-averaged model of the previously published runs is recovered
with

```bash
./hotCarrierMHP --nr_phonon_bins 1 --dq_bin 2e9 --use_qres 0 \
                --use_screening 0 --alpha 0
```

Unscreened q-resolved runs do not settle, because the 1/q coupling weight
diverges at the small wavevectors where the emitted excess collects;
combine `--use_qres 1` with `--use_screening 1`.

## Running the example

Every mechanism and material parameter is a command-line flag — there are
no compile-time switches. Run from the build folder:

```bash
cd build/examples/hotCarrierMHP

# Full default model (MAPbI3, all mechanisms on)
./hotCarrierMHP

# Equilibrium-phonon reference (HPB off)
./hotCarrierMHP --use_hot_phonon 0

# CsSnI3 presets with a selective device contact
./hotCarrierMHP --material sn --E_ex 0.2 --delta_E 0.05

# Acoustic-lifetime study, reproducible
./hotCarrierMHP --tau_ac 10e-12 --box 200e-9 --seed 1

# Dump the per-bin phonon spectrum N_q(q,t)
./hotCarrierMHP --dump_nq 1
```

At the end of a run the console prints the working `V_OC`, extraction
`yield`, and `PCE`.

## Key command-line flags

Defaults are for MAPbI₃; `--material` selects presets (`pb`, `sn`,
`fapbi3`, `cspbi3`, `cspbbr3`). Times are in seconds, energies in eV,
densities in m⁻³. Boolean flags take `0`/`1`.

| Flag | Default | Description |
|---|---|---|
| `--material` | `pb` | Presets: `pb` MAPbI₃, `sn` CsSnI₃, `fapbi3`, `cspbi3`, `cspbbr3` |
| `--use_hot_phonon` | `1` | Hot-phonon bottleneck (dynamic N_q) |
| `--use_acoustic` | `1` | Acoustic reservoir (three-temperature model) |
| `--use_holes` | `1` | Bipolar transport (electrons + holes) |
| `--use_recomb` | `1` | Radiative + Auger recombination |
| `--use_cc` | `1` | Carrier–carrier scattering |
| `--use_esc` | `1` | Energy-selective extraction |
| `--use_bf` | `0` | Pauli band filling (selects the large box) |
| `--use_screening` | `1` | Debye screening of the Fröhlich vertex |
| `--use_qres` | `1` | Coupling-weighted N_q for rate and angle |
| `--nr_phonon_bins` | `40` | Phonon bins; `1` = averaged (single-mode) limit |
| `--dq_bin` | `5e7` | Bin width [1/m] |
| `--alpha` / `--alpha_e` / `--alpha_h` | Kane | Non-parabolicity [1/eV]; per-species two-band values unless overridden, `0` = parabolic |
| `--use_multimode` | `0` | Two-branch (multi-mode) Fröhlich coupling |
| `--ridley_w` | `0` | Ridley LO→LA+TO branching fraction (0–1) |
| `--tau_lo` | `0.6e-12` (pb) / `0.1e-12` (sn) | LO (Klemens) lifetime; MAPbI₃ default is the measured value of Fu et al., *Nat. Commun.* **8**, 1300 (2017) |
| `--tau_lo_f0` / `--tau_lo_f1` | `1` / `1` | Linear τ_LO(q) profile endpoints (fractions of `--tau_lo`) |
| `--tau_ac` | `30e-12` | Acoustic-reservoir lifetime |
| `--screen_eps` | `eps_hi` | Background permittivity in the Debye wavevector |
| `--E_photon` | `3.1` | Pump photon energy |
| `--E_ex` / `--delta_E` | `0.4` / `0.05` | ESC window centre / width |
| `--density` | `1e24` | Carrier density (10²⁴ m⁻³ = 10¹⁸ cm⁻³) |
| `--box` | auto | Cubic box side [m]; `0` auto-selects (150 nm, or 464 nm with `--use_bf`) |
| `--total_time` | `10e-12` | Simulation time |
| `--dump_nq` | `0` | Write the per-bin spectrum `phononSpectrum*.txt` |
| `--seed` | `0` | RNG seed; `0` = wall-clock, nonzero = reproducible |

Further overrides: `--eps_hi`, `--eps_lo`, `--mass_e`, `--mass_h`,
`--omega_lo`, `--omega_to`, `--tau_to`, `--tau_ex`, `--E_gap`, `--T_lat`,
`--E_ex_e`, `--E_ex_h`.

## Output files

Each filename carries a suffix encoding the enabled mechanisms — `HPB` or
`EQ`, then `_AC` (acoustic), `_MM` (multi-mode), `_RID` (Ridley), `_BF`
(band filling), `_1C` (single-carrier), `_ESC` (extraction). The full
default model, for example, writes `carrierTempHPB_AC_ESC.txt`.

| File (before suffix) | Columns | Description |
|---|---|---|
| `avgEnergyElectrons` / `avgEnergyHoles` | `t  <E>[eV]` | Mean carrier energy vs time |
| `phononOccupation` | `t  N_LO  [N_ac  T_ac]  [N_TO  T_TO]  [q_s]` | LO occupation plus reservoir state |
| `nrCarriers` | `t  N_e  N_h  N_esc_e  N_esc_h` | Surviving and extracted carrier counts |
| `carrierTemp` | `t  T_MB_e  T_MB_h  T_FD_e  μ_e  T_FD_h  μ_h` | Maxwell–Boltzmann and Fermi–Dirac carrier temperatures / chemical potentials |
| `phononSpectrum` (with `--dump_nq 1`) | blocks of `q  N_q` per dump time | Per-bin LO occupation spectrum |

> **Reproducibility:** particle moves run over OpenMP threads, each with
> its own RNG stream, so the seed **and** the thread count together define
> a trajectory. Pin `OMP_NUM_THREADS` and pass a nonzero `--seed` for
> repeatable runs.
