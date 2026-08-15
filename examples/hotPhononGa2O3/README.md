# Hot-Phonon Transport in β-Ga₂O₃

Simulates hot-phonon-limited high-field electron transport in bulk
β-Ga₂O₃. Published β-Ga₂O₃ Monte Carlo studies evaluate the polar-optical
scattering rate with the phonon population held at its equilibrium
Bose–Einstein value. β-Ga₂O₃ is the wide-bandgap material in which that
assumption is least safe: it has the strongest Fröhlich coupling of the
family, low optical-phonon energies (high emission rate), and a thermal
conductivity 5–25× below GaN's, so the acoustic reservoir that absorbs
the decayed LO energy is itself slow to unload.

**Executable:** `build/examples/hotPhononGa2O3/hotPhononGa2O3`

## What it models

The driver runs the same transport problem two ways and reports the
velocity–field curve, the mean energy, and the phonon populations for
both:

- `--use_hpb 0` — equilibrium phonon occupation (the published
  assumption).
- `--use_hpb 1` — wavevector-resolved non-equilibrium occupation N_q
  evolved self-consistently with the carrier ensemble
  (`emcPhononBath`), optionally coupled to a finite-lifetime acoustic
  reservoir (three-temperature model, `--use_acoustic 1`).

Optional mechanisms, all off by default:

- `--use_multimode 1` — the five-mode ab-initio polar-phonon set instead
  of the single 44 meV effective mode.
- `--use_screening 1` — free-carrier (plasmon) screening of the Fröhlich
  vertex. `--screen_eps` sets the background permittivity entering the
  Debye wavevector; the static constant `epsLo` is the default, the
  optical constant `epsHi` the physically motivated alternative at the
  LO frequency. Both choices are deliberate modelling options.
- `--use_impurity 1` — ionized-impurity (Brooks–Herring) scattering.
- `--use_qresolved 1` — draw the phonon occupation entering the rate from
  the coupling-weighted mean over the accessible wavevector window;
  `--qres_angle` (default on) samples the final-state angle from the
  same distribution.
- `--use_pauli 1` — Pauli exclusion on a k-space occupancy grid.

Note that the LO and acoustic lifetimes are estimates; no measured
values exist for β-Ga₂O₃, so results should be swept over the stated
ranges.

## Running the example

```bash
cd build/examples/hotPhononGa2O3

# Velocity-field curve, equilibrium vs hot phonons
./hotPhononGa2O3 --fields 10,50,100,200,300,400 --use_hpb 0
./hotPhononGa2O3 --fields 10,50,100,200,300,400 --use_hpb 1 --tau_lo 5e-12
```

## Key command-line flags

| Flag | Default | Description |
|---|---|---|
| `--fields` | — | Comma-separated field list [kV/cm] |
| `--use_hpb` | `1` | Non-equilibrium (resolved) phonon occupation |
| `--use_acoustic` | `1` | Acoustic reservoir (three-temperature model) |
| `--tau_lo` | `5e-12` | LO lifetime [s]; estimated, sweep 0.5–20 ps |
| `--tau_ac` | `20e-12` | Acoustic lifetime [s]; estimated |
| `--temp` | `300` | Lattice temperature [K] |
| `--doping` | `1e23` | Doping density [m⁻³] (10¹⁷ cm⁻³) |
| `--box` | `3e-7` | Cubic box side [m] |
| `--dt` / `--time` | `1e-16` / `5e-12` | Time step / total time [s] |
| `--steady_frac` | `0.4` | Final fraction of the run averaged for the steady state |
| `--reinit_every` | `1` | Steps between scatter-table rebuilds; N_q evolves on the τ_LO scale, so a modest stride is harmless but must be converged against 1 |
| `--use_multimode` | `0` | Five-mode ab-initio polar set |
| `--use_screening` / `--screen_eps` | `0` / `epsLo` | Plasmon screening and its background permittivity |
| `--use_impurity` | `0` | Ionized-impurity scattering |
| `--use_qresolved` / `--qres_angle` | `0` / `1` | Coupling-weighted N_q in the rate / final-angle sampling |
| `--use_pauli` | `0` | Pauli exclusion (`--pauli_dk`, `--pauli_emax`) |
| `--seed` / `--threads` | `1` / `4` | RNG seed / OpenMP threads |
| `--outdir` | `.` | Output directory |

## Output files

Per run, into `--outdir` (`<tag>` encodes the enabled mechanisms):

| File | Contents |
|---|---|
| `ga2o3_vE_<tag>.txt` | field [kV/cm], drift velocity [cm/s], ⟨E⟩ [eV], N_LO, T_LO, T_ac |
| `ga2o3_trace_<tag>_F<field>.txt` | time trace of the same quantities at one field |
| `ga2o3_nq_<tag>_F<field>.txt` | per-bin phonon spectrum N_q(q) |
