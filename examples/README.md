# Examples: the calculations of the paper

Each example lists its GPT setup, the MAGNETO calculation, and the scripts that collect the data and draw the paper figure. Paths are relative to the repository root.

General workflow:

1. **GPT**: run `run.bat` in the GPT folder. This needs `GPTBIN` and `GPTLICENSE`, and the field maps in `field_map/` (see [`field_map/README.md`](../field_map/README.md)).
2. **MAGNETO**: run the `run_*.m` script in `matlab/...`. It writes `results/*.mat`.
3. **Paper data and figures**: run the `collect_*.m` script and then the `make_*.m` script in `paper/...`. They write `data/*.txt` and `figures/`.

Steps 1 and 2 produce large files that are not in the repository. The processed results of step 3 (`paper/*/data`, `paper/*/figures`) are included.

All MATLAB scripts run with `matlab -batch "<script name>"` from their own folder, or from any folder once it is on the path.

---

## 01 Single-electron undulator (paper Sec. III.A, Fig. 1)

240 MeV electron, planar undulator with 15 periods (λu = 30.47 mm, K = 0.98, RADIA field map). On-axis spectrum at (0, 0, 11 m). It is compared with SRW and with the analytical spectrum of an ideal undulator (Wiedemann, Eq. 26.46).

| Step | Location |
|---|---|
| GPT, field map | `gpt_sims/single_electron/E240MeV_und15p_paraxial/` (map: `field_map/und15p_K1_paraxial.gdf`) |
| GPT, ideal field (`unduplan`) | `gpt_sims/single_electron/E240MeV_unduplan_ideal/` |
| MAGNETO | `matlab/radiation_und15p_240MeV_single_electron/run_und240_spectra.m` |
| Analytical spectrum | `matlab/undulator_analytic_comparison/undulator240/undulator_analytic.m`; comparison `compare_undulator.m`; check `validate_analytic.m` |
| Ideal hard-edge trajectory (the trajectory assumed by the analytical model) | `matlab/undulator_analytic_comparison/undulator240/diagnostic_textbook_trajectory/diag_textbook.m` |
| SRW | `srw/undulator_240MeV_six_spectra_15_150nm/` |
| Paper data and figure | `paper/sec3_single_electron/make_sec3_single_electron.m`, `analytic/collect_sec3_analytic.m`, `make_sec3_figures_with_analytic.m` → `figures/fig_undulator_theta0_logy` |

Results: MAGNETO and SRW agree to within 0.06% of the local spectral envelope. With the ideal trajectory, MAGNETO and the analytical spectrum agree to within 1.6×10⁻³ at 11 m and 2.7×10⁻⁴ at 1 km.

## 02 Single-electron dipole (Sec. III.B, Fig. 2)

240 MeV electron, dipole of 0.1 m and 0.65 T with Enge fringe fields (1 cm gap), bending angle 81 mrad, λc = 49.6 nm. Two observers: on the incoming axis, at (0, 0, 11 m), and on the tangent at the arc midpoint, at (0, 0.407 m, 11 m).

| Step | Location |
|---|---|
| GPT | `gpt_sims/single_electron/E240MeV_dipole_L10cm/` (hard-edge variant: `E240MeV_dipole_L10cm_hardedge/`) |
| MAGNETO, on axis | `matlab/radiation_dipole240_single_electron/run_dipole240_onaxis.m` (step convergence: `tout_sweep.m`) |
| Analytical spectrum (Jackson, Eq. 14.79) | `matlab/dipole_analytic_comparison/dipole240_tangent/compare_dipole_tangent.m`, `dipole240_onaxis/compare_dipole_onaxis.m` |
| SRW | `srw/dipole_240MeV_onaxis_15_150nm/` (tangent observer: `dipole240_tangent/`) |
| Paper data and figure | `paper/sec3_single_electron/collect_dipole_trajectory.m`, `make_sec3_figures_with_analytic.m` → `figures/fig_dipole` |

Results: MAGNETO and SRW agree to within 0.35% (on axis). Along the arc-midpoint tangent, MAGNETO and the analytical spectrum agree to within 1.4×10⁻⁴.

## 03 Undulator and dipole together: relative phase (Sec. III.C, Fig. 3)

The same undulator, followed by the same dipole at a distance D (bending plane perpendicular to the wiggle plane). A single electron is tracked through both magnets for 213 values of D. The field is projected onto ±45° linear polarizations.

| Step | Location |
|---|---|
| GPT | `gpt_sims/single_electron/E240MeV_und15p_dipole_interference/` (`combined.in`, `dipole_only.in`, `undulator_only.in`) |
| GPT runs and MAGNETO for all D | `matlab/interference_und_dipole/run_scan.m` (calls GPT itself); checks `check_D0p5.m`; analysis `analyze_scan.m` |
| Paper data and figures | `paper/sec3_single_electron/interference/collect_sec5_interference.m`, `collect_sec5_fringe_vs_D.m`, `make_sec5_figures.m` → `figures/fig_interference_spectra_vs_D` |

Results: at 70 nm, the fringe period in D is 31.0 mm (2γ²λ) and the visibility is 0.98. The combined calculation agrees with the prediction from the separately calculated fields and the geometric delay to within 1.2×10⁻³ of the interference amplitude.

## 04 Coherent radiation of a 2000-electron bunch (Sec. IV, Figs. 4–6)

2000 electrons (one macro-particle per electron), 240 MeV, σx = σy = 20 µm, σz = 10 nm, no space charge. Observer on axis at 11 m.

| Step | Location |
|---|---|
| Initial distributions and GPT | `gpt_sims/multi_electron/E240MeV_und15p_N2000/` and `E240MeV_dipole_N2000/` (`make_distribution.m`, `run.bat`) |
| MAGNETO and form-factor comparison | `matlab/radiation_multi_electron/und15p_N2000_onaxis/` (`run_multi_onaxis.m`, `compare_form_factor.m`, `compare_translation.m`), `dipole_N2000_onaxis/` (`run_multi_onaxis.m`, `compare_three.m`, `check_shot_noise.m`) |
| Energy spread and emittance | GPT cases (map `field_map/und15p_K1_laplace_2mm.gdf`, GitHub Release asset): `gpt_sims/multi_electron/E240MeV_und15p_N2000_spread_emittance/` (`make_distributions.m`, one folder per case). MAGNETO: `matlab/radiation_multi_electron/und15p_N2000_spread_emittance/run_all.m`, `analyze_cases.m` |
| Paper data and figures | `paper/sec4_coherent/collect_sec4_coherent.m`, `collect_sec4_spread_emittance.m`, `make_fig4_1.m`, `make_fig4_2.m`, `make_fig4_3_lines_n1.m` |

Results: at the undulator fundamental (102 nm), S_b/(N² S₁) = 0.69, compared with |F|² = 0.68 for a Gaussian. The third and fifth harmonics are partially coherent and nearly incoherent, respectively. Energy spread and emittance reduce the coherent enhancement through the evolution of the bunch. A form factor of the initial distribution does not capture this.

## 05 Bunch-length diagnostic at FACET-II (Sec. V, Figs. 7–8)

7 GeV, 100 pC (6.2×10⁸ electrons, represented by 5×10⁴ macro-particles), 2% energy spread, normalized emittance 30 µm. The undulator is a 3-period Halbach array (NdFeB, 64 mm period, 6 mm gap, K ≈ 8). It is followed by a 1 m, 0.13 T dipole that bends the beam up by 5.6 mrad. The observer is at (0, −2 mm, 5.9 m) behind a horizontal polarizer. The coordinate origin is the bunch position at t = 0, 0.4 m upstream of the undulator center, so the energies per sr refer to R0 = 5.9 m. Detector bands: D1 50–280 nm, D2 300–1000 nm.

| Step | Location |
|---|---|
| Undulator field map | `field_map/gap6_laplace_0p5mm.gdf` (GitHub Release asset) |
| GPT (beam in 25 chunks of 2000) | `gpt_sims/multi_electron/E7GeV_und_dipole/` (`multi_electron.in`, `run.bat`); single reference electron: `ideal_electron/` |
| Beam sample | `matlab/radiation_multi_electron/E7GeV_gap6und_sigz_sweep/results/standard_sample.mat` (standard normal numbers, made by `run_sweep.m`) |
| MAGNETO, calibration (separate GPT run for each σz) | `matlab/radiation_multi_electron/E7GeV_und_dipole/run_direct_sigz.m`, `analyze_direct_sigz.m`, `check_direct_sigz.m` |
| MAGNETO, σz = 10 nm fields and expected value | `run_large_sigz.m`, `analyze_large_sigz.m` |
| Choice of the observation point (single electron) | `single_electron.m`, `scan_lines.m`, `ideal_contamination.m` |
| Paper data and figures | `paper/sec5_slac_application/compute_polarizer_bands.m`, `estimate_resolution_limit.m`, `make_sec5_figures.m` → `figures/fig_slac_polarizer`, `figures/fig_slac_calibration` |

Results: the D1 energy falls from 82 J/sr at σz = 10 nm to 0.29 J/sr at 50 nm, a relative sensitivity of 10–17% per nm. D2 remains sensitive up to about 150 nm. The macro-particle noise floors (2.5×10⁻² and 2.0×10⁻³ J/sr) limit the resolved range to about 70 nm (D1) and 150 nm (D2).

The `matlab/radiation_multi_electron/E7GeV_und_dipole/` folder also contains the exploratory checks made while choosing the setup: `diag_*`, `check_*`, `nokick_check`, `dipole_distance_*`, `noise_theory_check`, and `verify_timeshift_gpt`.
