# SRW benchmark

Independent spectra from the [Synchrotron Radiation Workshop (SRW)](https://github.com/ochubar/SRW) for the single-electron comparisons of examples 01 and 02. SRW tracks the electron through the same field map (or the same Enge dipole), with the same initial conditions and observation points as GPT and MAGNETO.

| Folder | Case | Spectrum used in the paper |
|---|---|---|
| `undulator_240MeV_six_spectra_15_150nm/` | 15-period undulator, 240 MeV, θx = 0–0.5 mrad at 11 m | `spectrum_theta_*.txt` |
| `dipole_240MeV_onaxis_15_150nm/` | 0.1 m, 0.65 T dipole, 240 MeV, on axis at 11 m | `dipole240_onaxis_spectrum.txt` |
| `dipole_240MeV_onaxis_15_150nm/dipole240_tangent/` | Same dipole, arc-midpoint tangent observer | `SRW_tangent_enge_11m.txt` |

`track.py` contains the shared tracking setup. The scripts record how the spectra were produced: `calculate.py`, `refine.py` and `export.py` for the spectra, and `analyze.py` / `analysis.py` for the comparisons. They expect a local SRW build (`srwpy`) next to the repository and the GPT output and field maps. They are therefore documentation of the benchmark, not a turnkey script. The environment variable `GPT_ROOT` (GPT installation folder) is used only to record file hashes for provenance.

The spectra files hold columns for both SRW and MAGNETO. The MATLAB scripts in `paper/sec3_single_electron/` read them directly.
