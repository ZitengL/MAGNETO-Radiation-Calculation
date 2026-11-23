# MAGNETO

MAGNETO calculates the radiation of electrons from their tracked trajectories. For each trajectory it evaluates the acceleration term of the Liénard–Wiechert field at a fixed observation point and integrates it directly over the emission time. For a bunch, the complex fields of all particles are summed with a common time reference before the squared magnitude is taken, so the interference between electrons and between different magnets is retained.

This repository contains the MATLAB implementation and the calculations of the paper

> Z. Liu, B. Schaap, P. Denham, C. Emma, K. K. Swanson, and P. Musumeci,
> *Forward modeling of a coherent undulator radiation diagnostic for attosecond electron bunches at FACET-II*
> (AAC proceedings, in preparation).

The trajectories are computed with [General Particle Tracer (GPT)](https://www.pulsar.nl/gpt/) through three-dimensional magnetic field maps.

## Method in brief

For an observer at a fixed position `x_obs` (in the same laboratory coordinates as the trajectory), with `R = |x_obs - r(t)|` and `n = (x_obs - r)/R`:

1. The normalized acceleration is computed from the fields recorded by GPT at each sample (relativistic Lorentz force),
   `dβ/dt = q/(γ m c) [E + c β×B − β (β·E)]`.
   This gives the instantaneous acceleration at the same time as `r`, `β`, and `γ`, without finite differences or interpolation.
2. Radiation field: `E_rad = q/(4π ε0 c) · n×[(n−β)×dβ/dt] / ((1 − n·β)^3 R)`.
3. Fourier transform, evaluated directly over the trajectory samples using `dt_obs = (1 − n·β) dt`:
   `Ẽ(ω) ≈ Σ_i w_i E_rad,i (1 − n_i·β_i) exp[−iω(t_i + R_i/c − t_ref)]`,
   where `w_i` are trapezoidal weights in laboratory time (`Δt`, and `Δt/2` at the two end points for a uniform step).
   Any set of frequencies can be evaluated. There is no observer-time grid, no interpolation, and no FFT.
4. Spectrum per unit solid angle as seen from the coordinate origin: `S(ω) = ε0 c R0² |Ẽ(ω)|² / π`, with `R0 = |x_obs|`.
5. Bunches: the field of particle `j` is multiplied by its weight `n_j` (physical electrons per macro-particle), and the complex fields are summed with the same `t_ref` before squaring.
   With `n_m` electrons per macro-particle, the incoherent background is `n_m` times the physical one. The calculation therefore resolves coherent radiation only where the coherent signal clearly exceeds this background.

## Repository layout

| Folder | Content |
|---|---|
| [`magneto/`](magneto/) | MAGNETO functions: the radiation calculation and the GPT readers |
| [`examples/`](examples/) | Index of the five calculations of the paper, which folders belong to each, and how to run them |
| [`gpt_sims/`](gpt_sims/) | GPT input files (`.in`) and run scripts (`run.bat`) of every simulation |
| [`matlab/`](matlab/) | MAGNETO calculations on the GPT output, analysis and checks |
| [`paper/`](paper/) | Data (`data/*.txt`), figures and the scripts that collect the data and draw the paper figures |
| [`field_map/`](field_map/) | GPT field maps of the two undulators (the two largest are GitHub Release assets) |
| [`srw/`](srw/) | Synchrotron Radiation Workshop (SRW) benchmark scripts and the SRW spectra used in the paper |

The folder names are kept as in the original working tree so that all relative paths in the scripts remain valid.

## Requirements

- MATLAB (developed with R2024b).
- GPT (commercial, not included). The run scripts read two environment variables:
  - `GPTBIN`: the GPT `bin` directory (containing `gpt.exe`, `gdfa.exe`, `asci2gdf.exe`);
  - `GPTLICENSE`: your GPT license number.
- For the SRW benchmark: Python 3 with SRW, NumPy, and SciPy.

The scripts were run on Windows (`.bat`, `cmd` syntax in `system` calls).

## What is not included

The repository contains inputs, scripts, and processed results, but not the large binary files:

- GPT output (`*.gdf`, up to 22 GB per run);
- two of the three field maps (0.5 and 0.9 GB), which exceed the GitHub file-size limit. They are attached to the GitHub Release of this repository, see [`field_map/README.md`](field_map/README.md);
- intermediate MATLAB results (`*.mat`).

Scripts that read these files therefore need the corresponding GPT runs first. The processed data of every paper figure are in `paper/*/data/*.txt`.

## Quick start

```matlab
addpath('magneto')
traj = read_gpt_trajectory('gpt_sims/single_electron/E240MeV_und15p_paraxial/result.gdf');  % after running GPT
obs  = [0 0 11];                                   % observer [m], laboratory coordinates
lam  = logspace(log10(15e-9), log10(150e-9), 5000);
spec = emission_time_spectrum(traj, obs, lam, 'fB');
loglog(spec.lambda*1e9, spec.d2W_dwdOmega.total)   % J/(rad/s)/sr, per sr seen from the origin
```

See [`magneto/README.md`](magneto/README.md) for the functions and [`examples/README.md`](examples/README.md) for the paper calculations.

## License

- Code: MIT License, see [`LICENSE`](LICENSE).
- Data and figures (`paper/*/data`, `paper/*/figures`, result files): CC BY 4.0, see [`DATA_LICENSE.md`](DATA_LICENSE.md).

If you use MAGNETO or these results, please cite the paper above (see [`CITATION.cff`](CITATION.cff)).
