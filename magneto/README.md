# MAGNETO functions

| Function | Purpose |
|---|---|
| `emission_time_spectrum.m` | Radiation spectrum of one trajectory at one observation point (the MAGNETO calculation) |
| `read_gpt_trajectory.m` | Reads one particle's trajectory from a GPT `tout` output file |
| `read_gpt_particles.m` | Streams selected particles (by GPT ID) from a large multi-particle GPT output, for batch processing |
| `particle_traj.m` | Converts one particle of a `read_gpt_particles` batch into the trajectory format of `emission_time_spectrum` |
| `harmonic_wavelengths.m` | On-axis resonant wavelengths of an undulator, estimated from a GPT trajectory |
| `load_gdf.m` | Low-level reader of the binary GPT GDF format (full double precision) |

Add this folder to the MATLAB path: `addpath('magneto')`.

## `emission_time_spectrum`

```matlab
spec = emission_time_spectrum(traj, obs_point, lambda, bdot_method, t_ref)
```

| Input | Meaning |
|---|---|
| `traj` | Trajectory struct (`t`, `r`, `b`, `G`, `fE`, `fB`, `q`, `m`, `nmacro`) from `read_gpt_trajectory` or `particle_traj` |
| `obs_point` | Observer position `[x y z]` [m], in the laboratory coordinates of the trajectory |
| `lambda` | Wavelengths [m], any vector |
| `bdot_method` | `'fB'` (default): acceleration from the Lorentz force with the fields recorded by GPT. `'central'` and `'forward'` (finite differences) are kept for comparison |
| `t_ref` | Common observer-time reference [s]. Required when the fields of several particles are summed; the paper calculations use `t_ref = norm(obs_point)/c` |

| Output | Meaning |
|---|---|
| `spec.Ew` | Fourier-domain field `Ẽ(ω)` [V s/m], 3 × Nω (x; y; z). Already includes the charge `q · nmacro` of the (macro-)particle |
| `spec.d2W_dwdOmega.x/.y/.z/.total` | `ε0 c R0² |Ẽ|² / π` [J/(rad/s)/sr], per sr as seen from the coordinate origin, `R0 = |obs_point|` |
| `spec.lambda`, `spec.omega`, `spec.R_obs` | Wavelength, angular frequency, `R0` |

### Bunches

For several particles, call `emission_time_spectrum` for every particle with the same `obs_point`, `lambda`, and `t_ref`. Then add the complex fields `spec.Ew` and compute

```matlab
S = eps0*c0*norm(obs_point)^2/pi * abs(Ew_sum).^2;      % per component
```

A polarizer is a projection of `Ew_sum` before squaring. Do not multiply `Ew` by `nmacro` again; it is already included.

### Accuracy

The accuracy is set by the GPT output step: the phase `ω (1 − n·β) Δt` between samples must be small. Check convergence by halving the step. The `tout_convergence/` folders in `gpt_sims/` do this for the paper runs.

## Provenance

`load_gdf.m` was copied from the local `load_gdf` function of an earlier group script (`ATF_IFEL.m`, Brookhaven ATF undulator project).
