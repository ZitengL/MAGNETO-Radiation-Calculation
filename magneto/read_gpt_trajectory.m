function traj = read_gpt_trajectory(gdf_filepath, particle_index)
%% read_gpt_trajectory
% Reads one particle's trajectory from a GPT tout output file (binary, full
% double precision via load_gdf).
%
% Syntax:
%   traj = read_gpt_trajectory(gdf_filepath)
%   traj = read_gpt_trajectory(gdf_filepath, particle_index)
%
% Inputs:
%   gdf_filepath    - Path to the GPT '.gdf' file written with tout.
%   particle_index  - Row of the particle inside each time group. Default: 1.
%
% Output (struct, all arrays are rows / 3 x N):
%   traj.t       - Time [s], 1 x N
%   traj.r       - Position [x; y; z] [m], 3 x N
%   traj.b       - Velocity/c [Bx; By; Bz], 3 x N
%   traj.G       - Lorentz factor, 1 x N
%   traj.fE      - Electric field at the particle [V/m], 3 x N
%   traj.fB      - Magnetic field at the particle [T], 3 x N
%   traj.m, traj.q, traj.nmacro - Particle mass [kg], charge [C], macro size
%   traj.source  - The file the data came from

if nargin < 2; particle_index = 1; end

% load_gdf.m lives in the project root
this_dir = fileparts(mfilename('fullpath'));
% load_gdf.m is in the same folder

if ~exist(gdf_filepath, 'file'); error('GDF file not found: %s', gdf_filepath); end
data = load_gdf(gdf_filepath);

% Keep only time groups (tout output)
is_time = arrayfun(@(e) isfield(e.p, 'time') && isfield(e.d, 'x'), data);
data = data(is_time);
N = numel(data);

traj.t  = zeros(1, N);
traj.r  = zeros(3, N);
traj.b  = zeros(3, N);
traj.G  = zeros(1, N);
traj.fE = zeros(3, N);
traj.fB = zeros(3, N);
k = particle_index;
for i = 1:N
    d = data(i).d;
    traj.t(i)    = data(i).p.time;
    traj.r(:, i) = [d.x(k); d.y(k); d.z(k)];
    traj.b(:, i) = [d.Bx(k); d.By(k); d.Bz(k)];
    traj.G(i)    = d.G(k);
    traj.fE(:, i) = [d.fEx(k); d.fEy(k); d.fEz(k)];
    traj.fB(:, i) = [d.fBx(k); d.fBy(k); d.fBz(k)];
end
traj.m      = data(1).d.m(k);
traj.q      = data(1).d.q(k);
traj.nmacro = data(1).d.nmacro(k);
traj.source = gdf_filepath;

if any(diff(traj.t) <= 0); error('Time groups are not strictly increasing in %s', gdf_filepath); end
end
