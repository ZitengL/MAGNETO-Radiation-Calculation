function traj = particle_traj(tr, k)
%% particle_traj
% One particle of a read_gpt_particles batch, in the trajectory format used by
% emission_time_spectrum.m (same fields as read_gpt_trajectory.m).
%
% Syntax:
%   traj = particle_traj(tr, k)     % k = row in the batch (not the GPT ID)

ok = ~isnan(tr.z(k, :));
traj.t  = tr.t(ok);
traj.r  = [tr.x(k, ok); tr.y(k, ok); tr.z(k, ok)];
traj.b  = [tr.Bx(k, ok); tr.By(k, ok); tr.Bz(k, ok)];
traj.G  = tr.G(k, ok);
if size(tr.fEx, 1) == 1 && size(tr.x, 1) > 1          % fE skipped by read_gpt_particles (zeros)
    traj.fE = zeros(3, nnz(ok));
else
    traj.fE = [tr.fEx(k, ok); tr.fEy(k, ok); tr.fEz(k, ok)];
end
traj.fB = [tr.fBx(k, ok); tr.fBy(k, ok); tr.fBz(k, ok)];
traj.q  = tr.q(k);  traj.m = tr.m(k);  traj.nmacro = tr.nmacro(k);
traj.ID = tr.ID(k);
traj.source = tr.source;
end
