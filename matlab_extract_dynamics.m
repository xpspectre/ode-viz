function [species, fluxes] = matlab_extract_dynamics(m, con, t, opts)
% Simulate model and export components for visualization. States and inputs
% are combined as [states; inputs]. NOTE: Requires the model be finalized
% and the model and experiment filled in with the desired parameter values.
%
% Inputs:
%   m [ kroneckerbio model ]
%       Finalized model with k's
%   con [ kroneckerbio experiment ]
%       Experiment with s,q,h's
%   t [ nt x 1 double vector ]
%       Times to calculate model results at. Should be monotonically
%       increasing.
%   opts [ struct ]
%       Options struct with fields:
%       .
%
% Outputs:
%   species [ nx+nu x nt double matrix ]
%       Species concs over time
%   fluxes [ nx+nu x nr x nt double matrix ]
%       Reaction fluxes in to/out of each species for each reaction over time

if nargin < 4
    opts = [];
end

% Default options
% opts_ = [];
% opts = mergestruct(opts_, opts); % this is included in kroneckerbio

% Checks
assert(m.Ready, 'matlab_ode_vis_export:unfinalized_model', 'Model must be finalized')

% Simulate model+experiment
sim = SimulateSystem(m, con, t(end));

% Assemble time traces of species concs
x = sim.x(t);
u = sim.u(t);
species = [x; u];

% Assemble fluxes at each time
nt = length(t);
S = m.S; % stoichiometry matrix, [ nx+nu x nr double matrix ] (often denoted N)
r = m.r; % reaction (rates), function (t,x,u) -> [ nr x 1 double vector ] (often denoted v)
fluxes = zeros(m.nx+m.nu, m.nr, nt);
for it = 1:nt
    rti = r(t(it), x(:,it), u(:,it))'; % Reaction rates at t, transposed to match rows of S
    fluxes(:,:,it) = bsxfun(@times, S, rti); % expands down the rows of nx+nu, preserves sparsity of S
end

end

