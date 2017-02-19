function [species, fluxes, stoich] = matlab_export_dynamics(m, con, t, name, opts)
% Simulate model+experiment and write out files for visualizer. Also
% returns species concs and fluxes for convenience.
%
% Inputs:
%   m [ kroneckerbio model ]
%       Finalized model with k's
%   con [ kroneckerbio experiment ]
%       Experiment with s,q,h's
%   t [ 1 x nt double vector ]
%       Times to calculate model results at. Should be monotonically
%       increasing.
%   name [ string ]
%       Base/prefix name of output files
%   opts [ struct ]
%       Options struct with fields:
%       .ReactionNames [ 1 x nr cell vector of strings {{}} ]
%           Override reaction names that will be added to metadata
% Outputs:
%   species [ nx+nu x nt double matrix ]
%       Species concs over time
%   fluxes [ nx+nu x nr x nt double matrix ]
%       Reaction fluxes in to/out of each species for each reaction over time
%   stoich [ nx+nu x nr sparse double matrix ]
%       Stoichiometry matrix
%
% Side Effects:
%   Outputs 1 csv data file for each time
%   where the 1st col is the species conc and the subsequent cols are the
%   fluxes of each species thru each rxn. Outputs 1 json file for
%   metadata. csv data files are indexed by time index. Outputs 1 json file
%   for the stoichiometry matrix in sparse form, which is needed to
%   represent all reactions, including those w/ 0 flux at some points.

if nargin < 5
    opts = [];
end

% Default options
opts_ = [];
opts_.ReactionNames = {};
opts = mergestruct(opts_, opts); % this is included in kroneckerbio

% Checks
t = reshape(t,1,length(t)); % Make sure this is a row vector

% Get dynamics
[species, fluxes] = matlab_extract_dynamics(m, con, t, opts);

% Checks
nr = size(fluxes,2);
if ~isempty(opts.ReactionNames)
    assert(length(opts.ReactionNames) == nr, 'matlab_export_dynamics:rxn_names_mismatch', 'Number of reaction names must match number of reactions (2nd dim in fluxes) if provided')
end

% Create output directory
dirName = [name '_out/'];
if ~exist(dirName, 'dir')
    mkdir(dirName);
end

% Slice and dice data to per-time matrices and save
[nx, nt] = size(species); % treat states+inputs as just states
for it = 1:nt
    dataName_i = [dirName 'data' num2str(it) '.csv'];
    data_i = zeros(nx,1+nr);
    data_i(:,1) = species(:,it);
    data_i(:,2:end) = fluxes(:,:,it);
    csvwrite(dataName_i, data_i);
end

% Write metadata
check_add_path('jsonlab-1.5');
metaFile = [dirName 'meta.json'];

meta = [];
meta.model_name = m.Name;
meta.species = [{m.States.Name}, {m.Inputs.Name}];
meta.reactions = {m.Reactions.Name};
meta.times = t; % must be a row vector

savejson('', meta, metaFile);

% Write stoichiometry matrix - coordinate form sparse matrix
stoichFile = [dirName 'stoich.csv'];
[i,j,v] = find(m.S);
csvwrite(stoichFile, [i,j,v]);
end

