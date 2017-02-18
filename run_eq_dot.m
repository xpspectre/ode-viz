%% Visualize simple equilibrium model
clear; close all; clc
rng('default');

name = '_equilibrium_1';
directory = [name '_out/'];

check_add_path('jsonlab-1.5');

metaFile = [directory 'meta.json'];
meta = loadjson(metaFile);

species = meta.species;
reactions = meta.reactions;
times = meta.times;

% Load stoichiometry matrix
stoichFile = [directory 'stoich.json'];
stoich = loadjson(stoichFile);
S = stoich.stoich;

%% Generate dot file for 1st timepoint
ti = 1;
t = times(ti);

dataFile = [directory 'data' num2str(ti) '.csv'];
data = csvread(dataFile);
species = data(:,1);
fluxes = data(:,2:end);

nx = length(species);
[nx_, nr] = size(fluxes);

assert(nx == nx_, 'run_eq_dot:species_mismatch', 'Number of species from species and fluxes datasets inconsistent')
assert(nx == length(species))
assert(nr == length(reactions))

%% Detect presence of Graphviz cmd line executables
% TODO

%% Calculate species ranges
% Options
logscaleNodeSize = false;

% Defaults
maxNodeSize_ = 40; 
minNodeSize_ = 5;

maxSpecies = max(species);
minSpecies = min(species);
if minSpecies == 0 % Get min nonzero species
    minSpecies = min(species(species > 0));
end

% Scale
speciesScaled = (species - minSpecies) / (maxSpecies - minSpecies);
speciesScaled(species==0) = 0;

%% Calculate flux ranges
% Options
logscaleEdgeWidth = false;

% Defaults
maxEdgeWidth_ = 5;
minEdgeWidth_ = 0.5;
% 0 flux represented by dotted lines

% TODO: Get fluxes using sparse representation thru stoich matrix

maxFlux = max(fluxes(:));
minFlux = min(fluxes(:));
if minFlux == 0 % Get min nonzero flux
    minFlux = min(fluxes(fluxes > 0));
end

% Scale
fluxesScaled = (species - minFlux) / (maxFlux - minFlux);


%% Make gv file
gvFile = [directory 'plot.gv'];
fid = fopen(gvFile);
fprintf(fid, 'digraph model {\n');

% Default node: shape=ellipse,width=.75,height=.5


% First make all the nodes


% Then make all the edges


fprintf(fit, '}');
fclose(fid);

