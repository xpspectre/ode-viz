%% Visualize simple equilibrium model
clear; close all; clc
rng('default');

name = '_equilibrium_1';
directory = [name '_out/'];

check_add_path('jsonlab-1.5');

metaFile = [directory 'meta.json'];
meta = loadjson(metaFile);

speciesNames = meta.species;
reactionNames = meta.reactions;
times = meta.times;

% Load stoichiometry matrix
stoichFile = [directory 'stoich.json'];
stoich = loadjson(stoichFile);
S = stoich.stoich;

%% Generate dot file for 1st timepoint
ti = 10;
t = times(ti);

dataFile = [directory 'data' num2str(ti) '.csv'];
data = csvread(dataFile);
species = data(:,1);
fluxes = data(:,2:end);

nx = length(species);
[nx_, nr] = size(fluxes);

assert(nx == nx_, 'run_eq_dot:species_mismatch', 'Number of species from species and fluxes datasets inconsistent')
assert(nx == length(speciesNames))
assert(nr == length(reactionNames))

%% Detect presence of Graphviz cmd line executables
% TODO

%% Calculate species ranges
% Options
logscaleNodeSize = false;

% Defaults
maxNodeSize_ = 40; 
minNodeSize_ = 5;
midNodeSize_ = 10;

maxSpecies = max(species);
minSpecies = min(species);

if minSpecies == 0 % Get min nonzero species
    minSpecies = min(species(species > 0));
end

% Scale to node size
if minSpecies == maxSpecies % handle corner case of all species the same conc
    speciesScaled = (species > 0) * midNodeSize_;
else
    speciesScaled = (maxNodeSize_ - minNodeSize_) * (species - minSpecies) / (maxSpecies - minSpecies) + minNodeSize_;
end
speciesScaled(species==0) = 0;

%% Calculate flux ranges (abs vals)
% Options
logscaleEdgeWidth = false;

% Defaults
maxEdgeWidth_ = 5;
minEdgeWidth_ = 0.5;
midEdgeWidth_ = 1;
% 0 flux should be represented by dotted lines

% Get fluxes using sparse representation thru stoich matrix
fluxesAbs = abs(fluxes);
maxFlux = -Inf;
minFlux = Inf;
for ix = 1:nx
    mask = S(ix,:) ~= 0;
    maxFlux = max(maxFlux, max(fluxesAbs(ix,mask)));
    minFlux = min(minFlux, min(fluxesAbs(ix,mask)));
end

if minFlux == 0 % Get min nonzero flux - warning: currently expensive
    minFlux = min(fluxesAbs(fluxesAbs > 0));
end

% Scale to edge width
if minFlux == maxFlux % Handle corner case of all fluxes the same
    fluxesScaled = (fluxes ~= 0) * midEdgeWidth_;
else
    fluxesScaled = (maxEdgeWidth_ - minEdgeWidth_) * (fluxesAbs - minFlux) / (maxFlux - minFlux) + minEdgeWidth_;
end

%% Make gv file
gvFile = [directory 'plot.gv'];
fid = fopen(gvFile, 'w');
fprintf(fid, 'digraph model {\n');

% Default node: shape=ellipse,width=.75,height=.5

% First make all the species/nodes
for ix = 1:nx
    conc = speciesScaled(ix);
    if conc == 0 % special handling for 0 species conc
        noconc = ',style=dotted';
        conc = minNodeSize_;
    else
        noconc = '';
    end
    fprintf(fid, '%s [shape=circle,fontsize=%g%s];\n', speciesNames{ix}, conc, noconc);
end

% Make all the reaction/edges
%   Note: that flux already includes stoichiometry
for ir = 1:nr
    fprintf(fid, '%s [shape=box,fontsize=6,width=0,height=0];\n', reactionNames{ir}); % small box for each rxn name
    
    % Reactants connect to box with no arrow
    reactants = find(S(:,ir) < 0);
    for i = 1:length(reactants)
        ind = reactants(i);
        flux = fluxesScaled(ind,ir);
        if flux == 0 % special handling for 0 flux
            noflux = ',style=dotted';
            flux = 1;
        else
            noflux = '';
        end
        fprintf(fid, '%s -> %s [arrowhead=none,penwidth=%g%s];\n', speciesNames{ind}, reactionNames{ir}, flux, noflux);
    end
    
    % Products connected from box with arrow
    products = find(S(:,ir) > 0);
    for i = 1:length(products)
        ind = products(i);
        flux = fluxesScaled(ind,ir);
        if flux == 0 % special handling for 0 flux
            noflux = ',style=dotted';
            flux = 1;
        else
            noflux = '';
        end
        fprintf(fid, '%s -> %s [penwidth=%g%s];\n', reactionNames{ir}, speciesNames{ind}, flux, noflux);
    end
end

fprintf(fid, '}');
fclose(fid);

%% Run Graphviz
outFile = [directory 'plot.pdf'];
s = system(['"C:\Program Files (x86)\Graphviz2.38\bin\neato.exe" -Tpdf ' gvFile ' -o ' outFile]);


