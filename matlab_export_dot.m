function [ output_args ] = matlab_export_dot(directory, opts)
% Turn data+metadata output into dot file for Graphviz
%
% Inputs:
%   directory [ string ]
%       Location of exported data and metadata files from model and
%       simulations. Trailing path separator is optional.
%   opts [ struct ]
%       Options struct with fields:
%       .OutputDir [ string {directory} ]
%           Directory to place output dot files. Default is the same as the
%           input directory. Trailing path separator is optional.
%       .tInds [ nti x 1 integer vector | {'all'} ]
%           Time indices to perform plots at. Default is the string 'all',
%           which makes an image for each time in the available data
%       .CombinedBounds [ {true} | false ]
%           Whether to use the max/min of species concs and fluxes of all
%           times for a global scaling, or not (each time gets its own scaling)
%       .LogscaleNodeSize [ true | {false} ]
%           Whether to use log scale on node/species conc size. Useful if
%           concs over many orders of magnitudes are significant. Default
%           false uses linear scale. TODO: Implement logscale.
%       .LogscaleEdgeWidth [ true | {false} ]
%           Whether to use log scale on edge/flux widths. Useful if
%           fluxes over many orders of magnitudes are significant. Default
%           false uses linear scale. TODO: Implement logscale.
%       .MakePlots [ {true} | false ]
%           Whether to run the Graphviz programs to make the graphs (true)
%           or just generate dot files.
%       .PlotFiletype [ {'pdf'} | 'png' | ... ]
%           Plot filetype that Graphviz supports.
%       .PlotFunction [ {'neato'} | 'dot' | ... ]
%           Graphviz plot executable
if nargin < 2
    opts = [];
end

% Default options
opts_ = [];
opts_.OutputDir = directory;
opts_.tInds = 'all';
opts_.CombinedBounds = true;
opts_.LogscaleNodeSize = false;
opts_.LogscaleEdgeWidth = false;
opts_.MakePlots = true;
opts_.PlotFiletype = 'pdf';
opts_.PlotFunction = 'neato';
opts = mergestruct(opts_, opts); % this is included in kroneckerbio

% Check inputs
if opts.LogscaleNodeSize || opts.LogscaleEdgeWidth
    error('Not implemented yet')
end

% Load metadata
metaFile = fullfile(directory, 'meta.json');
meta = loadjson(metaFile);

speciesNames = meta.species;
reactionNames = meta.reactions;
times = meta.times;
cats = meta.cats;
inputs = meta.inputs;

nx = length(speciesNames);
nr = length(reactionNames);

% Load stoichiometry matrix
stoichFile = fullfile(directory, 'stoich.csv');
stoich = csvread(stoichFile);
S = sparse(stoich(:,1), stoich(:,2), stoich(:,3), nx, nr);

%% Specify times to process
if ischar(opts.tInds) && strcmp(opts.tInds, 'all')
    % keep all times
else
    times = times(opts.tInds);
end
nt = length(times);

%% Collect bounds
speciesBounds = zeros(nt,2); % [max,min] for each time
fluxBounds = zeros(nt,2); % [max,min] for each time
for it = 1:nt
    dataFile = fullfile(directory, ['data' num2str(it) '.csv']);
    data = csvread(dataFile);
    species = data(:,1);
    fluxes = data(:,2:end);
    
    nx = length(species);
    [nx_, nr] = size(fluxes);
    
    assert(nx == nx_, 'run_eq_dot:species_mismatch', 'Number of species from species and fluxes datasets inconsistent')
    assert(nx == length(speciesNames))
    assert(nr == length(reactionNames))
    
    % Get max/min species
    maxSpecies = max(species);
    minSpecies = min(species);
    if minSpecies == 0 % Get min nonzero species
        minSpecies = min(species(species > 0));
    end
    
    % Get max/min fluxes using sparse representation thru stoich matrix
    fluxesAbs = abs(fluxes);
    maxFlux = -Inf;
    minFlux = Inf;
    for ix = 1:nx
        if all(S(ix,:) == 0)
            continue
        end
        
        mask = S(ix,:) ~= 0;
        maxFlux = max(maxFlux, max(fluxesAbs(ix,mask)));
        minFlux = min(minFlux, min(fluxesAbs(ix,mask)));
    end
    if minFlux == 0 % Get min nonzero flux - warning: currently expensive
        minFlux = min(fluxesAbs(fluxesAbs > 0));
    end
    
    speciesBounds(it,:) = [maxSpecies, minSpecies];
    fluxBounds(it,:) = [maxFlux, minFlux];
end

if opts.CombinedBounds
    speciesBounds(:,1) = max(speciesBounds(:,1));
    speciesBounds(:,2) = min(speciesBounds(:,2));
    fluxBounds(:,1) = max(fluxBounds(:,1));
    fluxBounds(:,2) = min(fluxBounds(:,2));
end

%% Make gv files
% Calculate scaled node sizes and edge widths
% Defaults
maxNodeSize_ = 40;
minNodeSize_ = 5;
midNodeSize_ = 10;
maxEdgeWidth_ = 5;
minEdgeWidth_ = 0.5;
midEdgeWidth_ = 1;

% Check if needed
if opts.MakePlots
    if ~check_graphviz
        error('run_eq_dot:graphviz_missing', 'Graphviz executables not found')
    end
end

for it = 1:nt
    dataFile = fullfile(directory, ['data' num2str(it) '.csv']);
    data = csvread(dataFile);
    species = data(:,1);
    fluxes = data(:,2:end);
    
    % Scale to node size
    maxSpecies = speciesBounds(it,1);
    minSpecies = speciesBounds(it,2);
    if minSpecies == maxSpecies % handle corner case of all species the same conc
        speciesScaled = (species > 0) * midNodeSize_;
    else
        speciesScaled = (maxNodeSize_ - minNodeSize_) * (species - minSpecies) / (maxSpecies - minSpecies) + minNodeSize_;
    end
    speciesScaled(species<0) = 0;
    
    % Scale to edge width
    maxFlux = fluxBounds(it,1);
    minFlux = fluxBounds(it,2);
    fluxesAbs = abs(fluxes);
    if minFlux == maxFlux % Handle corner case of all fluxes the same
        fluxesScaled = (fluxes ~= 0) * midEdgeWidth_;
    else
        fluxesScaled = (maxEdgeWidth_ - minEdgeWidth_) * (fluxesAbs - minFlux) / (maxFlux - minFlux) + minEdgeWidth_;
    end
    fluxesScaled(fluxesScaled<0) = 0;
    
    %% Make gv file
    % 0 flux should be represented by dotted lines
    gvFile = fullfile(opts.OutputDir, sprintf('plot%i.gv', it));
    fid = fopen(gvFile, 'w');
    fprintf(fid, 'digraph model {\n');
    
    % Default node: shape=ellipse,width=.75,height=.5
    
    % First make all the species/nodes
    for ix = 1:nx
        conc = speciesScaled(ix);
        if conc <= 0 % special handling for 0 species conc
            noconc = ',style=dotted';
            conc = minNodeSize_;
        else
            noconc = '';
        end
        if inputs(ix) % special handling for inputs
            inputStr = ',color=blue';
        else
            inputStr = '';
        end
        fprintf(fid, '"%s" [shape=circle,fontsize=%g%s%s];\n', speciesNames{ix}, conc, noconc, inputStr);
    end
    
    % Make all the reaction/edges
    %   Note: that flux already includes stoichiometry
    for ir = 1:nr
        fprintf(fid, '"%s" [shape=box,fontsize=6,width=0,height=0];\n', reactionNames{ir}); % small box for each rxn name
        
        % Reactants connect to box with no arrow
        reactants = find(S(:,ir) < 0);
        for i = 1:length(reactants)
            ind = reactants(i);
            flux = fluxesScaled(ind,ir);
            if flux <= 0 % special handling for 0 flux
                noflux = ',style=dotted';
                flux = midEdgeWidth_;
            else
                noflux = '';
            end
            fprintf(fid, '"%s" -> "%s" [arrowhead=none,penwidth=%g%s,splines=curved];\n', speciesNames{ind}, reactionNames{ir}, flux, noflux);
        end
        
        % Products connected from box with arrow
        products = find(S(:,ir) > 0);
        for i = 1:length(products)
            ind = products(i);
            flux = fluxesScaled(ind,ir);
            if flux <= 0 % special handling for 0 flux
                noflux = ',style=dotted';
                flux = midEdgeWidth_;
            else
                noflux = '';
            end
            fprintf(fid, '"%s" -> "%s" [penwidth=%g%s,splines=curved];\n', reactionNames{ir}, speciesNames{ind}, flux, noflux);
        end
        
        % Catalytic species - add to reactants and product
        cat = cats{ir}.x;
        if cat(1) ~= 0
            for ic = 1:length(cat)
                ind = cat(ic);
                fprintf(fid, '"%s" -> "%s" [arrowhead=none,color=blue,splines=curved];\n', speciesNames{ind}, reactionNames{ir});
            end
        end
    end
    
    fprintf(fid, '}');
    fclose(fid);
    
    if opts.MakePlots
        outFile = fullfile(opts.OutputDir, sprintf('plot%i.%s', it, opts.PlotFiletype));
        s = system([opts.PlotFunction ' -T' opts.PlotFiletype ' ' gvFile ' -o ' outFile]);
    end
end

end

