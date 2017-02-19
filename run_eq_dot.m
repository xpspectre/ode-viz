%% Visualize simple equilibrium model
clear; close all; clc
rng('default');

name = '_equilibrium_1';
directory = [name '_out/'];

outputDir = [name '_plots/'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

opts = [];
opts.OutputDir = outputDir;
opts.PlotFiletype = 'png';

matlab_export_dot(directory, opts);

