%% Benchmarking activation cascade model
clear; close all; clc
rng('default');

name = '_inputs';

if ~exist(name, 'dir')
    mkdir(name);
end

buildModel = true;

%% Build or load cached model since model finalization can be expensive
modelName = '_inputs_model.mat';
if buildModel
    m = inputs_model;
    save(modelName)
else
    loaded = load(modelName);
    m = loaded.m;
end

%% Simulate model and save results
% There's 2 timescales in this model:
%   1. Fast timescale for activation
%   2. Slower inactivation
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');
tF = 1; % final time
sim = SimulateSystem(m, con, tF);

% Extract states
t = linspace(0, tF, 100);
y = sim.y(t);

% Plot result
% figure
% plot(t, y)
% legend({m.Outputs.Name}, 'location','best')
% xlabel('Time')
% ylabel('Amount')
% title('Inputs Model Traces')

%% Export components for graph visualization
t = linspace(0, tF, 10);
matlab_export_dynamics(m, con, t, name);

opts = [];
opts.PlotFunction = 'neato';

matlab_export_dot(name, opts);
