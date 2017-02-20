%% Benchmarking activation cascade model
clear; close all; clc
rng('default');

N = 5;
name = sprintf('_cascade_%i', N);

if ~exist(name, 'dir')
    mkdir(name);
end

buildModel = true;

%% Build or load cached model since model finalization can be expensive
modelName = '_cascade_model.mat';
if buildModel
    m = cascade_model(N);
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
tF = 0.1; % final time
sim = SimulateSystem(m, con, tF);

% Extract states
t = linspace(0, tF, 100);
x = sim.x(t);

% Plot result
% figure
% plot(t, x)
% legend({m.States.Name}, 'location','best')
% xlabel('Time')
% ylabel('Amount')
% title('Cascade Model Traces')

%% Export components for graph visualization
t = linspace(0, tF, 10);
matlab_export_dynamics(m, con, t, name);

opts = [];
opts.PlotFunction = 'dot';

matlab_export_dot(name, opts);
