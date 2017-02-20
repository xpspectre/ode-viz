%% Setup, run, and visualize simple equilibrium model
clear; close all; clc
rng('default');

name = '_equilibrium_1';
buildModel = true;

%% Build or load cached model since model finalization can be expensive
modelName = '_equilibrium_model.mat';
if buildModel
    m = equilibrium_model;
%     m = equilibrium_model_analytic;
    save(modelName)
else
    loaded = load(modelName);
    m = loaded.m;
end

%% Simulate model and save results
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');
tF = 1; % final time
sim = SimulateSystem(m, con, tF);

% Extract states
t = linspace(0, tF, 100);
x = sim.x(t);

% Plot result
figure
plot(t, x)
legend('A','B','C')
xlabel('Time')
ylabel('Amount')
title('Equilibrium Model Traces')

%% Export components for graph visualization
t = linspace(0, tF, 20);
matlab_export_dynamics(m, con, t, name);

