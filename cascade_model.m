function m = cascade_model(N)
% Activation cascade model with user-defined number of stages. A species
% starts in an inactive form Xi. It's then activated into Xi* by X(i-1)*.
% There is an initial activator A which gets consumed. Each Xi* is also
% deactivated over time.
%
% Inputs:
%   N [ positive scalar integer ]
%       Number of stages in the cascade
%
% Outputs:
%   m [ kroneckerbio massaction model ]
% 

m = InitializeModelMassActionAmount(sprintf('Cascade_%i', N));

m = AddCompartment(m, 'v', 3, 1);

m = AddState(m, 'A', 'v', 'A_0');
for iN = 1:N
    m = AddState(m, sprintf('X%i', iN), 'v', sprintf('X%i_0', iN));
    m = AddState(m, sprintf('X%i*', iN), 'v', 0);
end

m = AddSeed(m, 'A_0', 10);
for iN = 1:N
    m = AddSeed(m, sprintf('X%i_0', iN), 10);
end

m = AddParameter(m, 'kf', 6);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, {'on1'}, {'X1', 'A'}, {'X1*'}, 'kf');
for iN = 2:N
    m = AddReaction(m, {sprintf('on%i', iN)}, {sprintf('X%i', iN), sprintf('X%i*', iN-1)}, {sprintf('X%i*', iN), sprintf('X%i*', iN-1)}, 'kf');
end
for iN = 1:N
    m = AddReaction(m, {sprintf('off%i', iN)}, {sprintf('X%i*', iN)}, {sprintf('X%i', iN)}, 'kr');
end

m = addStatesAsOutputs(m);

m = FinalizeModel(m);