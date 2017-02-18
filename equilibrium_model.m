function m = equilibrium_model
% Basic mass action model with seeds
m = InitializeModelMassActionAmount('Equilibrium');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddSeed(m, 'A_0', 1);
m = AddSeed(m, 'B_0', 2);
m = AddSeed(m, 'C_0', 0);

m = AddState(m, 'A', 'Solution', 'A_0');
m = AddState(m, 'B', 'Solution', 'B_0');
m = AddState(m, 'C', 'Solution', 'C_0');

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, {'binding', 'unbinding'}, {'A', 'B'}, {'C'}, 'kf', 'kr');

m = FinalizeModel(m);