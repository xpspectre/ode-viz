function m = inputs_model
% Basic mass action model with states and inputs, seeds
m = InitializeModelMassActionAmount('Inputs');

m = AddCompartment(m, 'v', 3, 1);

m = AddSeed(m, 'B_0', 2.5);
m = AddSeed(m, 'C_0', 0);

m = AddState(m, 'B', 'v', 'B_0');
m = AddState(m, 'C', 'v', 'C_0');

m = AddInput(m, 'A', 'v', 1);

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, {'binding', 'unbinding'}, {'A', 'B'}, {'C'}, 'kf', 'kr');

m = AddOutput(m, 'A', 'A');
m = addStatesAsOutputs(m);

m = FinalizeModel(m);