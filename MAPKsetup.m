Options.SaveFullEvery         = 200000;
Options.MaxIterations         = 1000;
Options.NumOfPosteriorSamples = 20000;
Options.StartingParameters    = [3.8713 0.9196]; 

%%% RMHMC Setup %%%
Options.NumOfLeapFrogSteps = 10;
Options.StepSize           = 6e-2;
Options.NumOfFixPointSteps = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors for parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%

mu=Options.StartingParameters;
sigma=2*ones(1,Options.NumOfParameters);
for i=1:Options.NumOfParameters
 Options.PriorInfo.Type{i}       = 'Normal';
 Options.PriorInfo.Para(i,:)     = [mu(i),sigma(i)];
end%for

U=fliplr({0.92,0.751,0.633,0.359,0.389,0.256,0.197,0.194,0.097});
NumOfObs=length(U);
Data=cell(2,NumOfObs);
measurements=fliplr({0.993,0.789,0.646,0.737,0.798,0.487,0.6,0.689,0.043});
[Data{1,:}]=deal(measurements{:});
[Data{2,:}]=deal(0.2);
