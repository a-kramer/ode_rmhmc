%Options.SaveFullEvery         = 10000;
Options.MaxIterations         = 1000;
Options.NumOfPosteriorSamples = 10000;


%%% Set starting values %%%

Options.InferICs                  = false;
Options.StartingInitialConditions = [10 0 0];
%{
Options.StartingParameters        = log([...
0.05,...
1e-6,...
1,...
0.5,...
1,...
30]);
%}
Options.StartingParameters        = [ 0.3460    0.4023    5.1190    2.3106   -9.3211   -5.4594];
%%% RMHMC Setup %%%

Options.NumOfLeapFrogSteps = 28;
Options.StepSize           = 2.3e-2;
Options.NumOfFixPointSteps = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors for parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%

mu=zeros(1,Options.NumOfParameters);
sigma=6*ones(1,Options.NumOfParameters);
for i=1:Options.NumOfParameters
 Options.PriorInfo.Type{i}       = 'Normal';
 Options.PriorInfo.Para(i,:)     = [mu(i),sigma(i)];
end%for

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measurements=[18.3851,29.1153,26.3649,28.2924,42.3078,105.0286,94.9714];
sd=[20.3221,29.8399,31.7904,32.9355,16.5959,23.3895,23.3895];

measurements=num2cell(measurements);
sd=num2cell(sd);

U={0,0.01,0.1,0.3,1,10,100};
NumOfObs=length(U);

Data=cell(2,NumOfObs);
[Data{1,:}]=deal(measurements{:});
[Data{2,:}]=deal(sd{:});
