%Options.SaveFullEvery         = 10000;
Options.MaxIterations         = 1000;
Options.NumOfPosteriorSamples = 10000;


%%% Set starting values %%%

Options.InferICs                  = false;%              14
Options.StartingParameters        = [-2.3152   -3.1816  -6.0268  -21.1434    8.5476   10.5344   -0.1127    0.5159   -0.7367   -0.9322    6.7211   -1.1982   -1.4844   21.1269];
%{
[...
-0.9441759354, ...
-4.4228486292, ...
-3.9120230054, ...
-10.0351144741, ...
2.9993756285, ...
5.6584328849, ...
1.2900587532, ...
0.5140205147, ...
-3.3242363405, ...
3.4841282224, ...
0.4879663296, ...
-2.1803674603, ...
-1.0106014113, ...
7.3656091479];
%}
%%% RMHMC Setup %%%

Options.NumOfLeapFrogSteps = 20;
Options.StepSize           = 2e-2; %1.8e-2
%Options.StepSize           = 1.8e-2; % for the last row in Table 1
Options.NumOfFixPointSteps = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors for parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%

mu=Options.StartingParameters;
sigma=4*ones(1,Options.NumOfParameters);
sigma([1,2,4,7,8,10,12,14])=0.01*ones(1,8);
for i=1:Options.NumOfParameters
 Options.PriorInfo.Type{i}       = 'Normal';
 Options.PriorInfo.Para(i,:)     = [mu(i),sigma(i)];
end%for

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measurements=[18.3851,29.1153,26.3649,28.2924,42.3078,105.0286,94.9714];
sd=[20.3221,29.8399,31.7904,32.9355,16.5959,+23.3895,23.3895];

measurements=num2cell(measurements);
sd=num2cell(sd);

U={0,0.01,0.1,0.3,1,10,100};
NumOfObs=length(U);

Data=cell(2,NumOfObs);
[Data{1,:}]=deal(measurements{:});
[Data{2,:}]=deal(sd{:});
