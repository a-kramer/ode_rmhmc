addpath(genpath('./'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set random numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%randn('state', sum(100*clock));
%rand('twister', sum(100*clock));

% Fixed
randn('state', 1);
rand('twister', 1);

make_expSAC_model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load required toolboxes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This part loads the toolbox for solving the ODEs quickly %%%

JobID = 1;
TaskID = 1;
JobName = 'expSAC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options to pass to MCMC file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Options.MaxIterations         = 1000;
Options.NumOfPosteriorSamples = 20000;

Options.GraphicalOutput       = false;

% Set name for saving results
Options.EquationName          = [JobName '_' num2str(JobID) '_' num2str(TaskID)];


%%% Set starting values %%%

Options.InferICs                  = false;
% neutral starting point:
% Options.StartingParameters        = [1.4099, -1.3627, 0.6556, 0.9393, 7.2602,-8.1886, -1.0579,-2.7875, 3.2600, 2.9585];
% upper steady state starting point (Data Cluster 1):
Options.StartingParameters        = [ 1.7766   -2.1769    0.4349    0.2171    1.5836    1.7244   -1.3622   -3.7200    4.4530    2.9478];
% lower steady state starting point (Data Cluster 2):
%Options.StartingParameters        = [1.4099, -1.3627, 0.6556, 0.9393, 7.2602,-8.1886, -1.7,-2.7875, 3.9, 2.9585];


%%% RMHMC Setup %%%

Options.NumOfLeapFrogSteps = 20;
Options.StepSize           = 3e-3;
Options.NumOfFixPointSteps = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors for parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%
%
mu=Options.StartingParameters;
sigma=4*ones(1,F.m);
for i=1:F.m
 Options.PriorInfo.Type{i}       = 'Normal';
 Options.PriorInfo.Para(i,:)     = [mu(i),sigma(i)];
end%for
%AB7=[-1.7930,0.5717];
%AB9=[-0.1116,3.9108];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U={[1;1],[0.5479;1],[0.2407;1],[0.308;1],[1;2],[0.5479;2],[0.2407;2],[0.308;2]};
NumOfObs=length(U);
Data=cell(2,NumOfObs);

% Data Cluster 1
%measurements={0.0017,0.0275,0.0061,0.0232,0.0057,0.0417,0.0226,0.0461};
%SD={0.0015, 0.0664, 0.0118, 0.0514, 0.0129, 0.0815, 0.0405, 0.0742};

%Data Cluster 2
measurements={16.2572   16.8878   16.8416   16.8216   16.5654   17.7860   17.1745   17.6061};
SD={0.2665    1.5472    0.3190    1.2973    0.4331    0.7030    0.4297    0.7012};

[Data{1,:}]=deal(measurements{:});
[Data{2,:}]=deal(SD{:});

       
% Call main population MCMC routine
steadystate_RMHMC_Single_MCMC(Data, U, F, Options);



