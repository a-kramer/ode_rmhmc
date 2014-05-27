addpath(genpath('./'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set random numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%randn('state', sum(100*clock));
%rand('twister', sum(100*clock));

% Fixed
randn('state', 1);
rand('twister', 1);

make_Mifa_model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load required toolboxes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This part loads the toolbox for solving the ODEs quickly %%%


JobID = 1;
TaskID = 1;
JobName = 'Mifa'; % for newton raphson type sampler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options to pass to MCMC file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Options

Options.NumOfPopulations      = 1;
Options.GraphicalOutput       = false;
Options.NumOfParameters       = F.m;
% Set name for saving results
Options.EquationName          = [JobName '_' num2str(JobID) '_' num2str(TaskID)];


MifaSetup

% Call main population MCMC routine
steadystate_SMMALA_Single_MCMC(Data, U, F, Options);
