addpath(genpath('./'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set random numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%randn('state', sum(100*clock));
%rand('twister', sum(100*clock));

% Fixed
randn('state', 1);
rand('twister', 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load required toolboxes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This part loads the toolbox for solving the ODEs quickly %%%


JobID = 1;
TaskID = 1;
JobName = 'SB_HMC_Mifa';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options to pass to MCMC file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Options

Options.NumOfPopulations      = 1;
Options.GraphicalOutput       = false;

% Set name for saving results
Options.EquationName          = [JobName '_' num2str(JobID) '_' num2str(TaskID)];

Options.NumOfSpecies          = 6;
Options.NumOfParameters       = 14;
Options.NumOfOutputs          = 1;
Options.OutputMatrix          = [0, 0, 0, 0, 13740, 0];
Options.SBModelName           = 'SBMifa_for_HMC';
Options.SBModelParameterNames = SBparameters(Options.SBModelName)';


MifaSetup
       
% Call main population MCMC routine
SB_HMC_Single_MCMC(Data, U, Options);


