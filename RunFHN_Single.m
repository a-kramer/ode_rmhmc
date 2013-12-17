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

% FOR MAC
%{
% Install SBToolbox2 'quick'
CurrentDir = cd;
cd('/Applications/Matlab_Addons/SBTOOLBOX2')
installSB('quick')
cd(CurrentDir)
% Install add-on for fast ODE solving
cd('/Applications/Matlab_Addons/SBPD')
installSBPD('quick')
cd(CurrentDir)
%}

JobID = 1;
TaskID = 1;
JobName = 'FHN';

% Job names should be LockeModel_<number of time points>
%
% e.g. "LockeModel" runs Locke model with 100 time points
%SplitJobName  = regexp(JobName, '_', 'split');

%NumOfTimePoints       = str2num(SplitJobName{2});




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options to pass to MCMC file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Options.NumOfPopulations      = 1;

Options.SaveFullEvery         = 200000;
Options.MaxIterations         = 100;
Options.NumOfPosteriorSamples = 1000;

Options.GraphicalOutput       = false;

% Set name for saving results
Options.EquationName          = [JobName '_' num2str(JobID) '_' num2str(TaskID)]


Options.NumOfParameters       = 3;
Options.ObservedSpecies       = [1 2];
Options.UnobservedSpecies     = [];
Options.SBModelName           = 'SBFitzHughNagumo_2nd'; %'SBFitzHughNagumo';
Options.SBModelParameterNames = SBparameters(Options.SBModelName)';


%%% Set starting values %%%

Options.InferICs                  = false;
Options.StartingInitialConditions = [-1 1];
Options.StartingParameters        = [0.2 0.2 3]; % True parameters


%%% RMHMC Setup %%%

Options.NumOfLeapFrogSteps = 5;
Options.StepSize           = 2/Options.NumOfLeapFrogSteps;
Options.NumOfFixPointSteps = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors for parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Options.PriorInfo.Type{1} = 'Gamma';
Options.PriorInfo.Para(1,:) = [1 2];

Options.PriorInfo.Type{2} = 'Gamma';
Options.PriorInfo.Para(2,:) = [1 2];

Options.PriorInfo.Type{3} = 'Gamma';
Options.PriorInfo.Para(3,:) = [1 2];
%}
%
Options.PriorInfo.Type{1}       = 'Uniform';
Options.PriorInfo.UpperBound(1) = 4;
Options.PriorInfo.LowerBound(1) = 0;

Options.PriorInfo.Type{2}       = 'Uniform';
Options.PriorInfo.UpperBound(2) = 4;
Options.PriorInfo.LowerBound(2) = 0;

Options.PriorInfo.Type{3}       = 'Uniform';
Options.PriorInfo.UpperBound(3) = 4;
Options.PriorInfo.LowerBound(3) = 0;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors for intial conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options.PriorInfoIC.Type{1}   = 'Normal';
Options.PriorInfoIC.Para(1,:) = [-1, 1];

Options.PriorInfoIC.Type{2}   = 'Normal';
Options.PriorInfoIC.Para(2,:) = [1, 1];
%{
Options.PriorInfoIC.Type{1}       = 'Uniform';
Options.PriorInfoIC.UpperBound(1) = 0;
Options.PriorInfoIC.LowerBound(1) = -2;

Options.PriorInfoIC.Type{2}       = 'Uniform';
Options.PriorInfoIC.UpperBound(2) = 2;
Options.PriorInfoIC.LowerBound(2) = 0;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate data
Model             = 'SBFitzHughNagumo';
ParameterNames    = SBparameters(Options.SBModelName)';
InitialConditions = [-1 1];
Parameters        = [0.2 0.2 3];
StartTime         = 0;
EndTime           = 20;
D                 = 100;
TimePoints        = StartTime:(EndTime-StartTime)/(D-1):EndTime;

simdata           = SBPDsimulate(Model, TimePoints, InitialConditions, ParameterNames, Parameters, {});
TimeData          = simdata.time;
TrueData          = simdata.statevalues';

% Choose observed species
NoisyData         = TrueData;

[NumOfSpecies,D]  = size(NoisyData);

% Add noise based on standard deviation of each species
NoisyData         = NoisyData + randn(NumOfSpecies, D).*repmat(std(NoisyData')',1,D)*0.01;

for i = Options.ObservedSpecies
    Options.NoiseCov{i} = eye(D)*var(NoisyData(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% Open saved data set
ExpData        = open(['Locke2005_NoLight_200dps.mat']);
NoisyData      = ExpData.LockeDataNoisy;
TimeData       = ExpData.LockeTimeData;
%}



        
% Call main population MCMC routine
RMHMC_Single_MCMC(NoisyData, TimeData, Options);


