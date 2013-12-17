addpath(genpath('./'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set random numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%randn('state', sum(100*clock));
%rand('twister', sum(100*clock));

% Fixed
randn('state', 1);
rand('twister', 1);

make_expApoptosis_model

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
JobName = 'expApoptosis';

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
Options.MaxIterations         = 1000;
Options.NumOfPosteriorSamples = 50000;

Options.GraphicalOutput       = false;

% Set name for saving results
Options.EquationName          = [JobName '_' num2str(JobID) '_' num2str(TaskID)]


%%% Set starting values %%%

Options.InferICs                  = false;
Options.StartingParameters        = log([5.8e-5,    1e-5,  5e-4, ...
                                           3e-4,  5.8e-3,  5.8e-3, ...
                                        1.73e-2, 1.16e-2,  3.9e-3, 3.9e-3,...
                                           5e-4,    1e-3, 1.16e-2, ...
                                           0.21,     464,     507,...
                                           81.9,    0.21,    40]);


%%% RMHMC Setup %%%

Options.NumOfLeapFrogSteps = 10;
Options.StepSize           = 9e-6;
Options.NumOfFixPointSteps = 10;


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
mu=Options.StartingParameters;
sigma=1*ones(1,F.m);
for i=1:F.m
 Options.PriorInfo.Type{i}       = 'Normal';
 Options.PriorInfo.Para(i,:)     = [mu(i),sigma(i)];
end%for

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U={[0],[1e3],[2e3],[3e3]};
Data = {[0],[5e3],[5e3],[5e3];[1e2],[1e2],[1e2],[1e2]};

       
% Call main population MCMC routine
steadystate_RMHMC_Single_MCMC(Data, U, F, Options);


