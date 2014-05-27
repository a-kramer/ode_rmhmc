function [] = SB_SMMALA_Single_MCMC(Data, U, Options )
% Copyright 2012-2013 Ben Calderhead <b.calderhead@imperial.ac.uk>, Andrei Kramer <andrei.kramer@ist.uni-stuttgart.de>

%warning off
Options
% Start the timer
tic

%rand('twister', sum(100*clock))
%randn('state', sum(100*clock))
rand('twister', 1)
randn('state', 1)
T=600; % t of steady state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise user options if not already specified                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get N (number of chemical species) and D (number of time points)
[N, D] = size(Data);
NumOfObs     = D;

% Get specified options
OutputMatrix          = Options.OutputMatrix;
NumOfSpecies          = size(Options.OutputMatrix,2);
NumOfParameters       = Options.NumOfParameters;
NumOfOutputs          = size(OutputMatrix,1);
SBModelName           = Options.SBModelName;

%SpeciesObserved       = Options.ObservedSpecies;
%SpeciesUnobserved     = Options.UnobservedSpecies;

%N = length(SpeciesObserved) + length(SpeciesUnobserved);

MaxIterations         = Options.MaxIterations;
NumOfPosteriorSamples = Options.NumOfPosteriorSamples;

PriorInfo             = Options.PriorInfo;

EquationName          = Options.EquationName;

Parameters            = Options.StartingParameters;

NumOfSens=NumOfSpecies*NumOfParameters;
SensRange=[NumOfSpecies+1:NumOfSpecies+NumOfSens];

StepSize           = Options.StepSize;
options.abstol     = 1e-8;
options.reltol     = 1e-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise non changeable stuff                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup initial parameters
if isempty(Parameters)
    for PopulationNum = 1:NumOfPopulations
        for n = 1:NumOfParameters
            Parameters(n) = GetPrior(PriorInfo, n, 'random');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup initial values for solving the ODE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
if isempty(InitialValues)
    for PopulationNum = 1:NumOfPopulations
        for n = 1:N
            InitialValues(n) = GetPrior(PriorInfoIC, n, 'random');
        end
    end
end
%}

% Set up proposal counters
AcceptedMutation  = 0;
AttemptedMutation = 0;

% Set up parameter history variable
ParaHistory         = [];
LLHistory           = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up initial noise for likelihood function in each population %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invert the noise covariance matrices
% for i = SpeciesObserved
%     NoiseCovInv{i}    = inv(NoiseCov{i} + eye(D)*1e-8);
%     NoiseCovLogDet(i) = 2*sum(log(diag(chol(NoiseCov{i} + eye(D)*1e-8))));
% end
NoiseCovLogDet = 2*sum(log(cat(1,Data{2,:})));


% Set monitor rate for adapting step sizes
MonitorRate = 100;


% Set up converged flag
ContinueIterations = true;
Converged          = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precalculate some values for speed                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Precalculate the log likelihoods for each chain and population
% For each population
%try
    ySensitivities=cell(1,NumOfObs);
    [ySensitivities{:}]=deal(zeros(NumOfOutputs,NumOfParameters));
    
    
    Y=cell(1,NumOfObs);
    GradL = zeros(1, NumOfParameters);
    G = zeros(NumOfParameters);
    for j=1:NumOfObs
        % [InitialValues zeros(1,NumOfParameters*NumOfSpecies) ...
        % zeros(1,NumOfSpecies*NumOfParameters*NumOfParameters)]
        simdata = feval(SBModelName, [0,T], [], [Parameters, U{j}'], options);
        
        % Get sensitivities of all species with respect to parameter j
        Sensitivities = reshape(simdata.statevalues(2,SensRange),NumOfSpecies,NumOfParameters);

        ySensitivities{j}=OutputMatrix*Sensitivities;
        % Get species time series
        DataTemp         = simdata.statevalues(2,1:NumOfSpecies)';
        XEstimates       = DataTemp;
        Y{j}=OutputMatrix*XEstimates;
        % Calculate gradients for each of the parameters i.e. d(LL)/d(Parameter)
        
        
        GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*ySensitivities{j};
        
        G = G + ySensitivities{j}'*bsxfun(@rdivide,ySensitivities{j},Data{2,j}.^2);
    end
     
    % calculate GradLL using finite differences to be sure
    s=[+1,-1];
    h=1e-8; % finite difference scale
    GL=zeros(1,NumOfParameters);
    for i=1:2        % sign loop: +h and -h
        for k=1:NumOfParameters
            o=Parameters;
            o(k)=o(k)+s(i)*h;
            %calculate logLikelihood
            ll=0;
            for j=1:NumOfObs
                %ySensitivities{j}
                CX=Y{j}+ySensitivities{j}*(o-Parameters)';
                %fprintf('norm(F.f(xs{j},theta,U{j}))=%g\n',norm(F.f(xs{j},theta,U{j})));
                %fprintf('norm(F.f(X,o,U{j}))=%g\n',norm(F.f(X,o,U{j})));
                
                ll = ll - 0.5*sum(((CX-Data{1,j})./Data{2,j}).^2);                
            end
            fprintf('[%i] ll(o(%i)%+i*h)=%.10f\n',i,k,s(i),ll);
            %calculate finite differences:
            GL(k)=GL(k)+s(i)*ll/(2*h);
        end
    end
    fprintf('relative difference (norm) between GradL and its finite difference approximation: %g\n',norm((GradL-GL)./GradL));    
    
    g={GradL,GL}
   
    for ParaNum  = 1:NumOfParameters
        GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, Parameters(ParaNum));
    end

    G = G - diag(GetPriorLog2ndDeriv(PriorInfo, Parameters));
    
    % Save current metric tensor
    CurrentG     = G;
    CurrentCholG = chol(G); 
    CurrentGradL = GradL;
    
    LL=zeros(1,NumOfObs);
    for j=1:NumOfObs
        %YDdifference_j=((Y{j}-Data{1,j})./Data{2,j}).^2
        LL(j) = - 0.5*sum(((Y{j}-Data{1,j})./Data{2,j}).^2);
    end
    %LL
    CurrentLL = - 0.5*(NoiseCovLogDet + NumOfOutputs*NumOfObs*log(2*pi)) + sum(LL);
    
    % make a small parameter step and verify GradL
    fprintf('some gradient(L) tests using finite difference approximations (h=%g):\n',h);
    
    for a=linspace(1,100,2)
        d=a*2e-6*GL'/norm(GL);
        o=Parameters'+d;
        ll=0;
        for j=1:NumOfObs
            CX=Y{j}+ySensitivities{j}*d; %Sf{j} is verifiably correct
            ll = ll - 0.5*sum(((CX-Data{1,j})./Data{2,j}).^2);
        end
        
        fprintf('(1) logLikelihood(theta+d)=%.10f\n',ll);
        fprintf('(2) logLikelihood(theta){%g} + GradL*d{%g} = %.10f\n',sum(LL),g{1}*d,sum(LL)+g{1}*d);
        fprintf('(3) logLikelihood(theta){%g} + GradL_h*d{%g} = %.10f\n',sum(LL),g{2}*d,sum(LL)+g{2}*d);
        fprintf('(2)-(1) = %g\n',sum(LL)+g{1}*d-ll);
        fprintf('(3)-(1) = %g\n',sum(LL)+g{2}*d-ll);
    end

    
    
%catch
    
%     for n=1:SpeciesObserved
%         CurrentLL(n) = -Inf;
%     end
    
%end
ProposedLogPrior=zeros(1,NumOfParameters);
CurrentLogPrior=zeros(1,NumOfParameters);

% define functions for convenience:


WaitbarHandle=waitbar(0,sprintf('Sample Size: %i',NumOfPosteriorSamples));
%onCleanup(@() close(WaitbarHandle));
disp('Initialisation Completed..');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population MCMC Algorithm to sample the parameters based on likelihood  %
% of X given current parameters                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise iteration number
IterationNum = 0;

% Main loop
while ContinueIterations
    
    % Increment iteration number
    IterationNum = IterationNum + 1;
    
    OriginalParas = Parameters;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate parameter values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 
    
    AttemptedMutation = AttemptedMutation + 1;
    % parameter update
             ProposalMean = OriginalParas + 0.5*StepSize^2*(CurrentG\CurrentGradL')';
    NewParas = ProposalMean + StepSize*(CurrentCholG\randn(NumOfParameters,1))';

    % recalculate outputs using model trajectories:
    G = zeros(NumOfParameters);
    for j=1:NumOfObs
        %                                                                %[InitialValues ...
        %                                                                %   zeros(1,NumOfParameters*NumOfSpecies) ...
        %                                                                %   zeros(1,NumOfSpecies*NumOfParameters*NumOfParameters)]
        simdata          = feval(SBModelName, [0,T], [], [NewParas,U{j}'], options);
        
        Sensitivities = reshape(simdata.statevalues(2,SensRange),NumOfSpecies,NumOfParameters);
        
        ySensitivities{j}=OutputMatrix*Sensitivities;
        % Get species time series
        DataTemp         = simdata.statevalues(2,1:NumOfSpecies)';
        XEstimates       = DataTemp;
        Y{j}=OutputMatrix*XEstimates;
        
        % Now calculate metric tensor
        G = G + ySensitivities{j}'*bsxfun(@rdivide,ySensitivities{j},Data{2,j}.^2);
        
    end
    
    G = G - diag(GetPriorLog2ndDeriv(PriorInfo, NewParas));
    % Calculate gradient
    GradL = zeros(1, NumOfParameters);
    for j=1:NumOfObs % NumOfObs is the number of performed experiments
        GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*ySensitivities{j};
    end
    for ParaNum  = 1:NumOfParameters
        GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, NewParas(ParaNum));
    end
    
    % Calculate the log prior for current hyperparameter value
    for a = 1:NumOfParameters
        ProposedLogPrior(a) = GetPrior(PriorInfo, a, NewParas(a));
    end
    
    if min(ProposedLogPrior) > 0 
        
        ProposedLogPrior = log(ProposedLogPrior);
        
        % calculate reverse transition probability
             ReverseMean = OriginalParas + 0.5*StepSize*StepSize*(G\GradL')';
        
        % Calculate proposed H value
        % Get species time series

        for j=1:NumOfObs
            LL(j) =  -0.5*sum(((Y{j}-Data{1,j})./Data{2,j}).^2);
        end


        ProposedLL = -0.5*(NoiseCovLogDet + NumOfObs*NumOfOutputs*log(2*pi)) + sum(LL);
        
    
        CholG          = chol(G);
        ProposedLogDet = 0.5*( log(2) + NumOfParameters*log(pi)) + sum(log(diag(CholG)));
        TransitionProbabilityNewGivenOld = -0.5*sum(((StepSize*CurrentCholG)\(NewParas - ProposalMean)').^2) - sum(log(diag(StepSize*CurrentCholG)));
        TransitionProbabilityOldGivenNew = -0.5*sum(((StepSize*CholG)\(OriginalParas - ReverseMean)').^2) - sum(log(diag(StepSize*CholG)));

        ProposedH      = -(ProposedLL + sum(ProposedLogPrior)) + ProposedLogDet;
        
        % Calculate current H value
        
       

        
        % Calculate the log prior for current hyperparameter value
        for a = 1:NumOfParameters
            %CurrentLogPrior(a) = ModelParameterPrior(a, Parameters(a));
            CurrentLogPrior(a) = log(GetPrior(PriorInfo, a, Parameters(a)));
        end
        
        CurrentLogDet = 0.5*( log(2) + NumOfParameters*log(pi)) + sum(log(diag(CurrentCholG))) ;
        CurrentH      = -(sum(CurrentLL) + sum(CurrentLogPrior)) + CurrentLogDet;
        
        % Accept according to ratio
        Ratio = -ProposedH + CurrentH + TransitionProbabilityOldGivenNew - TransitionProbabilityNewGivenOld;
        
        
        if Ratio > 0 || (Ratio > log(rand))
            % Accept proposal
            % Update variables
            Parameters                 = NewParas;
            CurrentLL                  = ProposedLL;
            AcceptedMutation           = AcceptedMutation + 1;
            
            CurrentG           = G;
            CurrentGradL       = GradL;
            CurrentCholG       = CholG;
            
            %disp('Accepted')
        else
            %disp('Rejected')
        
        end
    end
    
    %%%%%%%%%%%%%%%%%%%
    % Save parameters %
    %%%%%%%%%%%%%%%%%%%
    if Converged
        ParaHistory(IterationNum-ConvergenceIterationNum, :) = Parameters;
        LLHistory(IterationNum-ConvergenceIterationNum, :)   = CurrentLL;
        
        %if SaveMetricTensors
        %    MetricTensorHistory{IterationNum-ConvergenceIterationNum} = CurrentG;
        %end
        
    end
    
    
    
    % If not yet converged...
    if Converged == false
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust proposal widths %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Adjust parameter proposal widths
        if mod(IterationNum, MonitorRate) == 0
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['Iteration ' num2str(IterationNum)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(' ')
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust proposal width for parameter value inference %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if AttemptedMutation > 0
                disp([num2str(100*AcceptedMutation/AttemptedMutation) '% mutation acceptance']);
            end
            
            disp(' ')
            
            % Reset counters
            AttemptedMutation = 0;
            AcceptedMutation  = 0;
            
            
        end
        

        
        % Change converged tab if converged
        if IterationNum > MaxIterations && Converged == false
            Converged               = true;
            ConvergenceIterationNum = IterationNum;
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp(['Converged at iteration number ' num2str(ConvergenceIterationNum)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            
            BurnInTime = toc;
            tic;
            
        end
        
        
        
    else % Converged so decide how long to sample from posteriors
        waitbar((IterationNum- MaxIterations)/NumOfPosteriorSamples,WaitbarHandle);
        if IterationNum == ConvergenceIterationNum + NumOfPosteriorSamples
            % 5000 posterior samples have been collected so stop
            ContinueIterations = false;
        end
        
    end
    
    
end


PosteriorTime = toc;

CurTime = fix(clock);
RandTime = ceil(rand*10000);
c=clock;
% Save posterior
FileName = ['SB_SMMALA_' EquationName '_Obs' num2str(NumOfObs) ...
  '_date_' date '_time_' num2str(c(4)) '_' num2str(c(5))];
% 
% FileName = ['ODE_RMHMC_' EquationName '_' num2str(D) 'DPS_' num2str(floor(now)) '_' num2str(CurTime(4:6)) '_' num2str(RandTime)];
save(['./Results/' FileName], 'Options', 'ParaHistory', 'LLHistory', 'BurnInTime', 'PosteriorTime', 'Data');



end
