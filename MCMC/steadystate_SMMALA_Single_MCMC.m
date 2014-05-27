function [] = steadystate_SMMALA_Single_MCMC( Data, U, F, Options )
% Copyright 2012-2013 Ben Calderhead <b.calderhead@imperial.ac.uk>, Andrei Kramer <andrei.kramer@ist.uni-stuttgart.de>

%warning off
Options
% Start the timer
tic

%rand('twister', sum(100*clock))
%randn('state', sum(100*clock))
rand('twister', 1)
randn('state', 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise user options if not already specified                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get N (number of chemical species) and D (number of time points)
[N, D] = size(Data);
NumOfObs     = D;

% Get specified options
NumOfSpecies          = F.n;
NumOfParameters       = F.m;
NumOfOutputs          = F.l;

MaxIterations         = Options.MaxIterations;
NumOfPosteriorSamples = Options.NumOfPosteriorSamples;

PriorInfo             = Options.PriorInfo;

EquationName          = Options.EquationName;

Parameters            = Options.StartingParameters;
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


% Current Gradient
CurrentGradL = zeros(1,NumOfParameters);


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

    theta=Parameters';
    xs=cell(1,NumOfObs);
    for j=1:NumOfObs
        odeset('Jacobian',@(t,x) F.Jf(x,theta,U{j}));
        if isreal(F.x0)
           x0=F.x0;
	else
           x0=F.x0(U{j});
        end
        [T,X]=ode23s(@(t,x) F.f(x,theta,U{j}),[0,1000],x0);   % later on we do newton_raphson
        figure(j);
        plot(T,X);
        xs{j}=X(end,:)';
        % fprintf('xs{j}\n');
        % xs{j}
        % fprintf('log10(f{%i}):\n',j);
        % real(log10(F.f(xs{j},theta,U{j})))
    end
    % simdata = feval(SBModelName, U, [InitialValues zeros(1,NumOfParameters*NumOfSpecies) zeros(1,sum(1:NumOfParameters)*NumOfSpecies)], Parameters, options);

    % model outputs
    Y=cell(1,NumOfObs);
    for j=1:NumOfObs
     Y{j}         = F.C*xs{j};
    end

    % Get sensitivities of all species with respect to Observation j
     Sf=cell(1,NumOfObs);
     Sh=cell(1,NumOfObs);
     Jf=cell(1,NumOfObs);

    for j=1:NumOfObs
         Jf{j}=F.Jf(xs{j},theta,U{j});
         Sf{j}=F.Sf(xs{j},theta,U{j});
         Sh{j}=F.C*Sf{j};         
    end
    
    % Calculate gradients for each of the parameters i.e. d(LL)/d(Parameter)
    GradL = zeros(1, NumOfParameters);
    % Data{1,·} contains the actual data for each experiment, while
    % Data{2,·} contains the respective measurement error (standard deviation)
    
    for j=1:NumOfObs
            GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*Sh{j};
    end
    for ParaNum  = 1:NumOfParameters
        GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, Parameters(ParaNum));
    end
    CurrentGradL=GradL;
    
    % Now calculate metric tensor
    G = zeros(NumOfParameters);
    for j=1:NumOfObs
     G = G + Sh{j}'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2); % bsxfun does Sf(:,k)./Data{2,j}(:) for all k
    end%for
    G = G - diag(GetPriorLog2ndDeriv(PriorInfo, Parameters));
    
    % Save current metric tensor
    CurrentG     = G;
    CurrentCholG = chol(G); 
    
    LL=zeros(1,NumOfObs);
    for j=1:NumOfObs
        %YDdifference_j=((Y{j}-Data{1,j})./Data{2,j}).^2
        LL(j) = - 0.5*sum(((Y{j}-Data{1,j})./Data{2,j}).^2);
    end
    %LL
    CurrentLL = - 0.5*(NoiseCovLogDet + NumOfOutputs*NumOfObs*log(2*pi)) + sum(LL);
    
%catch
    
%     for n=1:SpeciesObserved
%         CurrentLL(n) = -Inf;
%     end
    
%end
ProposedLogPrior=zeros(1,NumOfParameters);
CurrentLogPrior=zeros(1,NumOfParameters);
for a = 1:NumOfParameters
     CurentLogPrior(a) = GetPrior(PriorInfo, a, Parameters(a));
end

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
ParaHistory=zeros(MaxIterations,F.m);
LLHistory=zeros(MaxIterations,1);

while ContinueIterations
    
    % Increment iteration number
    IterationNum = IterationNum + 1;
    
    OriginalParas = Parameters;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate parameter values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    AttemptedMutation = AttemptedMutation + 1;

    % parameter update
             ProposalMean = OriginalParas + 0.5*StepSize*StepSize*(CurrentG\CurrentGradL')';
        
    NewParas = ProposalMean + StepSize*(CurrentCholG\randn(NumOfParameters,1))';

    for j=1:NumOfObs
        xs{j}=xs{j}+Sf{j}*(NewParas - OriginalParas)'; % move steady states
    end
    
    theta=NewParas'; % for steady state calculations
    xs=newton_raphson(F,xs,theta,U);
    
    for j=1:NumOfObs
        Y{j}         = F.C*xs{j};
    end
    
    for j=1:NumOfObs
        Jf{j}=F.Jf(xs{j},theta,U{j});
        % fprintf('det(Jf{j=%i}=%g\nnorm(rho)=%g\nnorm(xs{j})=%g\nnorm(f(xs{j},rho,U{j})./xs{j})=%g\n',j,det(Jf{j}),norm(rho),norm(xs{j}),norm(F.f(xs{j},rho,U{j})./xs{j}));
        if (any(isnan(Jf{j})))
            Jf{j}
            error('on main iteration %i the jacobian turned out to be NaN for input %i.\n',IterationNum,j);
        end
        Sf{j}=F.Sf(xs{j},theta,U{j});
        if (any(isnan(Sf{j})))
            Sf{j}
            error('on main iteration %i the sensitivity turned out to be NaN for input %i.\n',IterationNum,j);
        end
        Sh{j}=F.C*Sf{j};
    end
   
    % Calculate gradient and G
    GradL = zeros(1, F.m);
    for j=1:NumOfObs
        GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*Sh{j};
    end
    for ParaNum  = 1:NumOfParameters
        GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, NewParas(ParaNum));
    end   
    
    % Now calculate metric tensor
    G = zeros(NumOfParameters);
    for j=1:NumOfObs
     G = G + Sh{j}'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2); % bsxfun does Sf(:,k)./Data{2,j}(:) for all k
    end%for
    G = G - diag(GetPriorLog2ndDeriv(PriorInfo, NewParas));
    
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

        ProposedH      = -(ProposedLL + sum(ProposedLogPrior)) + ProposedLogDet ;
        
        % Calculate current H value

        
        % Calculate the log prior for current hyperparameter value
        for a = 1:NumOfParameters
            %CurrentLogPrior(a) = ModelParameterPrior(a, Parameters(a));
            CurrentLogPrior(a) = log(GetPrior(PriorInfo, a, OriginalParas(a)));
        end
        
        CurrentLogDet = 0.5*( log(2) + NumOfParameters*log(pi)) + sum(log(diag(CurrentCholG)));
        CurrentH      = -(CurrentLL + sum(CurrentLogPrior)) + CurrentLogDet;
        
        % Accept according to ratio
        Ratio = -ProposedH + CurrentH + TransitionProbabilityOldGivenNew - TransitionProbabilityNewGivenOld;
        
        
        if Ratio > 0 || (Ratio > log(rand))
            % Accept proposal
            % Update variables
            Parameters                 = NewParas;
            CurrentLL                  = ProposedLL;
            CurrentG           = G;
            CurrentGradL       = GradL;
            CurrentCholG       = CholG;
            
            AcceptedMutation           = AcceptedMutation + 1;
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
FileName = ['steadystate_SMMALA_' EquationName '_Obs' num2str(NumOfObs) ...
  '_date_' date '_time_' num2str(c(4)) '_' num2str(c(5))];
% 
% FileName = ['ODE_RMHMC_' EquationName '_' num2str(D) 'DPS_' num2str(floor(now)) '_' num2str(CurTime(4:6)) '_' num2str(RandTime)];
save(['./Results/' FileName], 'Options', 'ParaHistory', 'LLHistory', 'BurnInTime', 'PosteriorTime', 'Data');



end
