function [] = SB_RMHMC_Single_MCMC(Data, U, Options )
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

%NoiseCov              = Options.NoiseCov;
Parameters            = Options.StartingParameters;
%InitialValues         = Options.StartingInitialConditions;

%InferICs              = Options.InferICs;

NumOfSens=NumOfSpecies*NumOfParameters;
NumOfdSens=NumOfSens*NumOfParameters;
SensRange=[NumOfSpecies+1:NumOfSpecies+NumOfSens];
dSensRange=[NumOfSpecies+NumOfSens+1:NumOfSpecies+NumOfSens+NumOfdSens];



% RMHMC Setup
NumOfLeapFrogSteps = Options.NumOfLeapFrogSteps;
StepSize           = Options.StepSize;
NumOfFixPointSteps = Options.NumOfFixPointSteps;


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
CurrentG     = zeros(NumOfParameters);
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

% Set up parameter history variable
ParaHistory         = [];
LLHistory           = [];

MetricTensorHistory = cell(1,NumOfPosteriorSamples);



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
    
    ySensitivities2=cell(1,NumOfObs);
    [ySensitivities2{:}]=deal(zeros(NumOfOutputs,NumOfParameters,NumOfParameters));
    
    
    Y=cell(1,NumOfObs);
    GradL = zeros(1, NumOfParameters);
    G = zeros(NumOfParameters);
    for j=1:NumOfObs
        % [InitialValues zeros(1,NumOfParameters*NumOfSpecies) ...
        % zeros(1,NumOfSpecies*NumOfParameters*NumOfParameters)]
        simdata = feval(SBModelName, [0,T], [], [Parameters, U{j}'], options);
        
        % Get sensitivities of all species with respect to parameter j
        Sensitivities = reshape(simdata.statevalues(2,SensRange),NumOfSpecies,NumOfParameters);
        % Get second order sensitivities of all species with respect to parameters j and k
        Sensitivities2 = reshape(simdata.statevalues(2,dSensRange),NumOfSpecies,NumOfParameters,NumOfParameters);
        
        ySensitivities{j}=OutputMatrix*Sensitivities;
        for k=1:NumOfParameters
         ySensitivities2{j}(:,:,k)=OutputMatrix*Sensitivities2(:,:,k);
        end
        % Get species time series
        DataTemp         = simdata.statevalues(2,1:NumOfSpecies)';
        XEstimates       = DataTemp;
        Y{j}=OutputMatrix*XEstimates;
        % Calculate gradients for each of the parameters i.e. d(LL)/d(Parameter)
        
        
        GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*ySensitivities{j};
        
        
        % Save current gradient
        CurrentGradL = GradL;
        
        % Now calculate metric tensor
           
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
    
    g={GradL,GL};
   
    
    
    
    
    for ParaNum  = 1:NumOfParameters
        GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, Parameters(ParaNum));
    end

    G = G - diag(GetPriorLog2ndDeriv(PriorInfo, Parameters));
    
    % Save current metric tensor
    CurrentG     = G;
    CurrentCholG = chol(G); 
    CurrentInvG  = inv(G + eye(NumOfParameters)*1e-8);
    
    
    % Now calculate the partial derivatives of the metric tensor
    CurrentGDeriv=cell(1,NumOfParameters);
    CurrentInvGdG=cell(1,NumOfParameters);
    CurrentTraceInvGdG=zeros(1,NumOfParameters);
    for k = 1:NumOfParameters
        CurrentGDeriv{k} = zeros(NumOfParameters);
        
        % Use sensitivities
        %
        % Standard derivatives of fisher expression for a gaussian
        for j=1:NumOfObs
        CurrentGDeriv{k} = CurrentGDeriv{k} + ...
                           ySensitivities2{j}(:,:,k)'*bsxfun(@rdivide,ySensitivities{j},Data{2,j}.^2) + ...
                           bsxfun(@rdivide,ySensitivities{j},Data{2,j}.^2)'*ySensitivities2{j}(:,:,k);
        end
        % Add prior to the FI: - 3rd derivative of log prior
        CurrentGDeriv{k} = CurrentGDeriv{k} - diag(GetPriorLog3rdDeriv(PriorInfo, Parameters));
        
        CurrentInvGdG{k}       = CurrentInvG*CurrentGDeriv{k};
        CurrentTraceInvGdG(k)  = trace(CurrentInvGdG{k});
        
    end
    
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
        d=a*2e-6*GL'/norm(GL)
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate parameter values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    OriginalParas = Parameters;
    NewParas      = OriginalParas;
    
    AttemptedMutation = AttemptedMutation + 1;
    
    
    
    GradL = CurrentGradL;
    G     = CurrentG;
    CholG = CurrentCholG;
    InvG  = CurrentInvG;
    
    
    ProposedMomentum = (randn(1,NumOfParameters)*CholG)';
    OriginalMomentum = ProposedMomentum;
    
    
    InvGdG      = CurrentInvGdG;
    TraceInvGdG = CurrentTraceInvGdG;
    
    if (randn > 0.5) TimeStep = 1; else TimeStep = -1; end
    
    %RandomSteps = ceil(rand*NumOfLeapFrogSteps);
    
    
    IntegrationErr = false;
    
    try
        
        % Perform leapfrog steps
        for StepNum = 1:NumOfLeapFrogSteps %RandomSteps

            %%%%%%%%%%%%%%%%%%%
            % Update momentum %
            %%%%%%%%%%%%%%%%%%%
            
            dJLdTheta = GradL;
            TraceTerm = 0.5*TraceInvGdG';

            % Multiple fixed point iteration
            PM = ProposedMomentum;
            for FixedIter = 1:NumOfFixPointSteps
                MomentumHist(FixedIter,:) = PM;

                InvGMomentum = InvG*PM;
                for d = 1:NumOfParameters
                    LastTerm(d)  = 0.5*(PM'*InvGdG{d}*InvGMomentum);
                end

                PM = ProposedMomentum + TimeStep*(StepSize/2)*(dJLdTheta' - TraceTerm + LastTerm');
            end
            ProposedMomentum = PM;




            %%%%%%%%%%%%%%%%%%%%%%%
            % Update w parameters %
            %%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Multiple Fixed Point Iteration %%%
            %
            %
            OriginalInvGMomentum  = G\ProposedMomentum;

            Pw = NewParas;
            for FixedIter = 1:NumOfFixPointSteps
                wHist(FixedIter,:) = Pw;

                InvGMomentum = G\ProposedMomentum;

                Pw = NewParas + (TimeStep*(StepSize/2))*OriginalInvGMomentum' + (TimeStep*(StepSize/2))*InvGMomentum';


                %%%%%%%%%%%%%%%
                % Calculate G %
                %%%%%%%%%%%%%%%
                %if FixedIter == NumOfFixPointSteps
                    % Calculate 2nd order sensitivities if last iteration
                G = zeros(NumOfParameters);
                for j=1:NumOfObs
%                                                                %[InitialValues ...
%                                                                %   zeros(1,NumOfParameters*NumOfSpecies) ...
%                                                                %   zeros(1,NumOfSpecies*NumOfParameters*NumOfParameters)]
                    simdata          = feval(SBModelName, [0,T], [], [Pw,U{j}'], options);

                    Sensitivities = reshape(simdata.statevalues(2,SensRange),NumOfSpecies,NumOfParameters);
                    Sensitivities2 = reshape(simdata.statevalues(2,dSensRange),NumOfSpecies,NumOfParameters,NumOfParameters);
                    
                    ySensitivities{j}=OutputMatrix*Sensitivities;
                    for k=1:NumOfParameters
                        ySensitivities2{j}(:,:,k)=OutputMatrix*Sensitivities2(:,:,k);
                    end
                    % Get species time series
                    DataTemp         = simdata.statevalues(2,1:NumOfSpecies)';
                    XEstimates       = DataTemp;
                    Y{j}=OutputMatrix*XEstimates;
                    
                    % Now calculate metric tensor
                    G = G + ySensitivities{j}'*bsxfun(@rdivide,ySensitivities{j},Data{2,j}.^2);
                    
                end

                G = G - diag(GetPriorLog2ndDeriv(PriorInfo, Pw));


            end
            NewParas = Pw;
          
            % Get inverse metric tensor
            InvG = inv(G + eye(NumOfParameters)*1e-8);


            % Calculate gradient and G
            GradL = zeros(1, NumOfParameters);
            for j=1:NumOfObs % NumOfObs is the number of performed experiments
                GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*ySensitivities{j};
            end
            for ParaNum  = 1:NumOfParameters
                GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, NewParas(ParaNum));
            end

            % Now calculate the partial derivatives of the metric tensor
            GDeriv=cell(1,NumOfParameters);
            for k = 1:NumOfParameters
                GDeriv{k} = zeros(NumOfParameters);

                % Use sensitivities
                %
                % Standard derivatives of fisher expression for a gaussian

                for j=1:NumOfObs
                    GDeriv{k} = GDeriv{k} + ...
                        ySensitivities2{j}(:,:,k)'*bsxfun(@rdivide,ySensitivities{j},Data{2,j}.^2) + ...
                        bsxfun(@rdivide,ySensitivities{j},Data{2,j}.^2)'*ySensitivities2{j}(:,:,k);
                end

                % Add prior to the FI: - 3rd derivative of log prior
                GDeriv{k} = GDeriv{k} - diag(GetPriorLog3rdDeriv(PriorInfo, NewParas));

                InvGdG{k}       = InvG*GDeriv{k};
                TraceInvGdG(k)  = trace(InvGdG{k});

            end



            %%%%%%%%%%%%%%%%%%%
            % Update momentum %
            %%%%%%%%%%%%%%%%%%%
            
            % Calculate last term in dH/dTheta
            dJLdTheta = GradL;
            TraceTerm = 0.5*TraceInvGdG';

            InvGMomentum = (InvG*ProposedMomentum);
            for d = 1:NumOfParameters
                LastTerm(d) = 0.5*((ProposedMomentum'*InvGdG{d}*InvGMomentum));
            end

            ProposedMomentum = ProposedMomentum + TimeStep*(StepSize/2)*(dJLdTheta' - TraceTerm + LastTerm');

        end
    catch    
       IntegrationErr = true;
       Err=lasterror();
       Err.message
    end
        
        
    % Calculate the log prior for current hyperparameter value
    for a = 1:NumOfParameters
        ProposedLogPrior(a) = GetPrior(PriorInfo, a, NewParas(a));
    end
    
    if min(ProposedLogPrior) > 0 && IntegrationErr == false
        
        ProposedLogPrior = log(ProposedLogPrior);
        
        
        % Calculate proposed H value
        % Get species time series

        for j=1:NumOfObs
            LL(j) =  -0.5*sum(((Y{j}-Data{1,j})./Data{2,j}).^2);
        end


        ProposedLL = -0.5*(NoiseCovLogDet + NumOfObs*NumOfOutputs*log(2*pi)) + sum(LL);
        
    
        CholG          = chol(G);
        ProposedLogDet = 0.5*( log(2) + NumOfParameters*log(pi) + 2*sum(log(diag(CholG))) );
        ProposedH      = -(ProposedLL + sum(ProposedLogPrior)) + ProposedLogDet + (ProposedMomentum'*InvG*ProposedMomentum)/2;
        
        
        % Calculate current H value
        % Calculate the log prior for current hyperparameter value
        for a = 1:NumOfParameters
            %CurrentLogPrior(a) = ModelParameterPrior(a, Parameters(a));
            CurrentLogPrior(a) = log(GetPrior(PriorInfo, a, Parameters(a)));
        end
        
        CurrentLogDet = 0.5*( log(2) + NumOfParameters*log(pi) + 2*sum(log(diag(CurrentCholG))) );
        CurrentH      = -(sum(CurrentLL) + sum(CurrentLogPrior)) + CurrentLogDet + (OriginalMomentum'*CurrentInvG*OriginalMomentum)/2;
        
        % Accept according to ratio
        Ratio = -ProposedH + CurrentH;
        
        
        if Ratio > 0 || (Ratio > log(rand))
            % Accept proposal
            % Update variables
            Parameters                 = NewParas;
            CurrentLL                  = ProposedLL;
            AcceptedMutation           = AcceptedMutation + 1;
            
            CurrentG           = G;
            CurrentGradL       = GradL;
            CurrentInvG        = InvG;
            CurrentCholG       = CholG;
            
            CurrentInvGdG      = InvGdG;
            CurrentTraceInvGdG = TraceInvGdG;
            
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
FileName = ['ODE_RMHMC_' EquationName '_Obs' num2str(NumOfObs) ...
 '_LFsteps_' num2str(NumOfLeapFrogSteps) '_FPsteps_' num2str(NumOfFixPointSteps) ...
 '_date_' date '_time_' num2str(c(4)) '_' num2str(c(5))];
% 
% FileName = ['ODE_RMHMC_' EquationName '_' num2str(D) 'DPS_' num2str(floor(now)) '_' num2str(CurTime(4:6)) '_' num2str(RandTime)];
save(['./Results/' FileName], 'Options', 'ParaHistory', 'LLHistory', 'BurnInTime', 'PosteriorTime', 'Data');



end
