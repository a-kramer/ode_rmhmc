function [] = steadystate_RMHMC_Single_MCMC( Data, U, F, Options )
%
% Y will have to be a cell array with one entry per experiment condition 
% (i.e. model input u)
% for the apoptosis model, the input is the initial concentration of one of
% the compounds.
% TimePoints serve no purpose here, since we are doing steady state
% analysis.
% So, instead, we pass a set of initial conditions U to the sampler
% Options will have to contain the ode-rhs function, its Jacobian and
% Sensitivities.



%warning off

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
NumOfSpecies = F.n;
NumOfObs     = length(U);

% Get specified options
NumOfParameters       = F.m;

% SBModelName           = Options.SBModelName;

% SpeciesObserved       = Options.ObservedSpecies;
% SpeciesUnobserved     = Options.UnobservedSpecies;

% N = length(SpeciesObserved) + length(SpeciesUnobserved);

MaxIterations         = Options.MaxIterations;
NumOfPosteriorSamples = Options.NumOfPosteriorSamples;

PriorInfo             = Options.PriorInfo;

EquationName          = Options.EquationName;

%NoiseCov              = Options.NoiseCov;
Parameters            = Options.StartingParameters;
%InitialValues         = Options.StartingInitialConditions;

%InferICs              = Options.InferICs;


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
    %for PopulationNum = 1:NumOfPopulations
        for n = 1:NumOfParameters
            Parameters(n) = GetPrior(PriorInfo, n, 'random');
        end
    %end
end


% Current Gradient
CurrentG     = zeros(NumOfParameters);
CurrentGradL = zeros(1,NumOfParameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup initial values for solving the ODE % not needed anymore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isempty(InitialValues)
%     %for PopulationNum = 1:NumOfPopulations
%         for n = 1:N
%             InitialValues(n) = GetPrior(PriorInfoIC, n, 'random');
%         end
%     %end
% end
% 

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
%for i = 1:F.l
%    NoiseCovInv{i}    = inv(NoiseCov{i} + eye(D)*1e-8);
    NoiseCovLogDet = 2*sum(log(cat(1,Data{2,:})));
%end



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
    rho=exp(reshape(Parameters,[F.m,1]));
    for j=1:NumOfObs
        odeset('Jacobian',@(t,x) F.Jf(x,rho,U{j}));
        if (isreal(F.x0))
            x0=F.x0;
        else
            x0=F.x0(rho,U{j});
        end%if
        [T,X]=ode23s(@(t,x) F.f(x,rho,U{j}),[0,30000],x0);   % later on we do newton_raphson(F,U,rho,U);
        figure(j);
        plot(T,X);
        xs{j}=X(end,:)';
%         fprintf('xs{j}\n');
%         xs{j}
        fprintf('log(f{j})\n');
        real(log(F.f(xs{j},rho,U{j})))
    end
    % simdata = feval(SBModelName, U, [InitialValues zeros(1,NumOfParameters*NumOfSpecies) zeros(1,sum(1:NumOfParameters)*NumOfSpecies)], Parameters, options);

    
    % Get species time series
    Y=cell(1,NumOfObs);
    for j=1:NumOfObs
     Y{j}         = F.C*xs{j};
    end

    % Get sensitivities of all species with respect to Observation j
     Sf=cell(1,NumOfObs);
     Jf=cell(1,NumOfObs);
    dSf=cell(1,NumOfObs);
    dSh=cell(1,NumOfObs);
   [dSh{:}]=deal(zeros(F.l,F.m,F.m));

    for j=1:NumOfObs
         Sf{j}=F.Sf(xs{j},rho,U{j});
         Sh{j}=F.C*Sf{j}*diag(rho);
         Jf{j}=F.Jf(xs{j},rho,U{j});
         fprintf('det(Jf{%i})=%f\n',j, det(Jf{j}));
        dSf{j}=F.dSf(xs{j},rho,U{j},Jf{j},Sf{j});
        for k=1:F.m % here I use the sparse function because it is convenient to build the apropriate matrix to add;
                                                                  % all rows, one column: k
         dSh{j}(:,:,k)=F.C*dSf{j}(:,:,k)*diag(rho)*rho(k) + sparse(1:F.l,k,Sh{j}(:,k)*rho(k),F.l,F.m);
        end%for
    end
%     Sf
%     Sh
%     Jf
%     dSf
%     dSh
%     Y
    %pause();
%     for j = 1:NumOfParameters
%         Sensitivities{j} = simdata.statevalues(:,(NumOfSpecies*j)+1:(NumOfSpecies*(j+1)));
%         
%         % Get second order sensitivities of all species with respect to parameters j and k
%         for k = j:NumOfParameters
%             CurrentStartIndex   = ( (NumOfSpecies*(NumOfParameters+1)) + (sum(1:NumOfParameters)-sum(1:NumOfParameters-(j-1)))*NumOfSpecies ) + (k-j)*NumOfSpecies + 1;
%             Sensitivities2{j,k} = simdata.statevalues(:, CurrentStartIndex:(CurrentStartIndex+(NumOfSpecies-1)));
%             Sensitivities2{k,j} = Sensitivities2{j,k};
%         end
%         
%     end

    
    % Calculate gradients for each of the parameters i.e. d(LL)/d(Parameter)
    GradL = zeros(1, NumOfParameters);
    % Data{1,·} contains the actual data for each experiment, while
    % Data{2,·} contains the respective measurement error (standard deviation)
    
    for j=NumOfObs
            GradL = GradL - 0.5*((Y{j}-Data{1,j})./(Data{2,j}.^2))'*Sh{j};
    end
    for ParaNum  = 1:NumOfParameters
        GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, Parameters(ParaNum));
    end
    
    % Save current gradient
    CurrentGradL = GradL;
    
    % Now calculate metric tensor
    G = zeros(NumOfParameters);
    for j=1:NumOfObs
     G = G + Sh{j}'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2); % bsxfun does Sf(:,k)./Data{2,j}(:) for all k
    end%for
    
    G = G - diag(GetPriorLog2ndDeriv(PriorInfo, Parameters));
    
    % Save current metric tensor
    
    % fprintf('det(G)=%f\n',det(G));
    CurrentG     = G;
    CurrentCholG = chol(G); 
    CurrentInvG  = inv(G + eye(NumOfParameters)*1e-8);
    
    
    % Now calculate the partial derivatives of the metric tensor
    CurrentGDeriv=cell(1,F.m);
    for k = 1:F.m 
        CurrentGDeriv{k} = zeros(NumOfParameters);
        for j=1:NumOfObs            
            CurrentGDeriv{k} = CurrentGDeriv{k} + ...
                dSh{j}(:,:,k)'*bsxfun(@rdivide,Sh{j},Data{2,j}) + ...
                bsxfun(@rdivide,Sh{j},Data{2,j})'*dSh{j}(:,:,k);
            % Add prior to the FI: - 3rd derivative of log prior
            CurrentGDeriv{k} = CurrentGDeriv{k} - diag(GetPriorLog3rdDeriv(PriorInfo, Parameters));
            
            CurrentInvGdG{k}       = CurrentInvG*CurrentGDeriv{k};
            CurrentTraceInvGdG(k)  = trace(CurrentInvGdG{k});
        end%for
        % fprintf('CurrentGDeriv{%i}\n',k);
        % CurrentGDeriv{k}
    end
    
    
    
    % Calculate the current likelihoods of the current parameters
    LL=zeros(1,NumOfObs);
    for j=1:NumOfObs
        %YDdifference_j=((Y{j}-Data{1,j})./Data{2,j}).^2
        LL(j) = - 0.5*sum(((Y{j}-Data{1,j})./Data{2,j}).^2);
    end
    %LL
    CurrentLL = - 0.5*(NoiseCovLogDet + F.l*NumOfObs*log(2*pi)) +sum(LL);

%catch
    %
    % CurrentLL = -Inf;
%end



disp('Initialisation Completed..');


WaitbarHandle=waitbar(0,sprintf('Sample Size: %i',NumOfPosteriorSamples));

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
    
    RandomSteps = ceil(rand*NumOfLeapFrogSteps);
    
    
    % IntegrationErr = false;
    
    
        
        % Perform leapfrog steps
        for StepNum = 1:RandomSteps

            %%%%%%%%%%%%%%%%%%%
            % Update momentum %
            %%%%%%%%%%%%%%%%%%%
            
            dJLdTheta = GradL;
            TraceTerm = 0.5*TraceInvGdG';

            % Multiple fixed point iteration
            PM = ProposedMomentum;
            for FixedIter = 1:NumOfFixPointSteps
                MomentumHist(FixedIter,:) = PM;
                % fprintf('leapfrog[%i], FixedPoint[%i]: det(G)=%g\n',StepNum,FixedIter,det(G))
                InvGMomentum = (G+eye(F.m)*1e-3)\PM;
                for d = 1:NumOfParameters
                    LastTerm(d)  = 0.5*(PM'*InvGdG{d}*InvGMomentum);
                end
                % LastTerm
                PM = ProposedMomentum + TimeStep*(StepSize/2)*(dJLdTheta' - TraceTerm + LastTerm');
            end
            ProposedMomentum = PM;




            %%%%%%%%%%%%%%%%%%%%%%%
            % Update w parameters %
            %%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Multiple Fixed Point Iteration %%%
            %
            %
            % ProposedMomentum
            % G
            % det(G)
            OriginalInvGMomentum  = (G + eye(NumOfParameters)*1e-3)\ProposedMomentum;
            % fprintf('norm(OriginalInvGMomentum)=%g\n',norm(OriginalInvGMomentum));
            Pw = NewParas;
            rho_base=exp(reshape(Pw,[F.m,1]));
            xs_base=xs;
            Sf_base=Sf;
            for FixedIter = 1:NumOfFixPointSteps
                wHist(FixedIter,:) = Pw;

                InvGMomentum = (G + eye(NumOfParameters)*1e-3)\ProposedMomentum;
                % fprintf('norm(InvGMomentum)=%g\n',norm(InvGMomentum));
                ParameterIncrease=(TimeStep*(StepSize/2))*OriginalInvGMomentum' + (TimeStep*(StepSize/2))*InvGMomentum';
 
                Pw = NewParas + ParameterIncrease;


                %%%%%%%%%%%%%%%
                % Calculate G %
                %%%%%%%%%%%%%%%
                
                for j=1:NumOfObs
                    xs{j}=xs_base{j}+Sf_base{j}*diag(rho_base)*reshape(ParameterIncrease,[F.m,1]); % move steady states
                end
                rho=exp(reshape(Pw,[F.m,1]));
                xs=newton_raphson(F,xs,rho,U); % correct steady states using newton raphson;  
                
                for j=1:NumOfObs                   
                    Y{j}         = F.C*xs{j};
                end
                
                for j=1:NumOfObs
                    Jf{j}=F.Jf(xs{j},rho,U{j});
                    % fprintf('det(Jf{j=%i}=%g\nnorm(rho)=%g\nnorm(xs{j})=%g\nnorm(f(xs{j},rho,U{j})./xs{j})=%g\n',j,det(Jf{j}),norm(rho),norm(xs{j}),norm(F.f(xs{j},rho,U{j})./xs{j}));
                    if (any(isnan(Jf{j})))
                        Jf{j}
                        error(sprintf('on main iteration %i, LeapFrogStep %i and FixPointStep %i the jacobian turned out to be NaN for input %i.\n',IterationNum,StepNum, FixedIter,j));
                    elseif (abs(det(Jf{j})) < 1e-4)
                        Jf{j}=Jf{j}+eye(F.n)*1e-2;
                        warning(sprintf('on main iteration %i, LeapFrogStep %i and FixPointStep %i the jacobian had a very low determinant for input %i.\n',IterationNum,StepNum, FixedIter,j))
                    end
                    Sf{j}=F.Sf(xs{j},rho,U{j});
                    if (any(isnan(Sf{j})))
                        Sf{j}
                        error(sprintf('on main iteration %i, LeapFrogStep %i and FixPointStep %i the sensitivity turned out to be NaN for input %i.\n',IterationNum,StepNum, FixedIter,j));
                    end
                    Sh{j}=F.C*Sf{j}*diag(rho);

                end
                % Sf
                % Sh
                % Jf
                % dSf
                % dSh
                % pause();
                % Now calculate metric tensor
                G = zeros(NumOfParameters);
                for j=1:NumOfObs
                    G = G + Sh{j}'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2); % bsxfun does Sf{j}(:,k)./Data{2,j}(:) for all k
                end%for
                
                G = G - diag(GetPriorLog2ndDeriv(PriorInfo, Pw));


            end
            dSf{j}=F.dSf(xs{j},rho,U{j},Jf{j},Sf{j});
            if (any(isnan(dSf{j})))
                dSf{j}
                error(sprintf('on main iteration %i, LeapFrogStep %i and FixPointStep %i the second order sensitivity turned out to be NaN for input %i.\n',IterationNum,StepNum, FixedIter,j));
            end
            
            for k=1:F.m % here I use the sparse function because it is convenient to build the apropriate matrix to add;
                dSh{j}(:,:,k)=F.C*dSf{j}(:,:,k)*diag(rho)*rho(k) + sparse(1:F.l,k,Sh{j}(:,k)*rho(k),F.l,F.m);
            end%for

            NewParas = Pw;


            % Get inverse metric tensor
            %InvG = inv(G + eye(F.m)*1e-8);


            % Calculate gradient and G
            GradL = zeros(1, F.m);
            for j=NumOfObs
                GradL = GradL - 0.5*((Y{j}-Data{1,j})./(Data{2,j}.^2))'*Sh{j};
            end
            for ParaNum  = 1:NumOfParameters
                GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, NewParas(ParaNum));
            end
            
            GDeriv=cell(1,F.m);
            for k = 1:F.m
                GDeriv{k} = zeros(NumOfParameters);
                for j=1:NumOfObs
                    GDeriv{k} = GDeriv{k} + ...
                        dSh{j}(:,:,k)'*bsxfun(@rdivide,Sh{j},Data{2,j}) + ...
                        bsxfun(@rdivide,Sh{j},Data{2,j})'*dSh{j}(:,:,k);
                    % Add prior to the FI: - 3rd derivative of log prior
                    GDeriv{k} = GDeriv{k} - diag(GetPriorLog3rdDeriv(PriorInfo, NewParas));
                    
                end%for
                %fprintf('main[%i], leapfrog[%i]: det(GDeriv{%i})=%g\n',IterationNum,StepNum,k,det(GDeriv{k}));
                InvGdG{k}       = (G+eye(F.m)*1e-6)\GDeriv{k};
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

        
    %fprintf('distance moved: %f\n',(norm(NewParas-Parameters)));    
    % Calculate the log prior for current hyperparameter value
    for a = 1:NumOfParameters
        ProposedLogPrior(a) = GetPrior(PriorInfo, a, NewParas(a));
    end
    
    if min(ProposedLogPrior) > 0 
        ProposedLogPrior = log(ProposedLogPrior);
        
        
        % Calculate proposed H value
        % Get species time series
        %XEstimates = simdata.statevalues(:,1:NumOfSpecies)';
    
        % Calculate the current likelihoods of the current parameters
        LL=zeros(1,NumOfObs);
        rho=exp(reshape(NewParas,[F.m,1]));
        xs=newton_raphson(F,xs,rho,U); % correct steady states using newton raphson;
        
        for j=1:NumOfObs
            %fprintf('acceptance step: norm(xs{%i})=%g\n',j,norm(xs{j}));
            Y{j}         = F.C*xs{j};
        end
        
        for j=1:NumOfObs
            %YDdifference_j=((Y{j}-Data{1,j})./Data{2,j}).^2;
            LL(j) =  -0.5*sum(((Y{j}-Data{1,j})./Data{2,j}).^2);
        end

        ProposedLL = -0.5*(NoiseCovLogDet + NumOfObs*F.l*log(2*pi)) + sum(LL);
        
    
        CholG          = chol(G);
        ProposedLogDet = 0.5*( log(2) + NumOfParameters*log(pi) + 2*sum(log(diag(CholG))) );
        %ProposedH      = -(sum(ProposedLL) + sum(ProposedLogPrior)) + ProposedLogDet + (ProposedMomentum'*InvG*ProposedMomentum)/2;
        ProposedH      = -(ProposedLL + sum(ProposedLogPrior)) + ProposedLogDet + (ProposedMomentum'*InvG*ProposedMomentum)/2;
        
        
        % Calculate current H value
        % Calculate the log prior for current hyperparameter value
        for a = 1:NumOfParameters
            %CurrentLogPrior(a) = ModelParameterPrior(a, Parameters(a));
            CurrentLogPrior(a) = log(GetPrior(PriorInfo, a, Parameters(a)));
        end%for
        
        CurrentLogDet = 0.5*( log(2) + NumOfParameters*log(pi) + 2*sum(log(diag(CurrentCholG))) );
        %CurrentH      = -(sum(CurrentLL) + sum(CurrentLogPrior)) + CurrentLogDet + (OriginalMomentum'*CurrentInvG*OriginalMomentum)/2;
        CurrentH      = -(CurrentLL + sum(CurrentLogPrior)) + CurrentLogDet + (OriginalMomentum'*CurrentInvG*OriginalMomentum)/2;
        
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

% Save posterior
FileName = ['ODE_RMHMC_' EquationName '_' num2str(NumOfObs) 'DPS_' num2str(floor(now)) '_' num2str(CurTime(4:6)) '_' num2str(RandTime)];
save(['./Results/' FileName], 'ParaHistory', 'LLHistory', 'BurnInTime', 'PosteriorTime', 'Y', 'Data','F');



end
