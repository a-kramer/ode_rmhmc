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
ParaHistory         = zeros(NumOfPosteriorSamples,NumOfParameters);
LLHistory           = zeros(NumOfPosteriorSamples,1);
JrcHistory=zeros(1,NumOfPosteriorSamples);
GrcHistory=zeros(1,NumOfPosteriorSamples);

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
    dSf=cell(1,NumOfObs);
    dSh=cell(1,NumOfObs);
   [dSh{:}]=deal(zeros(F.l,F.m,F.m));

    for j=1:NumOfObs
         Jf{j}=F.Jf(xs{j},theta,U{j});
         Sf{j}=F.Sf(xs{j},theta,U{j});
         Sh{j}=F.C*Sf{j};         
         %fprintf('det(Jf{%i})=%g\t\tlog10(rcond(Jf{%i}))=%g\n',j,det(Jf{j}),j,log10(rcond(Jf{j})));
        dSf{j}=F.dSf(xs{j},theta,U{j},Jf{j},Sf{j});
        for k=1:F.m
         dSh{j}(:,:,k)=F.C*dSf{j}(:,:,k);
        end%for
    end
    
    % Calculate gradients for each of the parameters i.e. d(LL)/d(Parameter)
    GradL = zeros(1, NumOfParameters);
    % Data{1,·} contains the actual data for each experiment, while
    % Data{2,·} contains the respective measurement error (standard deviation)
    
    for j=1:NumOfObs
            GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*Sh{j};
    end
   
    % calculate GradLL using finite differences to be sure
    s=[+1,-1];
    h=1e-8; % finite difference scale
    GL=zeros(1,F.m);
    for i=1:2        % sign loop: +h and -h
        for k=1:F.m
            o=theta;
            o(k)=o(k)+s(i)*h;
            %calculate logLikelihood
            ll=0;
            for j=1:NumOfObs
                %Sh{j}
                CX=Y{j}+Sh{j}*(o-theta); %Sf{j} is verifiably correct
                %second method to find the new steady state (same results):
                %odeset('Jacobian',@(t,x) F.Jf(x,o,U{j}));
                %[T,X]=ode23s(@(t,x) F.f(x,o,U{j}),[0,1000],xs{j});
                %X=X(end,:)';
                %CX=F.C*X;                
                ll = ll - 0.5*sum(((CX-Data{1,j})./Data{2,j}).^2);                
            end
            %fprintf('[%i] ll(o(%i)%+i*h)=%.10f\n',i,k,s(i),ll);
            %calculate finite differences:
            GL(k)=GL(k)+s(i)*ll/(2*h);
        end
    end
    fprintf('relative difference (norm) between GradL and its finite difference approximation: %g\n',norm((GradL-GL)./GradL));    
    
    g={GradL,GL};
    
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
    
    
    
    % Now calculate the partial derivatives of the metric tensor
    CurrentGDeriv=cell(1,F.m);
    for k = 1:F.m 
        CurrentGDeriv{k} = zeros(NumOfParameters);
        for j=1:NumOfObs            
            CurrentGDeriv{k} = CurrentGDeriv{k} + ...
                dSh{j}(:,:,k)'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2) + ...
                bsxfun(@rdivide,Sh{j},Data{2,j}.^2)'*dSh{j}(:,:,k);
            % Add prior to the FI: - 3rd derivative of log prior
            CurrentGDeriv{k} = CurrentGDeriv{k} - diag(GetPriorLog3rdDeriv(PriorInfo, Parameters));
            
            CurrentInvGdG{k}       = CurrentG\CurrentGDeriv{k};
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
    CurrentLL = - 0.5*(NoiseCovLogDet + F.l*NumOfObs*log(2*pi)) + sum(LL)

    % make a small parameter step and verify GradL
    fprintf('some gradient(L) tests using finite difference approximations (h=%g):\n',h);
    %g{1}
    %g{2}
    for a=linspace(1,100,2)
        d=a*2e-6*g{2}'/(1e-4+norm(g{2}))
        o=theta+d;
        ll=0;
        for j=1:NumOfObs
            CX=Y{j}+Sh{j}*d; %Sf{j} is verifiably correct
            %fprintf('norm(F.f(xs{j},theta,U{j}))=%g\n',norm(F.f(xs{j},theta,U{j})));
            %fprintf('norm(F.f(X,o,U{j}))=%g\n',norm(F.f(X,o,U{j})));
            %CX=F.C*X;
            ll = ll - 0.5*sum(((CX-Data{1,j})./Data{2,j}).^2);
        end
        
        fprintf('(1) logLikelihood(theta+d)=%.10f\n',ll);
        fprintf('(2) logLikelihood(theta){%g} + GradL*d{%g} = %.10f\n',sum(LL),g{1}*d,sum(LL)+g{1}*d);
        fprintf('(3) logLikelihood(theta){%g} + GradL_h*d{%g} = %.10f\n',sum(LL),g{2}*d,sum(LL)+g{2}*d);
        fprintf('(2)-(1) = %g\n',sum(LL)+g{1}*d-ll);
        fprintf('(3)-(1) = %g\n',sum(LL)+g{2}*d-ll);
    end
%catch
    %
    % CurrentLL = -Inf;
%end

ProposedLogPrior=zeros(1,F.m);

disp('Initialisation Completed..');


WaitbarHandle=waitbar(0,sprintf('Sample Size: %i',NumOfPosteriorSamples));
%onCleanup(@() close(WaitbarHandle));
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
    
    
    
    ProposedMomentum = (randn(1,NumOfParameters)*CholG)';
    OriginalMomentum = ProposedMomentum;
    
    InvGdG      = CurrentInvGdG;
    TraceInvGdG = CurrentTraceInvGdG;
    
    TimeStep = 1;
    %if (rand > 0.5) TimeStep = 1; else TimeStep = -1; end
    
    %RandomSteps = ceil(rand*NumOfLeapFrogSteps);
    old_xs=xs;
    old_Jf=Jf;
    old_Sf=Sf;
    old_Sh=Sh;
    old_dSh=dSh;
    old_Y=Y;
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
                % fprintf('leapfrog[%i], FixedPoint[%i]: det(G)=%g\n',StepNum,FixedIter,det(G))
                InvGMomentum = G\PM;
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
            OriginalInvGMomentum  = G\ProposedMomentum;
            % fprintf('norm(OriginalInvGMomentum)=%g\n',norm(OriginalInvGMomentum));
            Pw = NewParas;
            theta_base=Pw';
            xs_base=xs; % for FixedIter loop
            Sf_base=Sf; % for FixedIter loop
            for FixedIter = 1:NumOfFixPointSteps
                wHist(FixedIter,:) = Pw;

                InvGMomentum = G\ProposedMomentum;
                % fprintf('norm(InvGMomentum)=%g\n',norm(InvGMomentum));
                ParameterIncrease = (TimeStep*(StepSize/2))*OriginalInvGMomentum' + (TimeStep*(StepSize/2))*InvGMomentum';
                if (isfinite(ParameterIncrease))
                    Pw = NewParas + ParameterIncrease;
                else
                    warning('ParameterIncrease is not finite: %g (norm)',norm(ParameterIncrease));
                    Pw = NewParas;
                end

                %%%%%%%%%%%%%%%
                % Calculate G %
                %%%%%%%%%%%%%%%
                
                for j=1:NumOfObs
                    
                    xs{j}=xs_base{j}+Sf_base{j}*ParameterIncrease'; % move steady states
                end
                
                theta=Pw'; % for steady state calculations
                xs=newton_raphson(F,xs,theta,U); 
                
                for j=1:NumOfObs                   
                    Y{j}         = F.C*xs{j};
                end
                
                RC=zeros(1,NumOfObs);
                for j=1:NumOfObs
                    Jf{j}=F.Jf(xs{j},theta,U{j});
                    RC(j)=rcond(Jf{j});
                    % fprintf('det(Jf{j=%i}=%g\nnorm(rho)=%g\nnorm(xs{j})=%g\nnorm(f(xs{j},rho,U{j})./xs{j})=%g\n',j,det(Jf{j}),norm(rho),norm(xs{j}),norm(F.f(xs{j},rho,U{j})./xs{j}));
                    if (any(isnan(Jf{j})))
                        Jf{j}
                        error('on main iteration %i, LeapFrogStep %i and FixPointStep %i the jacobian turned out to be NaN for input %i.\n',IterationNum,StepNum, FixedIter,j);
                    elseif (abs(det(Jf{j})) < 1e-8)
                        Jf{j}=Jf{j}+eye(F.n)*1e-8;
                        warning('on main iteration %i, LeapFrogStep %i and FixPointStep %i the jacobian had a very low determinant for input %i.\n',IterationNum,StepNum, FixedIter,j)
                    end
                    Sf{j}=F.Sf(xs{j},theta,U{j});
                    if (any(isnan(Sf{j})))
                        Sf{j}
                        error('on main iteration %i, LeapFrogStep %i and FixPointStep %i the sensitivity turned out to be NaN for input %i.\n',IterationNum,StepNum, FixedIter,j);
                    end
                    Sh{j}=F.C*Sf{j};
                end
                
                % Now calculate metric tensor
                G = zeros(NumOfParameters);
                for j=1:NumOfObs
                    G = G + Sh{j}'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2); % bsxfun does Sf{j}(:,k)./Data{2,j}(:) for all k
                end%for
                
                G = G - diag(GetPriorLog2ndDeriv(PriorInfo, Pw));


            end
            %fprintf('norm(G)=%g\n',norm(G));
            
            for j=1:NumOfObs
                dSf{j}=F.dSf(xs{j},theta,U{j},Jf{j},Sf{j});
                if (any(isnan(dSf{j})))
                    dSf{j}
                    error('on main iteration %i, LeapFrogStep %i and FixPointStep %i the second order sensitivity turned out to be NaN for input %i.\n',IterationNum,StepNum, FixedIter,j);
                end
                for k=1:F.m
                    dSh{j}(:,:,k)=F.C*dSf{j}(:,:,k);
                end%for
            end
            
 
            NewParas = Pw;


            % Get inverse metric tensor
            %InvG = inv(G + eye(F.m)*1e-8);


            % Calculate gradient and G
            GradL = zeros(1, F.m);
            for j=1:NumOfObs
                GradL = GradL - ((Y{j}-Data{1,j})./(Data{2,j}.^2))'*Sh{j};
            end
            for ParaNum  = 1:NumOfParameters
                GradL(ParaNum) = GradL(ParaNum) + GetPriorLogDeriv(PriorInfo, ParaNum, NewParas(ParaNum));
            end
            
            GDeriv=cell(1,F.m);
            for k = 1:F.m
                GDeriv{k} = zeros(NumOfParameters);
                for j=1:NumOfObs
                    GDeriv{k} = GDeriv{k} + ...
                        dSh{j}(:,:,k)'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2) + ...
                        bsxfun(@rdivide,Sh{j},Data{2,j}.^2)'*dSh{j}(:,:,k);
                    % Add prior to the FI: - 3rd derivative of log prior
                    GDeriv{k} = GDeriv{k} - diag(GetPriorLog3rdDeriv(PriorInfo, NewParas));
                    
                end%for
                %fprintf('main[%i], leapfrog[%i]: det(GDeriv{%i})=%g\n',IterationNum,StepNum,k,det(GDeriv{k}));
                InvGdG{k}       = G\GDeriv{k};
                TraceInvGdG(k)  = trace(InvGdG{k});
            end
            

            %%%%%%%%%%%%%%%%%%%
            % Update momentum %
            %%%%%%%%%%%%%%%%%%%
            
            % Calculate last term in dH/dTheta
            dJLdTheta = GradL;
            TraceTerm = 0.5*TraceInvGdG';

            InvGMomentum = G\ProposedMomentum;
            for d = 1:NumOfParameters
                LastTerm(d) = 0.5*((ProposedMomentum'*InvGdG{d}*InvGMomentum));
            end
            
            Move=TimeStep*(StepSize/2)*(dJLdTheta' - TraceTerm + LastTerm');
            %fprintf('distance moved: %g\n',norm(Move));
            ProposedMomentum = ProposedMomentum + Move;
            
        end

        
    catch Err
        Err.getReport()
        xs=old_xs;
        Jf=old_Jf;
        Sf=old_Sf;
        Sh=old_Sh;
        dSh=old_dSh;
        Y=old_Y;
        G = zeros(NumOfParameters);
        for j=1:NumOfObs
            G = G + Sh{j}'*bsxfun(@rdivide,Sh{j},Data{2,j}.^2); % bsxfun does Sf{j}(:,k)./Data{2,j}(:) for all k
        end%for
        G = G - diag(GetPriorLog2ndDeriv(PriorInfo, OriginalParas));
        GDeriv=cell(1,NumOfParameters);
        [GDeriv{:}]=deal(zeros(F.m));        
        IntegrationErr=true;
    end
    
        
    %fprintf('distance moved: %f\n',(norm(NewParas-Parameters)));    
    % Calculate the log prior for current hyperparameter value
    for a = 1:NumOfParameters
        ProposedLogPrior(a) = GetPrior(PriorInfo, a, NewParas(a));
    end
    
    %CurrentLL
    %ProposedLogPrior
   
    if (IntegrationErr==false)
        ProposedLogPrior = log(ProposedLogPrior);
        
        
        % Calculate proposed H value
        % Get species time series
        %XEstimates = simdata.statevalues(:,1:NumOfSpecies)';
    
        % Calculate the current likelihoods of the current parameters
        LL=zeros(1,NumOfObs);
        %theta=NewParas';
        %xs=newton_raphson(F,xs,theta,U); % correct steady states using newton raphson;
        
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
        ProposedH      = -(ProposedLL + sum(ProposedLogPrior)) + ProposedLogDet + (ProposedMomentum'*(G\ProposedMomentum))/2;
        
        
        % Calculate current H value
        % Calculate the log prior for current hyperparameter value
        for a = 1:NumOfParameters
            %CurrentLogPrior(a) = ModelParameterPrior(a, Parameters(a));
            CurrentLogPrior(a) = log(GetPrior(PriorInfo, a, Parameters(a)));
        end%for
        
        CurrentLogDet = 0.5*( log(2) + NumOfParameters*log(pi) + 2*sum(log(diag(CurrentCholG))) );
        CurrentH      = -(CurrentLL + sum(CurrentLogPrior)) + CurrentLogDet + (OriginalMomentum'*(CurrentG\OriginalMomentum))/2;
        
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
            CurrentCholG       = CholG;
            
            CurrentInvGdG      = InvGdG;
            CurrentTraceInvGdG = TraceInvGdG;
            
            %disp('Accepted')
        else
            %disp('Rejected')
            xs=old_xs;
            Jf=old_Jf;
            Sf=old_Sf;
            Sh=old_Sh;
            dSh=old_dSh;
            Y=old_Y;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%
    % Save parameters %
    %%%%%%%%%%%%%%%%%%%
    if Converged
        ParaHistory(IterationNum-ConvergenceIterationNum, :) = Parameters;
        LLHistory(IterationNum-ConvergenceIterationNum)   = CurrentLL;
        JrcHistory(IterationNum-ConvergenceIterationNum)=min(RC);
        GrcHistory(IterationNum-ConvergenceIterationNum)=rcond(G);
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
    
    
end%while


PosteriorTime = toc;

CurTime = fix(clock);
RandTime = ceil(rand*10000);

% Save posterior
c=clock;
FileName = ['ODE_RMHMC_' EquationName '_Obs' num2str(NumOfObs) ...
 '_LFsteps_' num2str(NumOfLeapFrogSteps) '_FPsteps_' num2str(NumOfFixPointSteps) ...
 '_date_' date '_time_' num2str(c(4)) '_' num2str(c(5))];
save(['./Results/' FileName], 'Options', 'ParaHistory', 'LLHistory', 'BurnInTime', 'PosteriorTime', 'Y', 'Data','F','JrcHistory','GrcHistory');



end
