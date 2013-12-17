function [ output_args ] = PlotOutput( PosteriorSamples, TimePoints )

[NumOfSamples, NumOfParameters] = size(PosteriorSamples);

options.abstol = 1e-6;
options.reltol = 1e-6;

% Set up model
SBModelName           = 'SBFitzHughNagumo';
SBModelParameterNames = {'a' 'b' 'c'};

InitialValues = [-1 1];
%TimePoints    = 0:100/99:100;

figure(2)
hold on

for SampleNum = 1:NumOfSamples
    simdata    = feval(SBModelName, TimePoints, InitialValues, PosteriorSamples(SampleNum, :), options);
    t          = simdata.time;
    XEstimates = simdata.statevalues';

    plot(t, XEstimates)
    drawnow
    
end


end

