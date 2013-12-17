function [ output_args ] = GetSummary( )

Files      = dir('ODE*')
NumOfFiles = length(Files);

NumOfChains     = 3;
NumOfSamples    = 5000;
NumOfParameters = 3;
NumOfSpecies    = 2;

ESS       = zeros(NumOfChains*NumOfFiles, NumOfParameters);
NoiseESS  = zeros(NumOfChains*NumOfFiles, NumOfSpecies); % One noise term per species
MeanParas = zeros(NumOfChains*NumOfFiles, NumOfParameters);
VarParas  = zeros(NumOfChains*NumOfFiles, NumOfParameters);
StdParas  = zeros(NumOfChains*NumOfFiles, NumOfParameters);
VarNoiseParas  = zeros(NumOfChains*NumOfFiles, NumOfSpecies); % One noise term per species

for i = 1:NumOfFiles
    Data = open(Files(i).name);
    
    NumOfChains = length(Data.ParaHistory);
    
    for ChainNum = 1:NumOfChains
        MeanParas((i-1)*NumOfChains+ChainNum, :)      = mean(Data.ParaHistory{ChainNum});
        VarParas((i-1)*NumOfChains+ChainNum, :)       = var(Data.ParaHistory{ChainNum});
        StdParas((i-1)*NumOfChains+ChainNum, :)       = std(Data.ParaHistory{ChainNum});
        ESS((i-1)*NumOfChains+ChainNum, :)            = CalculateESS(Data.ParaHistory{ChainNum}, length(Data.ParaHistory{ChainNum})-1);
        NoiseESS((i-1)*NumOfChains+ChainNum, :)       = CalculateESS(Data.NoiseHistory{ChainNum}, length(Data.NoiseHistory{ChainNum})-1);
        VarNoiseParas((i-1)*NumOfChains+ChainNum, :)  = var(Data.NoiseHistory{ChainNum});
    end
    
    BurnInTime(i)    = Data.BurnInTime;
    PosteriorTime(i) = Data.PosteriorTime;
    
end

disp('%%%%% STATISTICAL ACCURACY RESULTS %%%%%')
disp(' ')
disp(['Mean Parameter Estimates: ' num2str(mean(MeanParas))])
disp(['S.D. in Parameter Estimates: ' num2str(std(MeanParas))])
disp(' ')
disp('%%%%% EFFICIENCY RESULTS %%%%%')
disp(' ')
disp(['Mean ESS for model parameters: ' num2str(mean(ESS))])
disp(['S.E. ESS for model parameters: ' num2str(std(ESS)/sqrt(NumOfChains*NumOfFiles))])
disp(['Mean marginal variance for model parameters: ' num2str(mean(VarParas))])
disp(' ')
disp(['Mean ESS for noise parameters: ' num2str(mean(NoiseESS))])
disp(['S.E. ESS for noise parameters: ' num2str(std(NoiseESS)/sqrt(NumOfChains*NumOfFiles))])
disp(['Mean marginal variance for noise parameters: ' num2str(mean(VarNoiseParas))])
disp(' ')
disp(['Average burn-in time taken per run: ' num2str(mean(BurnInTime))])
disp(['Average posterior time taken per run: ' num2str(mean(PosteriorTime))])

end

