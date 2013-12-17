function [ ] = CalculateStatistics( Method, DataSet )

% Get stats
Files = dir(['*' Method '*' DataSet '*.mat']);


for i = 1:length(Files)

    Data     = open(Files(i).name);
    
    ESS{i}   = CalculateESS(Data.betaPosterior, length(Data.betaPosterior)-1);
    
    MinESS(i)    = min(ESS{i});
    MaxESS(i)    = max(ESS{i});
    MedianESS(i) = median(ESS{i});
    MeanESS(i)   = mean(ESS{i});
    Times(i)     = Data.TimeTaken;
    
end

disp(['Results for ' Method ' with ' DataSet ' dataset.'])
disp(' ')

disp(['Min:    ' num2str(mean(MinESS)) ' +/- ' num2str(std(MinESS)/sqrt(length(MinESS)))])
disp(['Median: ' num2str(mean(MedianESS)) ' +/- ' num2str(std(MedianESS)/sqrt(length(MedianESS)))])
disp(['Mean:   ' num2str(mean(MeanESS)) ' +/- ' num2str(std(MeanESS)/sqrt(length(MeanESS)))])
disp(['Max:    ' num2str(mean(MaxESS)) ' +/- ' num2str(std(MaxESS)/sqrt(length(MaxESS)))])
disp(['Time:   ' num2str(mean(Times)) ' +/- ' num2str(std(Times)/sqrt(length(Times)))])

disp('')
disp(['Time per Min ESS: ' num2str(mean(Times)/mean(MinESS))])

end
