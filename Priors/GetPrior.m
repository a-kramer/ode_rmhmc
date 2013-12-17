function PP = GetPrior(PriorInfo, ParaNum, Value)

switch PriorInfo.Type{ParaNum}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uniform prior distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Para(1) = Lower bound
    % Para(2) = Upper bound
    case 'Uniform'
        LowerBound = PriorInfo.LowerBound(ParaNum);
        UpperBound = PriorInfo.UpperBound(ParaNum);
        
        if Value == 'random'
            % Produce random value from the prior
            PP = LowerBound + rand*(UpperBound-LowerBound);
        else
            % Calculate probability of value from the prior
            if (Value > UpperBound) || (Value < LowerBound)
                PP = 0;
            else
                PP = 1/(UpperBound-LowerBound);
            end
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normal prior distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Para(1) = Mean
    % Para(2) = SD
    case 'Normal'
        Para = PriorInfo.Para(ParaNum,:);
        
        if Value == 'random'
            % Produce random value from the prior
            PP = Para(1) + randn*Para(2);
        else
            % Calculate probability of value from the prior
            PP = normpdf(Value, Para(1), Para(2));
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gamma prior distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Para(1) = parameter 1
    % Para(2) = parameter 2
    case 'Gamma'
        Para = PriorInfo.Para(ParaNum,:);
        
        if Value == 'random'
            % Produce random value from the prior
            PP = gamrnd(Para(1), Para(2));
        else
            if (Value <= 0)
                PP = 0;
            else
                % Calculate probability of value from the prior
                PP = gampdf(Value, Para(1), Para(2));
            end
        end
        

end

        
        
end

