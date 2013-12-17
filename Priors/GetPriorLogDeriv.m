function PP = GetPriorLogDeriv(PriorInfo, ParaNum, Value)

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
        
        % Calculate derivative of log prior
        if (Value > UpperBound) || (Value < LowerBound)
            PP = -1e300;
        else
            PP = 0;
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normal prior distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Para(1) = Mean
    % Para(2) = SD
    case 'Normal'
        Para = PriorInfo.Para(ParaNum,:);
        
        % Calculate derivative of log prior
        PP = -(Value - Para(1))/Para(2)^2;
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gamma prior distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Para(1) = parameter 1
    % Para(2) = parameter 2
    case 'Gamma'
        Para = PriorInfo.Para(ParaNum,:);
        
        if (Value <= 0)
            PP = -1e300;
        else
            % Calculate derivative of log prior
            PP = (Para(1)-1)/Value - 1/Para(2);
        end
        

end

        
        
end

