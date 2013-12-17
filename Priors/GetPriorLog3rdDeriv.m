function PP = GetPriorLog3rdDeriv(PriorInfo, Values)

for ParaNum = 1:length(Values)

    Value = Values(ParaNum);

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
                PP(ParaNum) = -1e300;
            else
                PP(ParaNum) = 0;
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
            %PP(ParaNum) = -1/Para(2)^2; % 2nd derivative
            PP(ParaNum) = 0; % 3rd derivative


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gamma prior distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Para(1) = parameter 1
        % Para(2) = parameter 2
        case 'Gamma'
            Para = PriorInfo.Para(ParaNum,:);

            if (Value <= 0)
                PP(ParaNum) = -1e300;
            else
                % Calculate derivative of log prior
                %PP(ParaNum) = -(Para(1)-1)/Value^2; % 2nd derivative
                PP(ParaNum) = 2*(Para(1)-1)/Value^3; % 3rd derivative
            end


    end

end


end

