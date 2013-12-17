function [ Output ] = LogNormPDF( Values, Means, Variance)

% Values:     D by 1 Vector
% Means:      D by 1 Vector
% Variance:   D by 1 Matrix

if size(Values, 2) > 1
    Values = Values';
end

if size(Means, 2) > 1
    Means = Means';
end

D = length(Values);


Output = sum( -ones(D,1)*(0.5*log(2*pi*Variance)) - ((Values-Means).^2)./(2*(ones(D,1)*Variance)) );

end