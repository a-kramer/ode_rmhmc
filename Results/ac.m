function [ varargout ] = ac( Series , nLags , Q , nSTDs )

[rows , columns]  =  size(Series);

if (rows ~= 1) && (columns ~= 1) 
    error('econ:autocorr:NonVectorInput' , ' Input ''Series'' must be a vector.');
end

rowSeries   =  size(Series,1) == 1;

Series      =  Series(:);       % Ensure a column vector
n           =  length(Series);  % Sample size.
defaultLags =  20;              % BJR recommend about 20 lags for ACFs.

%
% Ensure the number of lags, nLags, is a positive 
% integer scalar and set default if necessary.
%

if (nargin >= 2) && ~isempty(nLags)
   if numel(nLags) > 1
      error('econ:autocorr:NonScalarLags' , ' Number of lags ''nLags'' must be a scalar.');
   end
   if (round(nLags) ~= nLags) || (nLags <= 0)
      error('econ:autocorr:NonPositiveInteger' , ' Number of lags ''nLags'' must be a positive integer.');
   end
   if nLags > (n - 1)
      error('econ:autocorr:LagsTooLarge' , ' Number of lags ''nLags'' must not exceed ''Series'' length - 1.');
   end
else
   nLags  =  min(defaultLags , n - 1);
end

%
% Ensure the hypothesized number of lags, Q, is a non-negative integer
% scalar, and set default if necessary.
%
if (nargin >= 3) && ~isempty(Q)
   if numel(Q) > 1
      error('econ:autocorr:NonScalarQ' , ' Number of lags ''Q'' must be a scalar.');
   end
   if (round(Q) ~= Q) || (Q < 0)
      error('econ:autocorr:NegativeInteger' , ' Number of lags ''Q'' must be a non-negative integer.');
   end
   if Q >= nLags
      error('econ:autocorr:QTooLarge' , ' ''Q'' must be less than ''nLags''.');
   end
else
   Q  =  0;     % Default is 0 (Gaussian white noise hypothisis).
end

%
% Ensure the number of standard deviations, nSTDs, is a positive 
% scalar and set default if necessary.
%

if (nargin >= 4) && ~isempty(nSTDs)
   if numel(nSTDs) > 1
      error('econ:autocorr:NonScalarSTDs' , ' Number of standard deviations ''nSTDs'' must be a scalar.');
   end
   if nSTDs < 0
      error('econ:autocorr:NegativeSTDs' , ' Number of standard deviations ''nSTDs'' must be non-negative.');
   end
else
   nSTDs =  2;     % Default is 2 standard errors (95% condfidence interval).
end

%
% Convolution, polynomial multiplication, and FIR digital filtering are
% all the same operation. The FILTER command could be used to compute 
% the ACF (by computing the correlation by convolving the de-meaned 
% Series with a flipped version of itself), but FFT-based computation 
% is significantly faster for large data sets.
%
% The ACF computation is based on Box, Jenkins, Reinsel, pages 30-34, 188.
%

nFFT =  2^(nextpow2(length(Series)) + 1);
F    =  fft(Series-mean(Series) , nFFT);
F    =  F .* conj(F);
ACF  =  ifft(F);
ACF  =  ACF(1:(nLags + 1));         % Retain non-negative lags.
ACF  =  ACF ./ ACF(1);     % Normalize.
ACF  =  real(ACF);

%
% Compute approximate confidence bounds using the Box-Jenkins-Reinsel 
% approach, equations 2.1.13 and 6.2.2, on pages 33 and 188, respectively.
%

sigmaQ  =  sqrt((1 + 2*(ACF(2:Q+1)'*ACF(2:Q+1)))/n);  
bounds  =  sigmaQ * [nSTDs ; -nSTDs];
Lags    =  [0:nLags]';

if nargout == 0                     % Make plot if requested.

%
%  Plot the sample ACF.
%
   lineHandles  =  stem(Lags , ACF , 'filled' , 'r-o');
   set   (lineHandles(1) , 'MarkerSize' , 4)
   grid  ('on')
   xlabel('Lag')
   ylabel('Sample Autocorrelation')
   title ('Sample Autocorrelation Function (ACF)')
   hold  ('on')
%
%  Plot the confidence bounds under the hypothesis that the underlying 
%  Series is really an MA(Q) process. Bartlett's approximation gives
%  an indication of whether the ACF is effectively zero beyond lag Q. 
%  For this reason, the confidence bounds (horizontal lines) appear 
%  over the ACF ONLY for lags GREATER than Q (i.e., Q+1, Q+2, ... nLags).
%  In other words, the confidence bounds enclose ONLY those lags for 
%  which the null hypothesis is assumed to hold. 
%

   plot([Q+0.5 Q+0.5 ; nLags nLags] , [bounds([1 1]) bounds([2 2])] , '-b');

   plot([0 nLags] , [0 0] , '-k');
   hold('off')
   a  =  axis;
   axis([a(1:3) 1]);

else

%
%  Re-format outputs for compatibility with the SERIES input. When SERIES is
%  input as a row vector, then pass the outputs as a row vectors; when SERIES
%  is a column vector, then pass the outputs as a column vectors.
%
   if rowSeries
      ACF     =  ACF.';
      Lags    =  Lags.';
      bounds  =  bounds.';
   end

   varargout  =  {ACF , Lags , bounds};

end

