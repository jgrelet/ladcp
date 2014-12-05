function [c,lags] = xcorr3(x,y,maxlag);
% function [c,lags] = xcorr3(x,y,maxlag);
%
% calculate cross correlation between two vectors
%
% input  : x            - first vector
%          y            - second vector
%          maxlag       - maximum lag in index space
%
% output : c            - correlation vector
%          lags         - lags for the correlation vector
%
% version 0.1  last change 12.07.2012

% Ref: Stearns, SD and David, RA (1988). Signal Processing Algorithms.
%      New Jersey: Prentice-Hall.

% get some sizes
n = max([length(x),length(y)]);
m = 2^nextpow2(n+maxlag);

x = x(:);
y = y(:);

if length(y) > length(x)
  x = [x;zeros(length(y)-length(x),1)];
else
  y = [y;zeros(length(x)-length(y),1)];
end

x = fft(x,2^nextpow2(2*m-1));
y = fft(y,2^nextpow2(2*m-1));

c = ifft(x.*conj(y));

lags = [-maxlag:maxlag];

if maxlag >= m,
  c = [zeros(maxlag-m+1,n^2);c(end-m+2:end,:);c(1:m,:);zeros(maxlag-m+1,n^2)];
else
  c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];
end

