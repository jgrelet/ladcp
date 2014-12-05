function [y,xm,ys] = meanmediannan(x,na)
% function [y,xm,ys] = meanmediannan(x,na)
%
% calculate a mean over -na:na central points
% 'central' here is in a median sense
%
% input  :	x		- data matrix
%		na		- range over which central data should
%				  be averaged
%
% output :	y		- mean over central points
%		xm		- values of central points (same size as x)
%		ys		- standard deviation of central points
%
% version 0.2	last change 13.02.2009

%       M. Visbeck

% handling of empty input              GK, 13.02.2009  0.1-->0.2

if isempty(x)
  y = [];
  xm = [];
  ys = [];
end
if nargin<2 
  na = 0; 
end
na = fix(na);
[m,n] = size(x);

if (m==1)
  x = x.';
end 

[m,n] = size(x);
xm = x+NaN;

for i=1:n
  ii = find(isfinite(x(:,i)));
  if length(ii)>0
    xs = sort(x(ii,i)); 
    indexav = round([-na:na] + length(xs)/2);
    ii = find(indexav>0 & indexav<=length(xs));
    indexav = indexav(ii);
    xm(indexav,i) = xs(indexav);
    y(i) = mean(xs(indexav));
    ys(i) = std(xs(indexav));
  else
   ys(i) = NaN;
   y(i) = NaN;
  end
end
