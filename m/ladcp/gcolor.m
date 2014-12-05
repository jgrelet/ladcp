function [h] = gcolor(x,y,xd,yd,c)
% function [h] = gcolor(x,y,xd,yd,c)
%
% special replacement for pcolor
%
% input  :	x	x-vector (center of faces)
%		y	y-vector (center of faces)
%		xd	vector of width of faces
%		yd	vector of height of faces
%		c	pseudo color matrix
%
% output :	h	handle to object
%
% version 0.1	last change 30.6.2006

% G.Krahmann, IFM-GEOMAR, June 2006

% check input arguments
if nargin~=5
  if nargin==3
    pcolor(x,y,xd)
  else
    error('wrong number of arguments')
  end
end

x = x(:);
y = y(:);
xd = xd(:);
yd = yd(:);

% calculate new x and y-vectors
xn = repmat(min(c(:)),[length(x)*3,1]);
for n=1:length(x)
  xn((n-1)*3+1) = x(n)-xd(n)/2;
  xn((n-1)*3+2) = x(n)+xd(n)/2;
  xn((n-1)*3+3) = x(n)+xd(n)/1.5;
end
yn = repmat(nan,[length(y)*3,1]);
for n=1:length(y)
  yn((n-1)*3+1) = y(n)-yd(n)/2;
  yn((n-1)*3+2) = y(n)+yd(n)/2;
  yn((n-1)*3+3) = y(n)+yd(n)/1.5;
end

% prepare new matrix
cn = repmat(nan,size(c)*3);
cn(1:3:end,1:3:end) = c;
cn(2:3:end,1:3:end) = c;
cn(1:3:end,2:3:end) = c;
cn(2:3:end,2:3:end) = c;

% display matrix
h = pcolor(xn,yn,cn);
