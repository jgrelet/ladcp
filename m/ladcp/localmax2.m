function [xout,yout,imaxout] = localmax2(X,Y);
% LOCALMAX2 - find local maxima by parabola fitting
% 		[xm,ym]=LOCALMAX2(x,y) returns maxima for every triplet of points (along 1st dim). 
% 		xm(i,j)=NaN means that there is no local maximum of y(:,j)=F{x(:,j)} between x(i,j) and x(i+2,j)
%		x may be the same size as y or as y(:,1) 
%  added to find the three highest points of array

% August 2001, ashcherbina@ucsd.edu
% April 2004, Visbeck visbeck@ldeo.columbia.edu

% fix array sizes
if size(X,1)==size(Y,1) & size(X,2)==1
   X = repmat(X,1,size(Y,2));
end

% find maximum
[dum,imax] = max(Y);
xout = imax+nan;
yout = imax+nan;
imaxout = imax+nan;

% select valid indices for fit
iok = find(imax>1 & imax<(size(Y,1)));

% restrict to valid values
imaxout(iok) = imax(iok);
imaxn = imax(iok);

if length(iok)>0
   disp(['    localmax2: found ',int2str(length(iok)),' valid values'])
else
  disp(['  localmax2 found no valid values'])
  return
end

for i=1:length(imaxn)
  j = iok(i);
  x1(i) = X(imaxn(i)-1,j);
  x2(i) = X(imaxn(i),j);
  x3(i) = X(imaxn(i)+1,j);
  y1(i) = Y(imaxn(i)-1,j);
  y2(i) = Y(imaxn(i),j);
  y3(i) = Y(imaxn(i)+1,j);
end




% compute fit to parabola 

a = (x3.*y2+x1.*y3-x1.*y2-y3.*x2-y1.*x3+y1.*x2)./...
    (x3.*x2.^2-x1.*x2.^2+x1.*x3.^2-x1.^2.*x3+x1.^2.*x2-x3.^2.*x2);
b = -(-x2.^2.*y3+x2.^2.*y1-y2.*x1.^2+y3.*x1.^2-x3.^2.*y1+y2.*x3.^2)./...
    ((-x3+x2).*(x2.*x3-x2.*x1+x1.^2-x3.*x1) );
c = (x2.^2.*y1.*x3-x2.^2.*x1.*y3-x3.^2.*y1.*x2+...
    y3.*x1.^2.*x2+x3.^2.*x1.*y2-y2.*x1.^2.*x3)./...
    ((-x3+x2).*(x2.*x3-x2.*x1+x1.^2-x3.*x1));

% check for maximum in range
xn = x1+nan;
ii = find(a<0);
xn(ii) = -b(ii)./a(ii)/2;

% compute max target strength if desired
xout(iok) = xn;
if nargout>1
  yn = xn.^2.*a+xn.*b+c;
  yout(iok) = yn;
end
