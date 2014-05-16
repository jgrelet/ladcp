function [lag,co]=besttlag(a1,a2,nlag,npoint,nsect)
% function [lag,co]=besttlag(a1,a2,nlag,npoint,nsect)
%
% function to find the best shift for two vectors
% uses median difference
% 
%  input:  a1: first vector (time,variable) should be the high resolution one
%          a2: second vector(time,variable)
%          nlag = number of point to try shifiting 
%          npoint = use how many points
%          nsect = number of chunks of data to use
%
% version 0.2	last change 18.6.2008

% M. Visbeck LDEO August-2002

% fixed warning with NaN in data	GK, June 2008	0.1-->0.2



l1 = size(a1,1);
l2 = size(a2,1);

if nargin<4 
  npoint = min(length(a1),3*nlag); 
end 
if nargin<5 
  nsect = ceil(l1/(npoint)); 
end 
if nargin<3 
  nlag = fix(npoint/8); 
end

% create an index vector dividing the whole data set into nsect chunks
ismed = round(linspace(1,l1,nsect+4));

% avoid beginning and end
ismed([1,2, end-1, end]) = [];

% add largest gradient to list
iok = [fix(l1*0.25) : l1*.75];
[dum,ii] = find(nmax(abs(diff(a1(iok,2)))));
ismed(end+1) = iok(ii);

% loop over chunks of data   
for n=1:length(ismed)

  % select data for section
  isect = [-fix(npoint/2):(npoint/2)]+ismed(n);
  iok = find(isect>0 & isect<l1);
  isect = isect(iok);
 
  % interpolate a2 on a1 with regards to time
  warning off
  good = find(~isnan(a2(:,2)));
  a22=interp1(a2(good,1),a2(good,2),a1(isect,1),'nearest');
  warning on
  a12=a1(isect,2);
  [lagv(n),i1,i2,cov(n)]=bestlag(a12,a22,nlag);
  if lagv(n)==0 & cov(n)==1 
    cov(n)=nan; 
  end

end

if nmax(cov)<0.97
  % try acceleration near large w-difference
  n = length(lagv)+1;
  [lagv(n),i1,i2,cov(n)] = bestlag(diff(a12),diff(a22),nlag);
  disp(['    acceleration lag: ',int2str(lagv(n)),'  correlation: ',...
	num2str(cov(n))])
end

if nmax(cov)<0.97
  % one scan over whole length of max correlation was not great
  n = length(lagv)+1;
  isect = [nlag:(length(a1)-nlag)];
  good = find(~isnan(a2(:,2)));
  a22 = interp1(a2(good,1),a2(good,2),a1(isect,1),'nearest');
  a12 = a1(isect,2);
  [lagv(n),i1,i2,cov(n)] = bestlag(a12,a22,nlag);
  disp(['    all data lag: ',int2str(lagv(n)),'  correlation: ',num2str(cov(n))])
end

if nmax(cov)<0.97
  % one scan over whole length of max correlation was not great
  n = length(lagv)+1;
  ii = find(~isfinite(a12)); 
  a12(ii)=0; 
  ii = find(~isfinite(a22)); 
  a22(ii)=0; 
  [lagv(n),i1,i2,cov(n)] = bestlag(cumsum(a12),cumsum(a22),nlag);
  disp(['    all data integral lag: ',int2str(lagv(n)),'  correlation: ',num2str(cov(n))])
end

% choose the most likely lag
ii = find(~isfinite(cov));
cov(ii) = [];
lagv(ii) = [];
if length(lagv)>2
  lag0 = median(lagv);
else
  lag0 = nan;
end
disp(['    median lag ',int2str(lag0)])
nnlag = [-nlag:nlag];
lagh = hist(lagv,nnlag);
% prefer 0 lag slightly
[nlag1,i1] = max(lagh-abs(nnlag)/(2*max(nnlag)));
lag1 = nnlag(i1);
iok = find(lag1==lagv);
co = nmean(cov(iok));
disp(['    most popular lag ',int2str(lag1)])
[cos,is] = sort(cov);
lags = lagv(is);
lag2 = lags(end);
disp(['    best correlated lag ',int2str(lag2)])
% decide which one to use
if co*sqrt(length(iok)) > cos(end)
  lag = lag1;
else
  lag = lag2;
end
iok = find(lag==lagv);
co = nmax(cov(iok));

disp(['    BESTTLAG:  lag is: ',int2str(lag),'  for which ',...
      int2str(length(iok)/nsect*100),'% of ',int2str(nsect),' lags agree'])

return

%===============================
% the following is a special median definition that does
% not average when there is an even number of elements GK
%================================
function y = median(x,dim)
%MEDIAN Median value.
%   For vectors, MEDIAN(X) is the median value of the elements in X.
%   For matrices, MEDIAN(X) is a row vector containing the median
%   value of each column.  For N-D arrays, MEDIAN(X) is the median
%   value of the elements along the first non-singleton dimension
%   of X.
%
%   MEDIAN(X,DIM) takes the median along the dimension DIM of X.
%
%   Example: If X = [0 1 2
%                    3 4 5]
%
%   then median(X,1) is [1.5 2.5 3.5] and median(X,2) is [1
%                                                         4]
%
%   See also MEAN, STD, MIN, MAX, COV.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.15 $  $Date: 2002/06/05 17:06:39 $

if nargin==1, 
  dim = min(find(size(x)~=1)); 
  if isempty(dim), dim = 1; end
end
if isempty(x), y = []; return, end

siz = [size(x) ones(1,dim-ndims(x))];
n = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),n,prod(siz)/n);

% Sort along first dimension
x = sort(x,1);

if rem(n,2) % Odd number of elements along DIM
  y = x((n+1)/2,:);
else % Even number of elements along DIM
  % y = (x(n/2,:) + x((n/2)+1,:))/2;
  y =  x(fix((n/2)+1),:);
end

% Check for NaNs
y(isnan(x(1,:)) | isnan(x(n,:))) = NaN;

% Permute and reshape back
siz(dim) = 1;
y = ipermute(reshape(y,siz(perm)),perm);
