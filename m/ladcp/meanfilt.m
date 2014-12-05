function [y,stdy] = meanfilt(x,k,dim,method)
% function [y,stdy] = meanfilt(x,k,[dim],[method])
% 
% Running-Mean Filter
% The K values before and after the running index are averaged.
% At the ends the filter-width is decreasing.
%
% input  : 	x			- data vector
%		k			- number of values before and after
%		dim	[1]		- dimension to filter
%		method  ['const']	- method may be 'const' , 'triang' or
%					  'const2'
%
% output :	y		- filtered vector	
%		stdy		- standard deviation of unfiltered data
%				  only available for method 'const'
%
% uses :	nstd.m  nmean.m  meanweight.m
%
% version 0.1.3		last change 17.02.2009

% C.Mertens
%
% changed G.Krahmann, IfM Kiel,   added MODE
% Mar 1994	added output of stdy, G.Krahmann
%
% rewritten for MATLAB 5	G.Krahmann, LODYC Paris Sep 1997, 0.1.0
% removed bug           	G.Krahmann, LODYC Paris Feb 1998, 0.1.0-->0.1.1
% removed bug           	G.Krahmann, LODYC Paris May 1998, 0.1.1-->0.1.2
% catch k==0                    GK, 17.02.2009  0.1.2-->0.1.3

if k==0
  y = x;
  stdy = nan*x;
  return
end

% extract options
if nargin<4
  method='const';
  if nargin==3
    if isstr(dim)
      method = dim;
      dim = 1;
    end
  end
end


% extract dimension
if nargin<3
  dim=1;
else
  if isempty(dim)
    dim=1;
  elseif dim>ndims(x)
    disp('ERROR: wanted dimension does not exist')
    return
  end
end
  

% permute dimensions
sd=size(x);
if sd(1)==1 & ndims(x)==2
  trans=1;
  x=x';
  sd=size(x);
else
  trans=0;
end
if dim>1 | ndims(x)>2
  perm=1;
  dd=[1:length(sd)];
  dd(dim)=nan;
  dd=[dim,dd(find(~isnan(dd)))];
  x=permute(x,dd);
  pd=sd(dd);
  x=reshape(x,[pd(1),prod(pd(2:length(pd)))]);
else
  perm=0;
  pd=sd;
end
%size(x)


% prepare arrays
[m,n]=size(x);
xn=[repmat(nan,[k,n]);x;repmat(nan,[k,n])];
y=zeros(m,n);
if nargout>1
  stdy=zeros(m,n);
end


% check for filter width
if (2*k>m)
  disp('ERROR: The filter-width is too large !!')
  disp('       The data is unchanged !!')
  y=x;
  return
end


% apply filter
if strcmp(method,'const')
  indarr = [0:2*k];
  for i=1:m
    y(i,:) = nmean(xn(i+indarr,:)) ;
  end
  if nargout==2
    for i=1:m
      stdy(i,:) = nstd(xn(i+indarr,:));
    end
  end
elseif strcmp(method,'const2')
  for i=[1:m]+k
    y(i,:) = mean(xn(i+[1:k],:)) ;
  end
  stdy = [];
  disp('ERROR: not tested')
elseif strcmp(method,'triang')
  tri=[[1:k+1],[k:-1:1]]'*ones(1,size(xn,2));
  for i=1:m
    y(i,:) = meanweight(xn(k+[i-k:i+k],:),tri,1,0) ;
  end
  stdy = [];
else
  disp('ERROR: unknown method')
  return
end


% reorder dimensions
y=reshape(y,[size(y,1),pd(2:length(pd))]);
if perm==1
  y=ipermute(y,dd);
end
if trans==1
  y=y';
  if exist('stdy')
    stdy=stdy';
  end
end

