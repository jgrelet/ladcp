function [reslag,i1,i2,co]=bestlag2(a1,a2,nlag)
% function [reslag,i1,i2,co]=bestlag2(a1,a2,nlag)
%
% function to find the best shift for two vectors
% uses median difference
% 
% input :  	a1      : first vector
%          	a2      : second vector
%          	nlag		: maximum lag
%
% output:   reslag  : resulting best lag
%           i1      : usable index vector of series a1
%           i2      : usable index vector of series a2
%           co      : correlation of overlapping parts
%
% version 0.4	last change 02.06.2011


% M. Visbeck LDEO August-2002
% rewrite Gerd Krahmann, IFM-GEOMAR, Sep 2007
% bug fix when handling data set lengths      GK, 03.07.2008  0.2-->0.3
% bug fix when to short data series           GK, 02.06.2011  0.3-->0.4

%
% parse arguments
%
if nargin<3 
  nlag = fix((length(a2)+length(a1))/5); 
end


%
% prepare input data
%
a1 = a1(:);
ind = (~isfinite(a1));
a1(ind) = 0;
a2 = a2(:);
ind = (~isfinite(a2));
a2(ind) = 0;
l1 = length(a1);
l2 = length(a2);
if abs(l2-l1)<50
  a1 = a1(1:min([l1,l2]));
  a2 = a2(1:min([l1,l2]));
elseif l1>l2
  a2 = [a2;zeros(l1-l2,1)];
elseif l2>l1
  a1 = [a1;zeros(l2-l1,1)];
end


%
% check proper setup of calculation
%
if max(abs(nlag)) > length(a1)
  disp(' not enough data to determine lag')
  reslag = 0;
  i1 = [];
  i2 = [];
  co = 1;
  return
end


%
% cross correlate data
%
if exist('xcorr')
  [res,reslag] = xcorr(a1-mean(a1),a2-mean(a2),nlag);
else
  [res,reslag] = xcorr3(a1-mean(a1),a2-mean(a2),nlag);
end
[dummy,resind] = nmax(res);
reslag = reslag(resind);



%
% compute correlation
%
if reslag>0
  aa1 = a1(1+reslag:end);
  aa2 = a2(1:end-reslag);
else
  aa1 = a1(1:end+reslag);
  aa2 = a2(1-reslag:end);
end
good = find( isfinite( aa1 + aa2 ) );
aa1 = aa1(good);
aa2 = aa2(good);
for n=1:4			% remove outliers that degrade correlation
				% this correlation gives warnings when too low
				% is thus not really important for the results
  daa = aa1 - aa2;
  good = find( abs(daa) < 5*nstd(daa) );
  aa1 = aa1(good);
  aa2 = aa2(good);
end
co = corrcoef(aa1,aa2);
co = co(1,2);


%
% prepare joint index vector
%
reslag = -reslag;
i1 = 1:l1;
i2 = i1+reslag;
ii = find(i2>0 & i2<=l2);
i1 = i1(ii); 
i2 = i2(ii); 
