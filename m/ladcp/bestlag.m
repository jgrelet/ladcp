function [lag,i1,i2,co]=bestlag(a1,a2,nlag)
% function [lag,i1,i2,co]=bestlag(a1,a2,nlag)
%
% function to find the best shift for two vectors
% uses median difference
% 
%  input:  a1: first vector
%          a2: second vector
%          nlag = number of point to try shifiting , or vector
%
% M. Visbeck LDEO August-2002

if nargin<3 
  nlag = fix((length(a2)+length(a1))/8); 
end
a1 = a1(:);
a2 = a2(:);

l1 = length(a1);
l2 = length(a2);
lmin = min([l1, l2]);

% set up lag index 
if length(nlag)>1
  ilag = nlag;
else
  ilag = [-nlag:nlag]; 
end

% set up target index 
lag1 = min(ilag);
lag2 = max(ilag);
ic = [max([1,-lag1])+1 : lmin-max([0,lag2])-1];

% reference matrix
c1 = meshgrid(a1(ic),ilag)';

% create lag matrix 
for n=1:length(ilag)
  icc = ic+ilag(n);
  c2(:,n) = a2(icc);
end

% remove NaN and other bad data
igood = find(isfinite(sum(c1')+sum(c2')));

if length(igood)<10

  disp(' not enough data to determine lag')
  lag = 0;
  co = 1;
  
else

  % use median as estimator
  cm = meanmediannan(abs(gradient(c2(igood,:)')'-gradient(c1(igood,:)')'),ceil(length(igood)/6));
  [m1,im1] = min(cm);

  % find best match by fitting parabola to result
  il = floor(length(ilag)/4);
  ifit = [-il:il]+im1;
  ii = find(ifit<1 | ifit>=length(ilag));
  ifit(ii) = [];
  if length(ifit>10);
    c = polyfit(ilag(ifit),cm(ifit),2);
    cmb = cm+mean(cm)*2;
    cmf = cm+nan;
    cmf(ifit) = polyval(c,ilag(ifit));
    cmb(ifit) = cmf(ifit)/6+5/6*cm(ifit);
    [m,im] = min(cmb);
  else
    m = m1;
    im = im1;
  end

  % check correlation of best lag 
  c3 = [c1(igood,1),c2(igood,im)];

  % remove bad spikes that degrade correlation
  nspike = 0;
  for n=1:4
    c3d = diff(c3');
    ii = find(abs(c3d)>(2*(1+n/8)*std(c3d)));
    if (length(c3)-length(ii))<(0.8*length(igood)) 
      break 
    end
    nspike = nspike+length(ii);
    c3(ii,:)=[];
  end
  if nspike>10
%    disp([' bestlag removed ',int2str(nspike),' spikes'])
  end
 
  % compute correlation
  co = corrcoef(c3);
  co = co(1,2);
  lag = ilag(im);
 
end

% prepare joint index vector
i1 = 1:l1;
i2 = i1+lag;
ii = find(i2>0 & i2<=l2);
i1 = i1(ii); 
i2 = i2(ii); 
