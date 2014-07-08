function [d,p] = outlier(d,p,values);
% function [d,p] = outlier(d,p,values);
%
% check for spurious signals in the data
% a similar pinging frequency (e.g. pinger or hydrosweep)
% noise at the end of the beams etc.
%
% input  :	d	LADCP analysis data structure
%		n       factor for outlier rejection
%                	   n(1)*rms(data) first sweep
%                	   n(2)*rms(data) second sweep etc.
%		values	LADCP value structure
%
% output :	d	modified LADCP analysis data structure
%
% version 0.2	last change 22.6.2006

% G.Krahmann, LDEO Jun 2000
% history:
% misc editing and cleanup	GK	Jun/2006	0.1->0.2

%
% mark outlier for blocks of 5 minute duration
%
nblock = ceil(5./nmean(diff(d.time_jul)*24*60));
p = setdefv(p,'outlier_n',nblock);

nblock = p.outlier_n;

rw = d.rw(d.izd,:);
rv = d.rv(d.izd,:);
ru = d.ru(d.izd,:);
if isfield(d,'ts')
  rt = d.ts(d.izd,:);
end
dummy = rw*0;

bvel = d.bvel;
dummyb = bvel*0;
if size(dummyb,2)==4, 
  ibvel = 1; 
else, 
  ibvel = 0; 
end

si = size(dummy);
sn = ceil(si(2)/nblock);

lob = length(find(isnan(dummyb)));
lo = length(find(isnan(dummy)));

for n=1:length(p.outlier)
  % calculate anomaly fields
  rwm = nmedian(rw);
  rw = rw - ones(size(rw,1),1)*rwm;
  ru = ru - ones(size(ru,1),1)*nmedian(ru);
  rv = rv - ones(size(rv,1),1)*nmedian(rv);
  if isfield(d,'ts')
    rt = rt - ones(size(rt,1),1)*nmedian(rt);
  end
  if ibvel, 
    bvel(:,3) = bvel(:,3)-rwm'; 
  end
  for m=1:sn
    ind = (m-1)*nblock+[1:nblock];
    ii = find( ind<=si(2) );
    ind = ind(ii);
    dummy2 = dummy(:,ind);
    rrw = rw(:,ind);
    badrw = find(abs(rrw)>p.outlier(n)*rms(rrw(:)));
    rru = ru(:,ind);
    badru = find(abs(rru)>p.outlier(n)*rms(rru(:)));
    rrv = rv(:,ind);
    badrv = find(abs(rrv)>p.outlier(n)*rms(rrv(:)));
    dummy2(badrw) = nan;
    dummy2(badru) = nan;
    dummy2(badrv) = nan;
    if isfield(d,'ts')
      rrt = rt(:,ind);
      badrt = find(abs(rrt)>p.outlier(n)*rms(rrt(:)));
      dummy2(badrt) = nan;
    end
    dummy(:,ind) = dummy2;
    if ibvel
     % bottom track
      bvel(ind,1) = bvel(ind,1)-nmedian(bvel(ind,1));
      bu = find(abs(bvel(ind,1))>p.outlier(n)*rms(bvel(ind,1)));
      dummyb(ind(bu),:) = nan;
      bvel(ind,2) = bvel(ind,2)-nmedian(bvel(ind,2));
      bv = find(abs(bvel(ind,2))>p.outlier(n)*rms(bvel(ind,2)));
      dummyb(ind(bv),:) = nan;
      bw = find(abs(bvel(ind,3))>p.outlier(n)*rms(bvel(ind,3)));
      dummyb(ind(bw),:) = nan;
      hbot(ind) = d.hbot(ind)-nmedian(d.hbot(ind));
      bh = find(abs(hbot(ind))>p.outlier(n)*rms(hbot(ind)));
      dummyb(ind(bh),:) = nan;
    end
  end
  rw = rw+dummy;
  rv = rv+dummy;
  ru = ru+dummy;
  if ibvel, 
    d.bvel = d.bvel + dummyb; 
    d.hbot = d.hbot + dummyb(:,1)'; 
  end
end  

d.weight(d.izd,:) = d.weight(d.izd,:)+dummy;
d.rw(d.izd,:) = d.rw(d.izd,:)+dummy;
d.ru(d.izd,:) = d.ru(d.izd,:)+dummy;
d.rv(d.izd,:) = d.rv(d.izd,:)+dummy;

disp(['    Outlier discarded ',int2str(length(find(isnan(dummy)))-lo),...
	' bins down looking'])
if ibvel
  disp(['    Outlier discarded ',int2str(length(find(isnan(dummyb)))-lob),...
	' bottom track'])
end

if values.up==1

  rw = d.rw(d.izu,:);
  rv = d.rv(d.izu,:);
  ru = d.ru(d.izu,:);
  dummy = rw*0;
  si = size(dummy);
  sn = ceil(si(2)/nblock);

  lo = length(find(isnan(dummy)));

  for n=1:length(p.outlier)
    % calculate anomaly fields
    rw = rw - ones(size(rw,1),1)*nmedian(rw);
    ru = ru - ones(size(ru,1),1)*nmedian(ru);
    rv = rv - ones(size(rv,1),1)*nmedian(rv);
    for m=1:sn
      ind = (m-1)*nblock+[1:nblock];
      ii = find( ind<=si(2) );
      ind = ind(ii);
      dummy2 = dummy(:,ind);
      rrw = rw(:,ind);
      badrw = find(abs(rrw)>p.outlier(n)*rms(rrw(:)));
      rru = ru(:,ind);
      badru = find(abs(rru)>p.outlier(n)*rms(rru(:)));
      rrv = rv(:,ind);
      badrv = find(abs(rrv)>p.outlier(n)*rms(rrv(:)));
      dummy2(badrw) = nan;
      dummy2(badru) = nan;
      dummy2(badrv) = nan;
      dummy(:,ind) = dummy2;
    end
    rw = rw+dummy;
    rv = rv+dummy;
    ru = ru+dummy;
  end  

  d.weight(d.izu,:) = d.weight(d.izu,:)+dummy;
  d.rw(d.izu,:) = d.rw(d.izu,:)+dummy;
  d.ru(d.izu,:) = d.ru(d.izu,:)+dummy;
  d.rv(d.izu,:) = d.rv(d.izu,:)+dummy;

  disp(['    Outlier discarded ',int2str(length(find(isnan(dummy)))-lo),...
	' bins up looking'])

end
