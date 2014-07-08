function [d,p,messages]=getbtrack(d,p,values,messages);
% function [d,p,messages]=getbtrack(d,p,values,messages);
%
% create own bottom track in addition to the one used before
%
% version 0.2  last change 28.05.2011

% Changed default bottom track range accoding to bin length    MV Jul 2008  0.7 - 0.8
% allow last but one bin for bottom (before it 
% was erroneously last but two bins)              GK, 28.05.2011  0.1-->0.2

%
% general function info
%
disp(' ')
disp('GETBTRACK:  create own bottom track data in addition to RDI')

disp(['    in: p.btrk_mode ',int2str(p.btrk_mode),' and p.btrk_used ',...
	int2str(p.btrk_used)])

nrdi = sum(isfinite(d.hbot));
disp(['    Found ',int2str(nrdi),' RDI bottom track values '])

% set range of acceptable ranges for bottom track data
p = setdefv(p,'btrk_range',[2 20]*(d.zd(2)-d.zd(1)));

%
%  Check if hbot values are mostly ==0
%
%  some instruments seem to have trouble reporting the distance of the bottom
%  despite giving reasonable bottom track values.
%  If this problem is diagnosed we make the distance ourselves.
%
%
ii1 = sum(isfinite(d.hbot));
ii0 = sum(d.hbot==0);
p.hbot_0 = ii0/(ii1+1)*100;

if p.hbot_0>80
  p.bottomdist = 1;
  disp(['>   Found ',int2str(p.hbot_0),'% of  hbot=0'])
  disp( '>   This means that this instrument does not report good bottom dist.')
end


%
% convert echo amplitude to relative target strength
% at=0.039; % attenuation dB/m for 150 kHz
% at=0.062; % attenuation dB/m for 300 kHz
%
if d.down.Frequency==300,
  p = setdefv(p,'ts_att_dn',0.062);
else
  p = setdefv(p,'ts_att_dn',0.039);
end
d.tg(d.izd,:) = targ(d.ts(d.izd,:),d.zd,p.ts_att_dn);

if values.up==1
  if d.up.Frequency==300,
    p = setdefv(p,'ts_att_up',0.062);
  else
    p = setdefv(p,'ts_att_up',0.039);
  end
  d.tg(d.izu,:) = targ(d.ts(d.izu,:),d.zu,p.ts_att_up);
end


%
% save bin number where down looker starts
%
ib1 = d.izd(1)-1;
nbin = length(d.izd);


if p.btrk_mode>=1
  if p.btrk_ts>0

    disp('    Using high bottom echo amplitudes to create bottom track')
    % don't use first bin
    fitb1 = 1;
    zd = d.zd(fitb1:end);
    ead = d.tg(d.izd(fitb1:end),:);
    % fit parabola to locate bottom
    [zmead,mead,imead] = localmax2(zd',ead);
    imead = imead+fitb1-1;
    dts = mead-ead(1,:);

    % decide which bin to use for bottom velocities
    dz = abs(diff(d.zd(1:2)));
    % imeadbv=round(imead+p.btrk_below);
    imeadbv = round((zmead-d.zd(1))/dz+1+p.btrk_below);

    if p.btrk_used==1
      if p.bottomdist==0
        % check RDI bottom track only if non zero
        
        ii = find(d.hbot<min(p.btrk_range));
        if length(ii)>0
          disp(['    Found ',int2str(length(ii)),...
		' bottom depth closer than btrk_range ',...
            	int2str(min(p.btrk_range))])
            if length(ii)>0.9*nrdi
                warn = ['>    Consider decreasing min of p.btrk_range '];
                disp(warn)
                messages.warn = strvcat(messages.warn,warn);
            end
          d.bvel(ii,:) = nan;
          d.hbot(ii) = nan;
        end
        ii = find(d.hbot>max(p.btrk_range));
        if length(ii)>0
          disp(['    Found ',int2str(length(ii)),...
		' bottom depth beyond btrk_range ',...
            	int2str(max(p.btrk_range))])
            if length(ii)>0.9*nrdi
                warn = ['>    Consider increasing max of p.btrk_range '];
                disp(warn)
                messages.warn = strvcat(messages.warn,warn);
            end
          d.bvel(ii,:) = nan;
          d.hbot(ii) = nan;
        end
      end
      % save RDI bottom track
      d.bvel_rdi = d.bvel;
      d.hbot_rdi = d.hbot;
      ii = find(isfinite(d.bvel(:,1)+d.bvel(:,2)));
      if length(ii)<10, 
        disp('    Found less than 10 RDI bottom track values, try own')
        p.btrk_used = 0; 
      end
    end

    disp(['    Using ',num2str(p.btrk_below),...
         ' bins below maximum target strength for own bottom track velocity'])

    % locate acceptable bottom tracks (don't accept first/last bin)
    ii = find(dts>p.btrk_ts & ...
              zmead>min(p.btrk_range) & zmead<max(p.btrk_range) & ...
               imeadbv<nbin & imeadbv>fitb1);
             
    if length(ii)>0

      % save bottom distance data
      d.hbot_own = d.hbot+NaN;
      d.hbot_own(ii) = zmead(ii);

      % force bottom distance if RDI mode fails to report distance
      if p.bottomdist | p.btrk_mode==2 | p.btrk_used~=1
        d.hbot = d.hbot_own;
        disp(['    Created ',int2str(length(ii)),' bottom distances'])
        if p.bottomdist 
	      p.btrk_used = 12; 
        end
      else
        disp(['    Created ',int2str(length(ii)),...
	' bottom distances keeping original'])
      end

      % make bottom velocity data
      bv = d.bvel+NaN;

      for j=1:length(ii)
        ji = ii(j);
        bv(ji,1) = nmedian(d.raw_u(ib1+imeadbv(ji)+[-1,0,0,1],ji));
        bv(ji,2) = nmedian(d.raw_v(ib1+imeadbv(ji)+[-1,0,0,1],ji));
        bv(ji,3) = nmedian(d.raw_w(ib1+imeadbv(ji)+[-1,0,0,1],ji));
        bv(ji,4) = nmedian(d.raw_e(ib1+imeadbv(ji)+[-1,0,0,1],ji));
      end

      % check for W-bot 
      wref = meanmediannan(d.rw(d.izd,:),2);
      ii = find(abs(wref'-bv(:,3))>p.btrk_wlim);
      disp(['    Removed ',int2str(length(ii)),...
           	' bottom track profiles W_btrk - W_ref difference > ',...
            num2str(p.btrk_wlim)])
      bv(ii,:) = nan;

      % check for outlier
      bv = boutlier(bv,d.hbot_own,p);
      d.bvel_own = bv;
      ii = find(isfinite(bv(:,1)+bv(:,2)));

      if (p.btrk_used~=1 & p.btrk_used~=12) | p.btrk_mode==2
        p.btrk_used = 2;
        d.bvel = d.bvel_own;
        disp(['    Created ',int2str(length(ii)),...
		' bottom track data from normal velocities'])
      else
        disp(['    Created ',int2str(length(ii)),...
		' bottom track velocities keeping original'])
      end

    else
      if (p.btrk_used~=1 & p.btrk_used~=12)
        p.btrk_used = -1;
      end
      disp('    Did not find any bottom echos to create own bottom track ')
    end
  else
    disp('    No valid own bottom track. Increase target strength difference ? ')
  end
else
  disp('    Forcing no use of bottom track data ')
  p.btrk_used = -1;
  d.bvel = d.bvel+nan;
  d.hbot = d.hbot+nan;
end


% summary output
disp(['    out: p.btrk_mode ',int2str(p.btrk_mode),' and p.btrk_used ',int2str(p.btrk_used)])

%=================================
function [bvel,p] = boutlier(bvel,hbot,p);
% function [bvel,p] = boutlier(bvel,hbot,p);
%
%
% input  : bvel bottom track velocity
%           lim       factor for outlier rejection
%                   lim(1)*rms(data) first sweep
%                   lim(2)*rms(data) second sweep etc.
%
% output :	d		changed LADCP analysis data structure
%
% version 0.1	last change 27.6.2000

lim = p.outlier;
nblock = p.outlier_n;

dummyb = bvel*0;
if size(dummyb,2)~=4 
  return 
end

si = size(dummyb);
sn = ceil(si(1)/nblock);

lob=length(find(isnan(dummyb)));

for n=1:length(lim)
  for m=1:sn
    ind = (m-1)*nblock+[1:nblock];
    ii = find( ind<=si(1) );
    ind = ind(ii);
    bvelt(ind,1)=bvel(ind,1)-nmedian(bvel(ind,1));
    bu = find(abs(bvelt(ind,1))>lim(n)*rms(bvelt(ind,1)));
    dummyb(ind(bu),:)=nan;
    bvelt(ind,2)=bvel(ind,2)-nmedian(bvel(ind,2));
    bv = find(abs(bvelt(ind,2))>lim(n)*rms(bvelt(ind,2)));
    dummyb(ind(bv),:)=nan;
    bw = find(abs(bvel(ind,3))>lim(n)*rms(bvel(ind,3)));
    dummyb(ind(bw),:)=nan;
    hbot(ind)=hbot(ind)-nmedian(hbot(ind));
    bh = find(abs(hbot(ind))>lim(n)*rms(hbot(ind)));
    dummyb(ind(bh),:)=nan;
  end
  bvel = bvel + dummyb; 
  hbot = hbot + dummyb(:,1)'; 
end  
disp(['    Boutlier removed ',int2str(length(find(isnan(dummyb)))-lob),...
      ' bottom track velocities '])
return


%================================================
function [ts,bcs]=targ(ea,dis,at,bl,eas,ap)
% function [ts,bcs]=targ(ea,dis,at,bl,eas,ap)
% Target strength of EA for volume scatterer
% ea  = echoamp in  dB
% dis = distance in  m
% at  = attenuation dB/m
% bl  = pulse/bin legth in  m
% eas = source level
% ap  = aperture in degree
% M. Visbeck 2004

% make distance matrix if needed
[lr,lc] = size(ea);

if size(dis,2)==1 | size(dis,1)==1
  dis = dis(:);
  dis = repmat(dis,[1,lc]);
end

if nargin<3
  at = 0.039; % attenuation dB/m for 150 kHz
  at = 0.062; % attenuation dB/m for 300 kHz
end

% binlength
if nargin<4 
  bl = median(abs(diff(dis(:,1)))); 
end

% source level in dB
if nargin <5 
  eas = 100; 
end

% beam aperature in DEGREE convert to radian
if nargin <6 
  ap = 2; 
end
al = ap*pi/180; 

% radius of top and bottom of each bin
r1 = tan(al)*(dis-bl/2);
r2 = tan(al)*(dis+bl/2);

% ensonified volume 
v = pi*bl/3*(r1.^2+r2.^2+r1.*r2);

% transmission loss
tl = 20*log10(dis)+at*dis;

% target strength
ts = ea-eas+2*tl-10*log10(v);

if nargout>1
  %backscatter cross section
  bcs = exp(ts./10);
end
