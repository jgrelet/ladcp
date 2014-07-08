function [di,p,data] = calc_ens_av(data,p,values)
% function [di,p,data] = calc_ens_av(data,p,values)
%
% LADCP-2 processing software version 9+
%
% calculate artificial ensemble averages to reduce noise
% and to make matrix sizes more managable 
% 
% depth resolution of the artificial ensembles is typically
% chosen so that no resolution is lost in the final profile
%
% version 0.3	last change xx.07.2008

% Martin Visbeck, LDEO, 6/10/99
% G.Krahmann, IFM-GEOMAR, Jul 2005	partially recoded
% small typos and other stuff               GK, Nov 2007	0.1-->0.2
% changed way shallow profiles are handled  MV, Jul 2008  0.2-->0.3


%
% general function info
%
disp(' ')
disp('CALC_ENS_AV: forming Super-Ensembles')


%
% if p.avdz<0  set default avdz to one bin length
% else assume p.avdz is the depth interval over which is to average
%
if p.avdz<=0
  p.avdz = nmedian(abs(diff(data.izm(:,1)))) * (-p.avdz);
end


%
% another one of those superfluous variables that are never changed  GK
%
avpercent = 100;


%
% hmm don't know, what this is good for  GK
%
if isfield(p,'override_Single_Ping_Err')
  data.down.Single_Ping_Err = p.override_Single_Ping_Err;
end
p = setdefv(p,'superens_std_min',...
	data.down.Single_Ping_Err/sqrt(data.down.Pings_per_Ensemble));


%
% catch special case of shallow profiles
%
if values.maxdepth<100 & isnan(p.avens)
  disp(['    Profile shallower than 100m use fixed ensemble averaging'])
  p.avens = 2;
  p.avdz = 0;
end


if isfinite(p.avens)
  disp(['    Averaging profiles over (p.avens) ',num2str(p.avens),' ensembles'])
else
  disp(['    Averaging profiles over (p.avdz) ',num2str(p.avdz),' m'])
end

if (isnan(p.avdz) | p.avdz<=0 ) & (isnan(p.avens) | p.avens < 2)
  disp('>   avdz=NAN or 0  => No pre-averaging done !!!')
  di.ru = data.ru;
  di.rv = data.rv;
  di.ruvs = data.ru*0+p.superens_std_min;
  di.rw = data.rw;
  di.re = data.re;
  di.ts = data.ts;
  di.tg = data.tg;
  di.weight = data.weight;
  di.bvel = data.bvel';
  di.hbot = data.hbot;
  % make up std
  di.bvels = data.bvel'*0+p.single_ping_accuracy;
  di.hdg = data.hdg;
  di.pit = data.pit;
  di.rol = data.rol;
  di.temp = data.temp;
  di.tsd = data.ts(data.izd(2),:);
  di.tsd_out = data.ts(data.izd(end),:);
  di.dtiv = data.z*0+1;
  di.time_jul = data.time_jul;
  di.z = data.z;
  di.izm = data.izm;
  di.slat = data.slat;
  di.slon = data.slon;
  di.izd = datat.izd;
  di.izu = data.izu;

else
  di.izr = p.all_trusted_i;

  % remove reference velocity and then average ensembles
  disp('    Removing reference velocity and average ensembles ')
  ilast = 1;
  il = length(data.izm);
  im = 0;
  ibin = [1:size(data.ru,1)];
 
  %
  % big loop
  %
  while ilast<il 

    im = im+1;
    i = ilast;

    %
    % set the method of how to select the ensembles which are
    % to be averaged
    % avens>1 averages a fixed number of ensembles
    % else all ensembles with the package within an avdz large
    % interval are averaged
    %
    if p.avens>0
      i1 = ilast+[1:p.avens];
    else
      ii = find(abs(data.izm(1,(ilast+1):il)-data.izm(1,ilast))>p.avdz);
      if length(ii)<1 
        ii = il-ilast; 
      end
      i1 = ilast+[1:ii(1)];
    end

    i1l = length(i1)/2;
    i1 = round(mean(i1)+[-i1l:i1l]);
    ii = find(i1<1); 
    i1(ii) = [];
    ii = find(i1>il); 
    i1(ii) = [];
    if length(i1)==1 
      i1 = [i1 i1]; 
    end
    ilast = max(i1);

    w = data.weight(p.all_trusted_i,i1)*0+1;
    w2 = data.weight(:,i1)*0+1;


    % U
    ur = nmedian(data.ru(p.all_trusted_i,i1).*w);
    ruav = nmean(ur);
    i3 = find(isnan(ur));
    ur(i3) = i3*0;
    iav = round(length(ur)/200 * avpercent);
    ur = meshgrid(ur,ibin);
    di.ru(:,im) = meanmediannan([data.ru(:,i1).*w2-ur]',iav)'+ruav;
    rus = nstd([data.ru(:,i1).*w2]')';
    % V
    vr = nmedian(data.rv(p.all_trusted_i,i1).*w);
    rvav = nmean(vr);
    i3 = find(isnan(vr));
    vr(i3) = i3*0;
    iav = round(length(vr)/200 * avpercent);
    vr = meshgrid(vr,ibin);
    di.rv(:,im) = meanmediannan([data.rv(:,i1).*w2-vr]',iav)'+rvav;
    % estimate mean STD of U and V
    di.ruvs(:,im) = sqrt(rus.^2+nstd([data.rv(:,i1).*w2]')'.^2);
    % W
    wr = nmedian(data.rw(p.all_trusted_i,i1).*w);
    rwav = nmean(wr);
    i3 = find(isnan(wr));
    wr(i3) = i3*0;
    iav = round(length(wr)/200 * avpercent);
    wr = meshgrid(wr,ibin);
    di.rw(:,im) = meanmediannan([data.rw(:,i1).*w2-wr]',iav)'+rwav;

    %EA
    di.re(:,im) = nmean(data.re(:,i1)')';

    %TS
    di.ts(:,im) = nmean(data.ts(:,i1)')';
    di.tg(:,im) = nmean(data.tg(:,i1)')';

    % weight
    di.weight(:,im) = nmean(data.weight(:,i1)')';

    % bottom track
    di.bvel(:,im) = nmean(data.bvel(i1,:))';
    bvel = data.bvel(i1,:);
    %  remove mean vertical velocity from bottom track w prior to STD 
    bvel(:,3) = bvel(:,3)-wr(1,:)';
    di.bvels(:,im) = nstd(bvel)';
    %  distance of bottom
    di.hbot(im) = nmean(data.hbot(i1));

    % bin depth
    di.izm(:,im) = mean(data.izm(:,i1)')';
    
    % heading
    dummy = exp(-sqrt(-1)*(data.hdg(1,i1))*pi/180);
    di.hdg(1,im) = -angle(mean(dummy))*180/pi;
    if values.up==1
      dummy = exp(-sqrt(-1)*(data.hdg(2,i1))*pi/180);
      di.hdg(2,im) = -angle(mean(dummy))*180/pi;
    end

    % pitch and roll
    di.pit(:,im) = mean(data.pit(:,i1),2);
    di.rol(:,im) = mean(data.rol(:,i1),2);

    % target strength
    di.tsd(im) = mean(data.ts(data.izd(2),i1),2);
    di.tsd_out(im) = mean(data.ts(data.izd(end),i1),2);

    % target strength
    di.temp(im) = mean(data.temp(i1));

    % ships position
    di.slon(im) = nmedian(data.slon(i1));
    di.slat(im) = nmedian(data.slat(i1));

    % number of ensembles
    di.dtiv(im) = length(i1);

    % time
    di.time_jul(im) = nmean(data.time_jul(i1));

    % depth
    di.z(im) = nmean(data.z(i1));
    % end of big loop
  end

  % adjust compass to RDI definition [0-360]
  di.hdg = di.hdg+(di.hdg<0)*360;

  % remove outlier
  di.izd = data.izd;
  di.izu = data.izu;
  [di,p] = outlier(di,p,values);

  % remove bottom track data with single ping or large wstd
  if isfinite(p.zbottom)
    % remove single ping bottom track ensembles
    ii = find(prod(di.bvels(1:3,:))==0);
    disp(['    Found ',int2str(length(ii)),...
	' bottom track std==0 set to 0.1 m/s'])
    di.bvels(:,ii) = 0.1;

    ii = find(di.bvels(3,:)>0);
    if length(ii)>0
      disp(['    Found ',int2str(length(ii)),' finite bottom track ensembles'])
      p = setdefv(p,'btrk_wstd',median(di.bvels(3,ii))*2);
      ii = find(di.bvels(3,:)>p.btrk_wstd | di.bvels(3,:)==0);
      disp(['    Discarded ',int2str(length(ii)),...
	' bottom tracks velocities because of wstd  > ',num2str(p.btrk_wstd)])
      di.bvel(:,ii) = NaN;
      di.bvels(:,ii) = NaN;
    else
      ind = find(isfinite(data.bvel(3,:)));
      disp(['    Found no valid bottom track ensemble from ',...
	int2str(length(ind)),' finite raw bottom tracks '])
    end 
  end

  

  %
  % check for small std's
  % blank out data using weight 
  %
  di.ruvs = di.ruvs+di.weight*0;
  ii = find(di.ruvs==0);
  di.weight(ii) = nan;
  disp(['    Setting ',int2str(length(ii)),...
	' weight values to NaN  because super ensemble std=0 '])
  di.ruvs = di.ruvs+di.weight*0;
  ii = find(di.ruvs<p.superens_std_min);
  di.ruvs(ii) = p.superens_std_min;
  disp(['    Setting ',int2str(length(ii)),...
	' values to minimum super ensemble std ',...
	num2str(p.superens_std_min)]);

  disp(['    Reduced ensemble size is ',int2str(length(di.z))])
  if length(di.z)<5
    error('not enough data to process station ')
  end
end

% remove ensembles without velocity data
ii = find(~isfinite(nmax(di.ru)));

if length(ii)>0
    disp(['    Removed ',int2str(length(ii)),...
  	  ' non finite super ensembles'])
    di.rw(:,ii) = [];
    di.ru(:,ii) = [];
    di.rv(:,ii) = [];
    di.ruvs(:,ii) = [];
    di.bvel(:,ii) = [];
    di.bvels(:,ii) = [];
    di.weight(:,ii) = [];
    di.izm(:,ii) = [];
    di.z(ii) = [];
    di.hdg(:,ii) = [];
    di.pit(:,ii) = [];
    di.rol(:,ii) = [];
    di.time_jul(ii)=[];
    di.dtiv(ii) = [];
    di.hbot(ii) = [];
    di.tsd(ii) = [];
    di.temp(ii) = [];
    di.slon(ii) = [];
    di.slat(ii) = [];
end


%
% time difference in seconds (these appear to be central differences GK)
%
dt = diff(di.time_jul)*86400;
di.dt = mean([dt([1,1:end]);dt([1:end,end])]);


%
% SADCP is already averaged
%
di.svel = data.svel;


bin_no = [0];
if length(data.zu) > 0 
  bin_no = [-length(data.zu):1 bin_no]; 
end
if length(data.zd) > 0 
  bin_no = [bin_no 1:length(data.zd)]; 
end

figload('tmp/16.fig',2)
subplot(3,1,2);
imagesc([1:size(di.ts,2)],bin_no,...
	[di.weight(1:length(data.zu),:); ...
	 ones(1,size(di.weight,2))*NaN; ...
	 di.weight(size(di.weight,1)-length(data.zd)+1:end,:)...
        ]);
csc = caxis;
colorbar
xlabel('Super Ensemble #');
ylabel('Bin #');
title('Weights based on various parameters')

streamer([p.name,' Figure 16']);
hgsave('tmp/16')


