function [data,values,messages,params,files] = rdiload(files,params,messages,values)
% function [data,values,messages,params] = rdiload(files,params,messages,values)
%
% RDILOAD Load and merge upward and downward looking ADCP raw data.
%
% data is the main output structure array with the following fields:
%
%    blen: bin length
%    nbin: number of bins
%    blnk: blank after transmit
%    dist: distance of bin 1 from transducer
%     tim: time axis
%     pit: pitch
%     rol: roll
%     hdg: heading
%       s: salinity
%       t: temperature
%      sv: sound velocity
%       u: east velocity
%       v: north velocity
%       w: vertical velocity
%       e: error velocity
%      ts: target strength
%      cm: correlation
%      hb: bottom track distance
%      ub: bottom track east velocity
%      vb: bottom track north velocity
%      wb: bottom track vertical velocity
%      eb: bottom track error velocity
%  instid: instrument serial numbers
%
% version 0.17	last change 14.02.2013

% originally Christian Mertens, IfM Kiel
% G. Krahmann, IFM-GEOMAR, June 2006

% added handling of multiple input files        GK, Dec 2006  	0.1-->0.2
% 3 beam solutions for up/down separately       GK, Jul 2007	0.2-->0.3
% 1200kHz pre-averaging, store raw		GK, Jul/Sep 07	0.3-->0.4
% changed ADCP-ADCP lag handling, more messages	GK, Sep 2007	0.4-->0.5
% warning for 0m blank and no edits, distance
% between up and downlooker                     GK, Mar 2008    0.5-->0.6
% error in calc of single_ping_err              GK, Jun 2008    0.6-->0.7
% added params.up2down to control resampling    GK, 11.07.2008  0.7-->0.8
% removed maxbinrange as parameter, introduced
% clear_ladcp_pressure as parameter             GK, 24.07.2008  0.8-->0.9
% read instrument serial number                 GK, 13.05.2011  0.9-->0.10
% use error velocity in bestlag preparation     GK, 14.05.2011  0.10-->0.11
% rename l to data to be consistent             GK, 28.05.2011  0.11-->0.12
% renamed cosd and sind to cos_d and sin_d      GK, 31.05.2011  0.12-->0.13
% bug in renaming l to data                     GK, 25.08.2011  0.13-->0.14
% time interpolation only up to 10 seconds      GK, 07.03.2012  0.14-->0.15
% recognize multiple files and change 'files'   
% always display % of 3-beam solutions, warn when
% >20%                                          GK, 08.11.2012  0.15-->0.16
% params.timoff_uplooker                        GK, 14.02.2013  0.16-->0.17

%
% general function info
%
disp(' ')
disp('RDILOAD:  load ADCP data')


%
% check whether there are multiple files
% 
count = 1;
while exist(files.ladcpdo(end,:),'file')
  nname = files.ladcpdo(end,:);
  nname(end-[6:-1:4]) = int2str0(count,3);
  if exist(nname,'file')
    files.ladcpdo(end+1,:) = nname;
    count = count+1;
  else
    break
  end
end
count = 1;
while exist(files.ladcpup(end,:),'file')
  nname = files.ladcpup(end,:);
  nname(end-[6:-1:4]) = int2str0(count,3);
  if exist(nname,'file')
    files.ladcpup(end+1,:) = nname;
    count = count+1;
  else
    break
  end
end


%
% parse input structures
%
fdown = files.ladcpdo;
fup = files.ladcpup;
pglim = params.pglim;
elim = params.elim;
tssave = params.ts_save;
cmsave = params.cm_save;
pgsave = params.pg_save;

% fixed leader
f.nbin = 1;   % number of depth cells
f.npng = 2;   % pings per ensemble
f.blen = 3;   % depth cell length
f.blnk = 4;   % blank after transmit
f.dist = 5;   % distance to the middle of the first depth cell
f.plen = 6;   % transmit pulse length
% why was the following implemented as a STRING ????   GK
f.serial = 7:14; % serial number of CPU board
f.bandwidth = 15; % bandwidth
f.serial_instrument = 16:19; % serial number of the instrument

% variable leader
v.tim = 1;    % true time (Julian days)
v.pit = 2;    % pitch
v.rol = 3;    % roll
v.hdg = 4;    % heading
v.t   = 5;    % temperature
v.s   = 6;    % salinity
v.sv  = 7;    % sound velocity
v.xmc = 8;    % transmit current
v.xmv = 9;    % transmit volt
v.tint = 10;  % internal temperature
v.pres = 11;  % internal pressure
v.presstd = 12; % internal pressure standard deviation


%
% time_matching preset, will be set to one, if multiple files are
% encountered
%
%time_matching = 0;


%
% load downward looking ADCP
%
[fid_dn,message] = fopen(fdown(1,:),'r','l');
if fid_dn == -1
  message = sprintf('%s: %s',fdown(1,:),message);
  disp('>   Was not able to open down looking RDI file')
  disp(message)
  error('>     terminating LADCP processing')
end
disp(['    Loading down-data ',fdown(1,:)])


%
% call the loading functions
%
% depending on the type of instrument this will be
% nbread.m    or   whread.m
%
if isbb(fid_dn)

  disp('    Detected BB or WH data')
  if size(fdown,1)==1
    [fd,vd,veld,cmd,ead,pgd,btd] = whread(fid_dn);
  else
    fclose(fid_dn);
    [fd,vd,veld,cmd,ead,pgd,btd] = whread_multi(fdown);
    fid_dn = fopen(fdown(1,:),'r','l');
   % time_matching = 1;
  end    
  data.bbadcp = 1;


  % check for short blank and no masked bins
  if fd(f.blnk)==0 && isempty(params.edit_mask_dn_bins)
    disp('>   Found 0m blank length and no masking of first downlooker bin.')
    disp('>   Recommend setting  p.edit_mask_dn_bins  .')
  end


else

  %
  % close and open file since the file format changed
  %
  fclose(fid_dn);
  fid_dn = fopen(fdown,'r','b');
  disp('    Detected NB data')
  [fd,vd,veld,swd,ead,pgd] = nbread(fid_dn);
  data.bbadcp = 0;


  %
  % add some dummy values for data that does not exist in NB instruments
  %
  cmd = ead*0+100;
  btd = NaN*ones(size(vd,1),1,16);


  %
  % the following appears to be a detection of ALL ZEROS records
  % they are being replaced by NaNs
  % 
  ok = double( prod(veld,3)~=sum(veld,3) );
  [ok,ind] = replace( ok, ok==0, nan); 
  if ~isempty(ind)
    disp('>   Special problem with NB data:')
    disp(['      replaced ',int2str(length(ind)),' records because of all 0s'])
    veld = veld.*repmat(ok,[1,1,4]);
  end

end
l1 = size(veld,1);
l2 = size(veld,2);
disp(['    Read ',int2str(l1),' ensembles with ',int2str(l2),' bins each']) 


%
% do the full bin removal, if requested
%
if ~isempty(params.edit_hardremove_mask_dn_bins);
  veld(:,params.edit_hardremove_mask_dn_bins,:) = nan;
  pgd(:,params.edit_hardremove_mask_dn_bins,:) = nan;
  ead(:,params.edit_hardremove_mask_dn_bins,:) = nan;
  cmd(:,params.edit_hardremove_mask_dn_bins,:) = nan;
end


%
% check for data stored in beam coordinates and rotate them to
% earth coordinates
%
dd = rditype(fid_dn);
data.down = dd;
data.down.rdi_std = calc_rdi_std(data.down);
values.rdi_std(1) = data.down.rdi_std;
if dd.Coordinates==0
  disp('>   Detected BEAM coordinates: rotating to EARTH coordinates')
  if dd.Frequency==1200
    data.d1200.veld_beam = veld;
  end
  veld = b2earth(veld,vd,dd);
end

% check for 1200kHz WH
if dd.Frequency==1200
  disp(['>   DETECTED 1200kHz WH: averaging bins, NAV=',int2str(params.nav_1200)])
  zd = fd(f.dist) + fd(f.blen)*([1:fd(f.nbin)] - 1);
  nbins = fd(1);
  count = 1;

  % backup unaveraged data
  data.d1200.veld = veld;
  data.d1200.cmd = cmd;
  data.d1200.ead = ead;
  data.d1200.pgd = pgd;
  data.d1200.fd = fd;

  % blank out first few meters when being raised
  nearbins = find(zd<params.extra_blank_1200);
  vertvel = nmean( squeeze( veld(:,:,3) )' );
  upvel = (vertvel<0);  % removed a find GK 20.05.2011
  veld(upvel,nearbins,:) = nan;
   
  % a 1200kHz measures very close to the rosette
  % high error velocities result when measured in the eddy tail of the rosette
  % blank out such data
  bad = find(isnan(veld));
  cmd(bad) = nan;
  ead(bad) = nan;
  pgd(bad) = nan;
  for n=1:params.nav_1200:nbins 
    ind = [1:params.nav_1200]+n-1;
    ind = ind( ( ind<=nbins ) );   % removed a find  GK, 20.05.2011
    nveld(:,count,:) = nmean(veld(:,ind,:),2);
    ncmd(:,count,:) = nmean(cmd(:,ind,:),2);
    nead(:,count,:) = nmean(ead(:,ind,:),2);
    npgd(:,count,:) = nmean(pgd(:,ind,:),2);
    count = count+1;
  end
  veld = nveld;
  cmd = ncmd;
  ead = nead;
  pgd = npgd;
  fd(f.nbin) = ceil( fd(f.nbin)/params.nav_1200 );
  fd(f.blen) = fd(f.blen)*params.nav_1200;
  fd(f.dist) = (fd(f.dist)-fd(f.blen)/params.nav_1200/2)+fd(f.blen)/2;
end    

%
% remove extra (?) bottom track dimension
%
btd = squeeze(btd);
if ndims(btd)>2
  disp('>   Removal of extra bottom track dimension failed !!!')
end


%
% median echoamplitude and correlation
%
ead_m = squeeze(nmedian(ead));
cmd_m = squeeze(nmedian(cmd));
if tssave(1)~=0
  ead_all = ead(:,:,tssave);
end
if cmsave(1)~=0
  cmd_all = cmd(:,:,cmsave);
end
if pgsave(1)~=0
  pgd_all = pgd(:,:,pgsave);
end
ead = median(ead,3);
cmd = median(cmd,3);

if sum(isfinite(btd(:)))>0
  data.btrk_used = 1;
else
  data.btrk_used = 0;
end


%
% transform to earth coordinates
%
if data.btrk_used == 1
  db = rditype(fid_dn);
  if db.Coordinates==0
    for n=1:4
      velb(:,1,n) = btd(:,4+n);
      velb(:,2,n) = NaN;
    end
    disp('>   DETECTED BEAM bottom track coordinates, rotating to EARTH coord.')
    db.use_binremap = 0;
    velb = b2earth(velb,vd,db);
    for n=1:4
      btd(:,4+n) = velb(:,1,n);
    end
  end
end


%
% finally we can close the file
%
fclose(fid_dn);


%
% apply percent-good threshold
%
pgd = pgd(:,:,4);
ind = (pgd < pglim);
pgd(ind) = NaN;
pgd(~ind) = 1;
if any(ind)
  disp(sprintf('    Removed %d downlooker values because of percent good < %g',...
		length(find(ind)),pglim));
  veld = veld.*repmat(pgd,[1,1,4]);
end


% 
% load upward looking ADCP
%
up = ~isempty(fup);

if up
  fid_up = fopen(fup(1,:),'r','l');
  if fid_up == -1
    disp(['>   Got file name for up-looker but could not open file'])
    disp(['>     Filename was: ',fup(1,:)])
    disp(['>     Continuing with just down-looker.'])
    up = 0;
    values.up = 0;
  else
    values.up = 1;
  end
else
  values.up = 0;
end

if values.up==1
  disp(['    Loading up-data ',fup(1,:)])
  if size(fup,1)==1
    [fu,vu,velu,cmu,eau,pgu,btu] = whread(fid_up);
  else
    fclose(fid_up);
    [fu,vu,velu,cmu,eau,pgu,btu] = whread_multi(fup);
    [fid_up,message] = fopen(fup(1,:),'r','l');
%    time_matching = 1;
  end    


  % check for short blank and no masked bins
  if fu(f.blnk)==0 && isempty(params.edit_mask_up_bins)
    disp('>   Found 0m blank length and no masking of first uplooker bin.')
    disp('>   Recommend setting  p.edit_mask_up_bins  .')
  end


  % add distance between uplooker and downlooker to the first
  % bin distance of the uplooker
  % this will do a HARD change of the data loaded from the instrument !!!
  fu(f.dist) = fu(f.dist) + params.dist_up_down;
    

  % check for beam coordinates
  du = rditype(fid_up);
  fclose(fid_up);
  data.up = du;
  data.up.rdi_std = calc_rdi_std(data.up);
  values.rdi_std(2) = data.up.rdi_std;
  if du.Coordinates==0
    disp('>   DETECTED BEAM coordinates: rotating to EARTH coordinates')
    if du.Frequency==1200
      data.d1200.velu_beam = velu;
    end
    velu = b2earth(velu,vu,du);
  end
  %
  
  % check for 1200kHz WH
  if du.Frequency==1200
    disp(['>   DETECTED 1200kHz WH: averaging bins, NAV=',...
	int2str(params.nav_1200)])
    zu = fu(f.dist) + fu(f.blen)*([1:fu(f.nbin)] - 1);
    nbins = fu(1);
    count = 1;

    % backup unaveraged data
    data.d1200.velu = velu;
    data.d1200.cmu = cmu;
    data.d1200.eau = eau;
    data.d1200.pgu = pgu;
    data.d1200.fu = fu;
    data.d1200.vu = vu;

    % blank out first few meters when being lowered
    nearbins = find(zu<params.extra_blank_1200);
    vertvel = nmean( squeeze( velu(:,:,3) )' );
    downvel = find(vertvel>0);
    velu(downvel,nearbins,:) = nan;
    
    bad = find(isnan(velu));
    cmu(bad) = nan;
    eau(bad) = nan;
    pgu(bad) = nan;
    if params.nav_1200<=nbins
      for n=1:params.nav_1200:nbins 
        ind = [1:params.nav_1200]+n-1;
        ind = ind( find( ind<=nbins ) );
        nvelu(:,count,:) = nmean(velu(:,ind,:),2);
        ncmu(:,count,:) = nmean(cmu(:,ind,:),2);
        neau(:,count,:) = nmean(eau(:,ind,:),2);
        npgu(:,count,:) = nmean(pgu(:,ind,:),2);
        count = count+1;
      end
      fu(f.nbin) = ceil( fu(f.nbin)/params.nav_1200 );
      fu(f.blen) = fu(f.blen)*params.nav_1200;
      fu(f.dist) = (fu(f.dist)-fu(f.blen)/params.nav_1200/2)+fu(f.blen)/2;
    else
      nvelu(:,1,:) = nmean(velu,2);
      ncmu(:,1,:) = nmean(cmu,2);
      neau(:,1,:) = nmean(eau,2);
      npgu(:,1,:) = nmean(pgu,2);
      nvelu(:,2,:) = nmean(velu,2);
      ncmu(:,2,:) = nmean(cmu,2);
      neau(:,2,:) = nmean(eau,2);
      npgu(:,2,:) = nmean(pgu,2);
      nvelu(:,3,:) = nmean(velu,2);
      ncmu(:,3,:) = nmean(cmu,2);
      neau(:,3,:) = nmean(eau,2);
      npgu(:,3,:) = nmean(pgu,2);
      fu(f.dist) = (fu(f.blen)*fu(f.nbin))/2;
      fu(f.nbin) = 3;
    end 
    velu = nvelu;
    cmu = ncmu;
    eau = neau;
    pgu = npgu;
  end    

  l1 = size(velu,1);
  l2 = size(velu,2);
  btu = squeeze(btu);
  disp(['    Read ',int2str(l1),' ensembles with ',int2str(l2),' bins each']) 
  % median echoamplitude and correlation
  %eau = targs(mean(eau,3)',z(:))';
  eau_m = squeeze(nmedian(eau));
  cmu_m = squeeze(nmedian(cmu));
  if tssave(1)~=0
    eau_all = eau(:,:,tssave);
  end
  if cmsave(1)~=0
    cmu_all = cmu(:,:,cmsave);
  end
  if pgsave(1)~=0
    pgu_all = pgu(:,:,pgsave);
  end
  eau = median(eau,3);
  cmu = median(cmu,3);

  % apply percent-good threshold
  pgu = pgu(:,:,4);
  ind = (pgu < pglim);
  if any(ind)
    disp(sprintf('    Removed %d uplooker values because of percent good < %g',...
	 	  length(find(ind)),pglim));
  end
  pgu(ind) = NaN;
  pgu(~ind) = 1;
  for k = 1:4
    velu(:,:,k) = velu(:,:,k).*pgu;
  end

end

% distance vectors
z = fd(f.dist) + fd(f.blen)*([1:fd(f.nbin)] - 1);
idb = [1:fd(f.nbin)];

if values.up==1
  iub = [1:fu(f.nbin)];

  %
  % In case the time of the instruments differs
  % we are going to shift the up-looker.
  % This part is rather nasty as a proper check would 
  % require quite a bit of extra code.
  %
  if params.timoff_uplooker~=0
    vu(:,1,v.tim) = vu(:,1,v.tim)+params.timoff_uplooker;
    figure(3)
    clf
    plot(vd(:,1,v.tim),veld(:,2,3))
    hold on
    plot(vu(:,1,v.tim),velu(:,2,3),'r')
    disp(' ')
    disp('PAUSED')
    disp('check whether params.timoff_uplooker brings vertical velocities to alignment')
    disp(' ')
    pause
    close(3)
  end

  % 
  % in case instruments are not synchronized and one instruments starts late, then
  % the following routines will fail
  % to solve the problem we need to prepad the data set of the late starting 
  % instrument
  %
  if vu(1,1,v.tim)-vd(1,1,v.tim)>10/86400
    disp('data at beginning of uplooker is missing')
    dt = nmean(diff(vu(:,1,v.tim)));
    extra_ind = round((vu(1,1,v.tim)-vd(1,1,v.tim))/dt);
    dummy = nan*vu(ones(extra_ind,1),1,:);
    vu = [dummy;vu];
    vu(1:extra_ind,1,v.tim) = vu(extra_ind+1,1,v.tim)-[extra_ind:-1:1]*dt;
    vu(1:extra_ind,1,v.pit) = 0;
    vu(1:extra_ind,1,v.rol) = 0;
    vu(1:extra_ind,1,v.hdg) = 0;
    velu = velu([ones(extra_ind,1);[1:size(velu,1)]'],:,:);
    velu(1:extra_ind,:,:) = nan;
    pgu = pgu([ones(extra_ind,1);[1:size(pgu,1)]'],:,:);
    pgu(1:extra_ind,:,:) = nan;
    eau = eau([ones(extra_ind,1);[1:size(eau,1)]'],:,:);
    eau(1:extra_ind,:,:) = nan;
    cmu = cmu([ones(extra_ind,1);[1:size(cmu,1)]'],:,:);
    cmu(1:extra_ind,:,:) = nan;
    btu = btu([ones(extra_ind,1);[1:size(btu,1)]'],:,:);
    btu(1:extra_ind,:,:) = nan;
  elseif vd(1,1,v.tim)-vu(1,1,v.tim)>20/86400
    disp('data at beginning of downlooker is missing')
    dt = nmean(diff(vd(:,1,v.tim)));
    extra_ind = round((vd(1,1,v.tim)-vu(1,1,v.tim))/dt);
    dummy = nan*vd(ones(extra_ind,1),1,:);
    vd = [dummy;vd];
    vd(1:extra_ind,1,v.tim) = vd(extra_ind+1,1,v.tim)-[extra_ind:-1:1]*dt;
    vd(1:extra_ind,1,v.pit) = 0;
    vd(1:extra_ind,1,v.rol) = 0;
    vd(1:extra_ind,1,v.hdg) = 0;
    veld = veld([ones(extra_ind,1);[1:size(veld,1)]'],:,:);
    veld(1:extra_ind,:,:) = nan;
    pgd = pgd([ones(extra_ind,1);[1:size(pgd,1)]'],:,:);
    pgd(1:extra_ind,:,:) = nan;
    ead = ead([ones(extra_ind,1);[1:size(ead,1)]'],:,:);
    ead(1:extra_ind,:,:) = nan;
    cmd = cmd([ones(extra_ind,1);[1:size(cmd,1)]'],:,:);
    cmd(1:extra_ind,:,:) = nan;
    btd = btd([ones(extra_ind,1);[1:size(btd,1)]'],:,:);
    btd(1:extra_ind,:,:) = nan;
%  elseif vd(1,1,v.tim)-vu(1,1,v.tim)>20/86400
%    disp('data at beginning of downlooker is missing')
%    error('this is not yet coded. Contact Gerd Krahmann for help.')
  end


  %
  % check if ping rate is the same for both instruments
  %
  % removed 'timediff' because I could not see the real difference
  % to 'pingdiff'  GK
  %
  timd = vd(:,1,v.tim)';
  timu = vu(:,1,v.tim)';

  if abs(nmax(diff(timd))+nmax(-diff(timd))) > 0.05/24/3600
    warn=('>   Ping rate varies in down instrument');
    disp(warn)
    disp(['>   Min down ping rate :',num2str(-86400*nmax(-diff(timd))),...
       '  max down ping rate :',num2str(86400*nmax(diff(timd)))])

  end

  if  abs(nmedian(diff(timd))-nmedian(diff(timu))) > 0.05/24/3600  |  params.force_resample_uplooker==1
    warn=('>   Average ping rates differ between instruments ');
    disp(warn)
    messages.warn = strvcat(messages.warn,warn);
    disp(['>   Avg down ping rate :',num2str(86400*nmean(diff(timd))),...
         '  avg up ping rate :',num2str(86400*nmean(diff(timu)))])

    if params.up2down==1
      disp(['>   Resampling up instrument to down instrument''s timing.'])
      iu = [1:length(timd)];
      ii = (iu>length(timu));  % removed a find  GK, 20.05.2011
      iu(ii) = length(timu);
      for i=find(isfinite(timd))
        [m,iu(i)] = min(abs(timu-timd(i)));
      end
      ilast = min(length(iu),length(timu));
      disp(['>   Cast ends now differ by ',num2str(iu(ilast)-ilast),...
          ' ensembles'])
      id = [1:length(timd)];
    else
      disp(['>   Resampling down instrument to up instrument''s timing.'])
      id = [1:length(timu)];
      ii = (id>length(timd));  % removed a find  GK, 20.05.2011
      id(ii) = length(timd);
      for i=find(isfinite(timu))
        [m,id(i)] = min(abs(timd-timu(i)));
      end
      ilast = min(length(id),length(timd));
      disp(['>   Cast ends now differ by ',num2str(id(ilast)-ilast),...
          ' ensembles'])
      iu = [1:length(timu)];
      if any(isnan(id))
        1
      end
    end
  else
    id = [1:length(timd)];
    iu = [1:length(timu)];
  end


  % 
  % find best lag between two instruments by matching vertical velocity 
  %
  wu = squeeze(velu(iu,:,3));
  dummy = squeeze(velu(iu,:,4));
  bad = find(isnan(dummy));
  if ~isempty(bad)
    wu(bad) = nan;
  end
  wd = squeeze(veld(id,:,3));
  dummy = squeeze(veld(id,:,4));
  bad = find(isnan(dummy));
  if ~isempty(bad)
    wd(bad) = nan;
  end
  wb2u = nmedian(wu');
  wb2d = nmedian(wd');

  if params.up2down==1
    [lag,iiu,id,co] = bestlag2(wb2u,wb2d);
  else
    [lag,iu,iid,co] = bestlag2(wb2u,wb2d);
  end
  if params.bestlag_testing_on==1 | co<0.9
    [lag1,iiu1,id1,co1] = bestlag(wb2u,wb2d,params.maxlag);
    disp('BESTLAG confirmation using old routine:')
    disp(['old result : ',int2str(lag1),' time steps   correlation : ',...
    	num2str(co1)])
    disp(['new result : ',int2str(lag1),' time steps   correlation : ',...
      num2str(co)])
    if lag1~=lag
      figure(3)
      clf
      plot(wb2u)
      hold on
      plot(wb2d,'r')
      warn = 'DIFFERENT LAG RESULTS !!!  CHECK NEW ROUTINE';
      disp(warn)
      title(warn)
      xlabel('time')
      ylabel('vertical velocity')
      messages.warn = strvcat(messages.warn,warn);
    end
  end
  if lag>20
    warn = ['>   Found LARGE timing difference between ADCPs !'];
    messages.warn = strvcat(messages.warn,warn);
  end 

  disp(['    Shifting ADCP timeseries by ',num2str(lag),' ensembles']);
  disp(['    Time-lag: ',num2str(lag),' time steps    Correlation : ',num2str(co)])
  if params.up2down==1
    iu = iu(iiu);
  else
    id = id(iid);
  end

  disp(['    Number of joint ensembles is : ',num2str(length(iu))]);


  %
  % parse instrument serial numbers
  %
  data.serial_cpu_u = fu(f.serial);
  data.serial_cpu_d = fd(f.serial);
  values.instid(1) = prod(data.serial_cpu_d+1) + sum(data.serial_cpu_d);
  values.instid(2) = prod(data.serial_cpu_u+1) + sum(data.serial_cpu_u);
  data.serial_inst_u = sum( fu(f.serial_instrument).*[1,256,256^2,256^3]);
  data.serial_inst_d = sum( fd(f.serial_instrument).*[1,256,256^2,256^3]);
  values.inst_serial(1) = data.serial_inst_d;
  values.inst_serial(2) = data.serial_inst_u;
  if abs(values.inst_serial(1))>1e6 | abs(values.inst_serial(2))>1e6
    disp('>   Can not interpret serial number in raw data file.')
    disp('>   This might happen because of an old firmware version')
  else
    params.down_sn = values.inst_serial(1);
    params.up_sn = values.inst_serial(2);
  end
  if ~isnan(values.inst_serial(1))
    disp(['    Down looking instrument ID : ',int2str(values.inst_serial(1))])
  end
  if ~isnan(values.inst_serial(2))
    disp(['    Up looking instrument ID   : ',int2str(values.inst_serial(2))])
  end

  %
  % merge upward and downward
  %
  data.zu=[0:(fu(f.nbin)-1)]*fu(f.blen)+fu(f.dist);
  data.zd=[0:(fd(f.nbin)-1)]*fd(f.blen)+fd(f.dist);
  data.serial_cpu_u = fu(f.serial);
  data.serial_cpu_d = fd(f.serial);
  data.npng_u = fu(f.npng);
  data.npng_d = fd(f.npng);
  data.nens_u = size(vu,1);
  data.nens_d = size(vd,1);
  data.blen = fd(f.blen);
  data.nbin = fd(f.nbin);
  data.blnk = fd(f.blnk);
  data.dist = fd(f.dist);
  data.bandwidth = [fd(f.bandwidth),fu(f.bandwidth)];
  data.tim = [vd(id,1,v.tim),vu(iu,1,v.tim)]';
  data.pit = [vd(id,1,v.pit),vu(iu,1,v.pit)]';
  data.rol = [vd(id,1,v.rol),vu(iu,1,v.rol)]';
  data.hdg = [vd(id,1,v.hdg),vu(iu,1,v.hdg)]';
  data.s = [vd(id,1,v.s),vu(iu,1,v.s)]';
  data.t = [vd(id,1,v.t),vu(iu,1,v.t)]';
  data.sv = [vd(id,1,v.sv),vu(iu,1,v.sv)]';
  data.xmc = [vd(id,1,v.xmc),vu(iu,1,v.xmc)]';
  data.xmv = [vd(id,1,v.xmv),vu(iu,1,v.xmv)]';
  data.tint = [vd(id,1,v.tint),vu(iu,1,v.tint)]';
  data.pres = [vd(id,1,v.pres),vu(iu,1,v.pres)]';
  data.presstd = [vd(id,1,v.presstd),vu(iu,1,v.presstd)]';
  data.u = [fliplr(velu(iu,iub,1)) veld(id,idb,1)]';
  data.v = [fliplr(velu(iu,iub,2)) veld(id,idb,2)]';
  data.w = [fliplr(velu(iu,iub,3)) veld(id,idb,3)]';
  data.e = [fliplr(velu(iu,iub,4)) veld(id,idb,4)]';
  data.ts = [fliplr(eau(iu,iub)) ead(id,idb)]';
  data.cm = [fliplr(cmu(iu,iub)) cmd(id,idb)]';
% No reason to keep this since pgu and pgd don't mean much anymore  
%l.pg = [fliplr(pgu(iu,iub)) pgd(id,idb)]';
  if tssave(1)~=0
    data.ts_all_u = eau_all(iu,iub,:);
    data.ts_all_d = ead_all(id,idb,:);
  end
  if cmsave(1)~=0
    data.cm_all_u = cmu_all(iu,iub,:);
    data.cm_all_d = cmd_all(id,idb,:);
  end
  if pgsave(1)~=0
    data.pg_all_u = pgu_all(iu,iub,:);
    data.pg_all_d = pgd_all(id,idb,:);
  end
% distance to surface
  hs = median(btu(iu,1:4),2)';
  if sum(isfinite(hs))>1 & nsum(hs)>0
    data.hs = hs;
  else
    % try to use targetstrength to find surface
    if sum(isfinite(eau))>1
      disp('    Using target strength of up looking ADCP to find surface ')
      eaum = nmedian(eau);
      eaua = eau-meshgrid(eaum,eau(:,1));
      [eam,hsb] = max(eaua(iu,:)');
      data.hs = data.zu(hsb);
      ii = (eam<20);  % removed a find  GK, 20.05.2011
      l.hs(ii) = NaN;
      ii = find(hsb==1 | hsb==size(eau,2));
      data.hs(ii) = NaN;
    end
  end

  data.hb = median(btd(id,1:4),2)';
  data.hb4 = btd(id,1:4)';
  data.ub = btd(id,5)';
  data.vb = btd(id,6)';
  data.wb = btd(id,7)';
  data.eb = btd(id,8)';
  data.tsd_m = ead_m;
  data.cmd_m = cmd_m;
  data.tsu_m = eau_m;
  data.cmu_m = cmu_m;
  for n = 1:4
    [dummy,ir] = min(abs(data.cmu_m(:,n)-max(data.cmu_m(1,:))*0.3));
    params.up_range(n) = data.zu(ir);
    [dummy,ir] = min(abs(data.cmd_m(:,n)-max(data.cmd_m(1,:))*0.3));
    params.dn_range(n) = data.zd(ir);
  end
  
else

  data.zd=[0:(fd(f.nbin)-1)]*fd(f.blen)+fd(f.dist);
  data.serial_cpu_d = fd(f.serial);
  data.blen = fd(f.blen);
  data.npng_d = fd(f.npng);
  data.nens_d = length(vd(v.tim));
  data.nbin = fd(f.nbin);
  data.blnk = fd(f.blnk);
  data.dist = fd(f.dist);
  data.bandwidth = fd(f.bandwidth);
  data.tim = vd(:,1,v.tim)';
  data.pit = vd(:,1,v.pit)';
  data.rol = vd(:,1,v.rol)';
  data.hdg = vd(:,1,v.hdg)';
  data.s = vd(:,1,v.s)';
  data.t = vd(:,1,v.t)';
  data.sv = vd(:,1,v.sv)';
  data.xmc = vd(:,1,v.xmc)';
  data.xmv = vd(:,1,v.xmv)';
  data.tint = vd(:,1,v.tint)';
  data.pres = vd(:,1,v.pres)';
  data.presstd = vd(:,1,v.presstd)';
  data.hdg = vd(:,1,v.hdg)';
  data.u = veld(:,idb,1)';
  data.v = veld(:,idb,2)';
  data.w = veld(:,idb,3)';
  data.e = veld(:,idb,4)';
  data.ts = ead(:,idb)';
  if tssave(1)~=0
    data.ts_all_d = ead_all(:,idb,:);
  end
  data.cm = cmd(:,idb)';
  if cmsave(1)~=0
    data.cm_all_d = cmd_all(:,idb,:);
  end
  data.pg = pgd(:,idb)';
  if pgsave(1)~=0
    data.pg_all_d = pgd_all(:,idb,:);
  end
  % fix to reduce funny bottom track dimension
  id = [1:length(data.tim)];
  data.hb = median(btd(id,1:4),2)';
  data.hb4 = btd(id,1:4)';
  data.ub = btd(id,5)';
  data.vb = btd(id,6)';
  data.wb = btd(id,7)';
  data.eb = btd(id,8)';
  data.tsd_m = ead_m;
  data.cmd_m = cmd_m;
  for n = 1:4
    [dummy,ir] = min(abs(data.cmd_m(:,n)-max(data.cmd_m(1,:))*0.3));
    params.dn_range(n) = data.zd(ir);
  end


  %
  % parse instrument serial numbers
  %
  data.serial_cpu_d = fd(f.serial);
  values.instid(1) = prod(data.serial_cpu_d+1) + sum(data.serial_cpu_d);
  data.serial_inst_d = sum( fd(f.serial_instrument).*[1,256,256^2,256^3]);
  values.inst_serial(1) = data.serial_inst_d;
  if abs(values.inst_serial(1))>1e6
    disp('>   Can not interpret serial number in raw data file.')
    disp('>   This might happen because of an old firmware version')
  else
    params.down_sn = values.inst_serial(1);
  end
  if ~isnan(values.inst_serial(1))
    disp(['    Down looking instrument ID : ',int2str(values.inst_serial(1))])
  end
  
end

if data.btrk_used == 1
  good = find(isfinite(data.wb));
  disp(['    Found ',int2str(length(good)),' finite RDI bottom velocities'])
end


%
% discard dummy error velocities
%
bad = find(data.eb==-32.768);
if ~isempty(bad)
  disp(['    Found ',int2str(length(bad)),...
	' NaN bottom error velocities and discarded them'])
  data.eb(bad) = nan;
end


%
% check for 3-beam solution
%
if isfield(data,'zu')
  ind = [1:length(data.zu)];
  jok_up = cumprod(size(find(~isnan(data.w(ind,:)))));
  j_up = cumprod(size(find(isnan(data.e(ind,:)) & ~isnan(data.w(ind,:)))));
  ind = length(data.zu)+[1:length(data.zd)];
  if j_up/jok_up > 0.2
    warn=(['>   Detected  ',int2str(j_up*100/jok_up),' %  3 BEAM solutions up-looking']);
    messages.warn = strvcat(messages.warn,warn);
  else
    warn=(['    Detected  ',int2str(j_up*100/jok_up),' %  3 BEAM solutions up-looking']);
  end
  disp(warn)
else
  ind = [1:length(data.zd)];
end 
jok_dn = cumprod(size(find(~isnan(data.w(ind,:)))));
j_dn = cumprod(size(find(isnan(data.e(ind,:)) & ~isnan(data.w(ind,:)))));
if j_dn/jok_dn > 0.2
  warn=(['>   Detected  ',int2str(j_dn*100/jok_dn),' %  3 BEAM solutions down-looking']);
  messages.warn = strvcat(messages.warn,warn);
else
  warn=(['    Detected  ',int2str(j_dn*100/jok_dn),' %  3 BEAM solutions down-looking']);
end
disp(warn)


%
% apply error velocity threshold
%
data.raw_u = data.u;
data.raw_v = data.v;
data.raw_w = data.w;
data.raw_e = data.e;
j = find(abs(data.e) > elim);
disp(['    Removed ',int2str(length(j)),...
	' values because of high error velocity'])
data.u(j) = NaN;
data.v(j) = NaN;
data.w(j) = NaN;
j = find(abs(data.eb) > elim);
disp(['    Removed ',int2str(length(j)),...
	' bottom values because of high error velocity'])
data.ub(j) = NaN;
data.vb(j) = NaN;
data.wb(j) = NaN;


%
% save all bottom track data together and remove bad data
%
data.bvel = [data.ub',data.vb',data.wb',data.eb'];
ind = find(data.bvel<-30);
data.bvel(ind) = NaN;
params.btrk_used = data.btrk_used;
if isfield(data,'hs')==1
  data.hsurf = data.hs;
end


%
% stuff to make output look like Martin's d structure
%
data.rw = data.w;
data.re = data.e;
data.temp = data.t;
data.weight = data.cm;
data.weight = data.weight./nmedian(nmax(data.weight));
data.tilt = real(asin(sqrt(sin_d(data.pit(1,:)).^2 + sin_d(data.rol(1,:)).^2)))/pi*180;
rold = mean(abs(diff([0,data.rol(1,:);data.rol(1,:),0]'))');
pitd = mean(abs(diff([0,data.pit(1,:);data.pit(1,:),0]'))');
data.tiltd = sqrt(rold.^2+pitd.^2);
data.hbot = data.hb;
data.hbot4 = data.hb4;
data.izd = [1:length(data.zd)];
params.bins_d = data.izd;
if values.up==1
  data.izu = fliplr([1:length(data.zu)]);
  params.bins_u = [data.izu,data.izd*0];
  params.bins_d = [data.izu*0,data.izd];
  data.izd = data.izd + length(data.izu);
  params.all_trusted_i = findany(params.bins_u+params.bins_d,params.trusted_i);
else
  data.izu = [];
  data.zu = [];
  params.bins_d = [data.izd];
  params.all_trusted_i = findany(params.bins_d,params.trusted_i);
end
data.soundc = 0;
data.ru = data.u;
data.rv = data.v;
for j=1:size(data.xmv,1)
  values.xmc(j) = meanmediannan(data.xmc(j,:),size(data.xmc,2)/4);
  values.xmv(j) = meanmediannan(data.xmv(j,:),size(data.xmv,2)/4);
  values.tint(j) = meanmediannan(data.tint(j,:),size(data.tint,2)/4);
end


%
% don't know what this is
% n_iz
%
n_iz = min([length(data.izd),6]);
sw = nstd(data.w(data.izd(2:n_iz),:));
ii = (sw>0);   % removed a find GK, 20.05.2011
sw = nmedian(sw(ii));
data.down.Single_Ping_Err = sw/tan(data.down.Beam_angle*pi/180)*...
                       sqrt(data.down.Pings_per_Ensemble);
values.down_beam_angle = data.down.Beam_angle;
values.nping_total = data.npng_d*data.nens_d;
if values.up==1
  n_iz = min([length(data.izu),6]);
  sw = nstd(data.w(data.izu(2:n_iz),:));
  ii = (sw>0);   % removed a find GK, 20.05.2011
  sw = nmedian(sw(ii));
  data.up.Single_Ping_Err = sw/tan(data.up.Beam_angle*pi/180)*...
                       sqrt(data.up.Pings_per_Ensemble);
  values.up_beam_angle = data.up.Beam_angle;
  values.nping_total(2) = data.npng_u*data.nens_u;
end
ii = isfinite(data.tim(1,:));
data = cutstruct(data,ii);


%
% check time for NB, which stores only day of year not year
%
for n=1:size(data.tim,1)
  time_greg = gregoria(data.tim(n,:));
  if time_greg(1)==1900
    if isfield(params,'correct_year')
      time_greg(:,1) = params.correct_year;
    else
      disp('>   Found Narrowband with wrong year.')
      error('Please set params.correct_year in cast_params.m')
    end
    disp('>   Narrowband : found year 1900, correcting to year of nav-data')
  end
  data.tim(n,:) = julian(time_greg)';
end


%
% get start and end times of data set
% these times will be used to extract other data
% such as navigation and shipboard ADCP
%
% the only way to override them is to give 
% start and end time explicitly in positions.dat
%
values.start_time = min(data.tim(1,:));
values.end_time = max(data.tim(1,:));


%
% clear pressure records, if requested
%
if params.clear_ladcp_pressure==1
  data.pres = data.pres*0;
  data.presstd = data.presstd*0;
end


%-------------------------------------------------------------------------------
function res = isbb(fid)
%ISBB True if broad-band ADCP.

% check header and data source identification bytes
hid = 127;
sid = 127;
id = fread(fid,2,'uint8');
if length(id)<2
  error('ISBB: ****** can not read file id *****')
else
  res = (id(1) == hid & id(2) == sid);
end

% rewind file
fseek(fid,0,'bof');


%-------------------------------------------------------------------------------



function [vele]=b2earth(velb,v,p)
% 
% convert beam ADCP data to earth velocities
%
% input velb:  beam coordinates
%          v:  attitude vector
%          p:  ADCP information
% 
% output vele: earth coordinates
%
% hard wired for LADCP systems
% M. Visbeck  Jan 2004

if p.Coordinates~=0
  disp('>   Data are not in beam coordinates ! Not rotating.')
  vele = velb;
  return
end

%p.use_tilt = 1;			% removed usage GK
%p = setdefv(p,'use_heading',1);	% removed usage GK
p = setdefv(p,'use_binremap',1);	% this seems to be used
p = setdefv(p,'beams_up',p.Up);
p.sensor_config = 1;
p.convex = 1;

% head = value for the ensemble, in degrees
% pitch = value for the ensemble, in degrees
% roll = value for the ensemble, in degrees
% earth = resultant matrix
% beams_up = 1 for upward looking (default), 0 for downward 
% convex = 1 for convex (default) or 0 for concave xducers
% sensor_config = sensor configuration, default = 1, fixed
% BeamAngle = default = 20 degrees

% Written by Marinna Martini for the 
% U.S. Geological Survey
% Branch of Atlantic Marine Geology
% Thanks to Al Pluddeman at WHOI for helping to identify the 
% tougher bugs in developing this algorithm


% precompute some constants
d2r=pi/180; % conversion from degrees to radians
C30=cos_d(p.Beam_angle);
S30=sin_d(p.Beam_angle);

if p.beams_up == 1, % for upward looking
  ZSG = [+1, -1, +1, -1];
else % for downward looking
  ZSG = [+1, -1, -1, +1];
end


% size of problem
%nb = size(velb,2);
ne = size(velb,1);

% disp([' converting ',int2str(ne),' profiles from beam- to earth coordinates'])

vele = velb*nan;

%big loop over profiles
for ii=1:ne

  roll = v(ii,1,3);
  pitch = v(ii,1,2);
  head = v(ii,1,4);
  beam = squeeze(velb(ii,:,:));

  % Step 1 - determine rotation angles from sensor readings
  % fixed sensor case
  % make sure everything is expressed in radians for MATLAB
  RR = roll.*d2r;
  KA = sqrt(1.0 - (sin_d(pitch).*sin_d(roll)).^2);
  PP = asin(sin_d(pitch).*cos_d(roll)./KA);
  HH = head.*d2r;

  % Step 2 - calculate trig functions and scaling factors
%  if p.use_tilt
    CP = cos(PP); 
    CR = cos(RR); 
    SP = sin(PP); 
    SR = sin(RR); 
%  else
%    CP = 1; 
%    CR = 1;
%    SP = 0; 
%    SR = 0;
%  end

%  if p.use_heading
    CH = cos(HH); 
    SH = sin(HH);
%  else
%    CH = 1; 
%    SH = 0; 
%  end

  % fixed sensor case
  M(1) = -SR.*CP;
  M(2) = SP;
  M(3) = CP.*CR;

  % compute scale factor for each beam to transform depths
  % in a tilted frame to depths in a fixed frame
  % RDI version of code
  %SC(1) = C30./(M(3).*C30 + ZSG(1).*M(1).*S30);
  %SC(2) = C30./(M(3).*C30 + ZSG(2).*M(1).*S30);
  %SC(3) = C30./(M(3).*C30 + ZSG(3).*M(2).*S30);
  %SC(4) = C30./(M(3).*C30 + ZSG(4).*M(2).*S30);
  % Pluddeman version of code 
  % changed for his difference in beam angle convention
  SC(1) = (M(3).*C30 + ZSG(1).*M(1).*S30);
  SC(2) = (M(3).*C30 + ZSG(2).*M(1).*S30);
  SC(3) = (M(3).*C30 + ZSG(3).*M(2).*S30);
  SC(4) = (M(3).*C30 + ZSG(4).*M(2).*S30);


  % form the transducer to instrument coordinate system
  % scaling constants
  % RDI version
  %VXS = (C*100.0)/(xfreq*4*S30);
  %VYS = (C*100.0)/(xfreq*4*S30);
  %VZS = (C*100.0)/(xfreq*8*C30);
  %VES = (C*100.0)/(xfreq*8);

  % form the transducer to instrument coordinate system
  % scaling constants
  % original Al Plueddeman version where theta is the 
  % beam angle from the horizontal
  % sthet0=sin(d2r*(90-BeamAngle));
  % cthet0=cos(d2r*(90-BeamAngle));
  % VXS = VYS = SSCOR / (2.0*cthet0);
  % VZS = VES = SSCOR / (2.0*sthet0);
  % correct for speed of sound using ADCP sound speed
  % based on thermistor measurements, where 1500 was the
  % assumed sound speed.
  %SSCOR = ssnd/1500;
  %SSCOR = ssnd/ECssnd;

  SSCOR = 1;
  % my version of Al's scaling constant, using RDI's
  % convention for theta as beam angle from the vertical
  VXS = SSCOR/(2.0*S30);
  VYS = VXS;
  VZS = SSCOR/(4.0*C30);
  VES = VZS;

  [NBINS] = size(beam);
  earth = zeros(size(beam));
  J = zeros(1,4);

  for IB=1:NBINS,
    % Step 3:  correct depth cell index for pitch and roll
    for i=1:4, 
      if p.use_binremap
        J(i) = fix(IB.*SC(i)+0.5); 
      else
        J(i) = IB;
      end
    end

    % Step 4:  ADCP coordinate velocity components
    if all(J > 0) && all(J <= NBINS),
      if any(isnan(beam(IB,:))),
        earth(IB,:)=ones(size(beam(IB,:))).*NaN;
      else
        %  if p.convex ,
        if p.beams_up ,
          % for upward looking convex
          VX = VXS.*(-beam(J(1),1)+beam(J(2),2));
          VY = VYS.*(-beam(J(3),3)+beam(J(4),4));
          VZ = VZS.*(-beam(J(1),1)-beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
          VE = VES.*(+beam(J(1),1)+beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
        else
          % for downward looking convex
          VX = VXS.*(+beam(J(1),1)-beam(J(2),2));
          VY = VYS.*(-beam(J(3),3)+beam(J(4),4));
          VZ = VZS.*(+beam(J(1),1)+beam(J(2),2)+beam(J(3),3)+beam(J(4),4));
          VE = VES.*(+beam(J(1),1)+beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
        end
        % else,
        %  if p.beams_up ,
        %   % for upward looking concave
        %   VX = VXS.*(+beam(J(1),1)+beam(J(2),2));
        %   VY = VYS.*(+beam(J(3),3)+beam(J(4),4));
        %   VZ = VZS.*(-beam(J(1),1)-beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
        %   VE = VES.*(+beam(J(1),1)+beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
        %  else
        %   % for downward looking concave
        %   VX = VXS.*(-beam(J(1),1)+beam(J(2),2));
        %   VY = VYS.*(+beam(J(3),3)-beam(J(4),4));
        %   VZ = VZS.*(+beam(J(1),1)+beam(J(2),2)+beam(J(3),3)+beam(J(4),4));
        %   VE = VES.*(+beam(J(1),1)+beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
        %  end  
        % end

        % Step 5: convert to earth coodinates
        VXE =  VX.*(CH*CR + SH*SR*SP) + VY.*SH.*CP + VZ.*(CH*SR - SH*CR*SP);
        VYE = -VX.*(SH*CR - CH*SR*SP) + VY.*CH.*CP - VZ.*(SH*SR + CH*SP*CR);
        VZE = -VX.*(SR*CP)            + VY.*SP     + VZ.*(CP*CR);
        earth(IB,:) = [VXE, VYE, VZE, VE];
      end % end of if any(isnan(beam(IB,:))),
    else
      earth(IB,:) = ones(size(beam(IB,:))).*NaN;
    end % end of if all(J > 0) && all(J < NBINS),
  end % end of IB = 1:NBINS

  % save results
  vele(ii,:,:) = earth;

end % Big Loop
