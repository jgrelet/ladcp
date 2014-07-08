function [data,params,values,messages]=...
	misc_cut_profile(data,params,values,messages)
%function [data,params,values,messages]=...
%	misc_cut_profile(data,params,values,messages)
%
% Cut the raw LADCP data to only the profile part
%
% version 0.5	last change 05.03.2008

% Martin Visbeck  December 2002, LDEO
% G.Krahmann, IFM-GEOMAR, Jul 2005

% flagging only downlooker for 'below bottom'	GK, Sep 2007	0.1-->0.2
% chgd reaction to deep first CTD depth, this happens when 
% CTD cast was started late by accident		GK, Sep 2007	0.2-->0.3
% small bug when not reaching bottom		GK, Mar 2008	0.3-->0.4
% another small bug when not reaching bottom	GK, Mar 2008	0.4-->0.5
% save indices for detph inversion  MV Jul 2008

%
% general function info
%
disp(' ')
disp('MISC_CUT_PROFILE:  cut raw LADCP data to only the profile part')


[values.maxdepth,ibot] = max(-data.z);


%
% assign a depth to each bin
%
if values.up==1
  [izm1,izm] = meshgrid([fliplr(data.zu),-data.zd],data.z);
else
  [izm1,izm] = meshgrid(-data.zd,data.z);
end
izm1 = izm1';
if data.soundc==1
  % make sound speed correction for depth vector
  sc = meshgrid(data.ss./data.sv(1,:),data.izd);
  izm1(data.izd,:) = izm1(data.izd,:).*sc;
  if p.up==1
    sc = meshgrid(data.ss./data.sv(2,:),data.izu);
    izm1(data.izu,:) = izm1(data.izu,:).*sc;
  end
  disp('    Correcting bin length for sound speed')
end
% add the two parts
data.izm = izm'+izm1;


%
% prepare flag matrix
%
data.izmflag = data.ru*0;


%
% flag all data below bottom as bad
%
% blanking distance is one bin length
%
% flag only downlooking instrument
%
if ~isnan(params.zbottom)
  if values.up==1
    ind = find( data.izm(data.izd,:) < -params.zbottom+data.blen );
    dummy = data.izmflag(data.izd,:);
    dummy(ind) = NaN; 
    data.izmflag(data.izd,:) = dummy;
  else
    ind = find( data.izm < -params.zbottom+data.blen );
    data.izmflag(ind) = NaN; 
  end
end


%
% flag all data close to the surface as bad
%
if values.up==1
  ind = find( data.izm > -(data.zu(2)-data.zu(1))/2 );
  data.izmflag(ind) = NaN; 
end
  

%
% set weights for these values to NaN too
%
bad = find( isnan(data.izmflag) & isfinite(data.ru) );
if ~isempty(bad)
  disp(['    Removing ',int2str(length(bad)),...
	' values below bottom or above surface'])
  data.weight = data.weight + data.izmflag;
end


%
% compute pressure from depth
%
data.p = sw_pres(abs(data.z)',values.lat)';

if isfield(data,'wctd')==1
  data.wm = data.wctd;
else
  data.wm = -gradient(data.z)./gradient(data.time_jul*86400);
end


%
% warn, if start depth appears too deep
%
if data.z(1)<-50, 
  warn = ['>   Starting LADCP depth is ',int2str(data.z(1))];
  disp(warn)
  messages.warnp = strvcat(messages.warnp,warn);
end


%
% warn, if end depth appears too deep
%
if data.z(end)<-50, 
  warn = ['>   Ending LADCP depth is ',int2str(data.z(end))];
  disp(warn)
  messages.warnp = strvcat(messages.warnp,warn);
end


%
% cut raw data to only include profile
%
z_all = data.z_ladcp;
time_all = data.time_jul;
i1 = find(data.z(1:ibot)>-params.cut);

%
% the following 2 checks have been modified so that not the whole cast
% was used when the CTD depth data started below params.cut
% instead now all non-NaN data is used
%
if isempty(i1)
  i1 = find(~isnan(data.z(1:ibot)));
  i1 = i1(1);
  disp(['>   Found no early depth record shallower than ',...
	int2str(params.cut),' m'])
  disp(['>   Using first depth record of ',int2str(data.z(i1(1))),...
	' m at ensemble ',int2str(i1(1))])
end
i2 = find(data.z(ibot:end)>-params.cut)+ibot-1;
if isempty(i2)
  i2 = find(~isnan(data.z(ibot:end)))+ibot-1;
  i2 = i2(end);
  disp(['>   Found no late depth record shallower than ',...
	int2str(params.cut),' m'])
  disp(['>   Using first depth record of ',int2str(data.z(i2(end))),...
	' m at ensemble ',int2str(i2(end))])
end
ii = data.z*0;
ic = [i1(end):i2(1)];
ii(ic) = 1;
if (sum(ii)~=length(ii)) & params.cut>0 
  disp('    Removing data at beginning and end of cast')
  disp('      adjusting start and end time ')
  values.firstlastindex = values.firstlastindex(1)-1+[ic(1),ic(end)];
  data = cutstruct(data,ii);
  params.zpar([1 3]) = params.cut;
  subplot(221) 
  ax=axis; 
  plot([1 1]*ic(1),ax(3:4),'--k')
  subplot(222) 
  ax=axis; 
  plot([1 1]*ic(end),ax(3:4),'--k')
  subplot(212) 
  ax=axis; 
  plot([1 1]*ic(1),ax(3:4),'--k')
  plot([1 1]*ic(end),ax(3:4),'--k')
  axis(ax)
  
  params.time_start = gregoria(data.time_jul(1));
  params.time_end = gregoria(data.time_jul(end));

  % save start and end depth
  params.zpar(1) = -data.z(1);
  params.zpar(3) = -data.z(end);
end


%
% reset the time and position values
%
if ~all(diff(data.tim(1,:))>0)
  error('FOUND THAT TIME VECTOR IS NOT MONOTONOUS, CHECK DATA')
end
values.start_time = data.tim(1,1);
values.end_time = data.tim(1,end);
values.start_pos = [data.slat(1),data.slon(1)];
values.end_pos = [data.slat(end),data.slon(end)];
values.lat = mean([values.start_pos(1),values.end_pos(1)]);
values.lon = mean([values.start_pos(2),values.end_pos(2)]);


%
% compute ship drift
%
params.dt_profile = (values.end_time-values.start_time)*86400;

dlat = values.end_pos(1) - values.start_pos(1);
dlon = values.end_pos(2) - values.start_pos(2);

values.xdisp = dlon * cosd(values.lat) * 60.0 * 1852.0;
values.ydisp = dlat * 60.0 * 1852.0;

values.uship= values.xdisp / params.dt_profile;
values.vship= values.ydisp / params.dt_profile;

streamer([params.name,'   Figure 4']);
hgsave('tmp/4')


%
% add information of the cut profile to figure 8
%
figload(['tmp/8.fig'],2)

if ~isempty(data.ctdtime_data)
  subplot(211)
  plot((time_all-data.ctdtime_time(1))*24,-z_all,'b')
  hold on
  plot((data.ctdtime_time-data.ctdtime_time(1))*24,data.ctdtime_data(:,1),'r')
  title('Cut CTD profile    CTD (r)  LADCP (b)')
  xlabel('Time in hours')
  ylabel('CTD pressure [dBar]')
  set(gca,'ydir','reverse')

  ax = axis;
  plot((data.time_jul(1)-data.ctdtime_time(1))*24*[1,1], ax(3:4),'r--')
  plot((data.time_jul(end)-data.ctdtime_time(1))*24*[1,1], ax(3:4),'r--')
  axis tight
end

streamer([params.name,'   Figure 8']);
hgsave('tmp/8')

