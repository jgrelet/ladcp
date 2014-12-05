function [data,params,messages,values] = mergedata(data,params,messages,values)
% function [data,params,messages,values] =...
%	 mergedata(data,params,messages,values)
%
% merge navigational, CTD, and LADCP data
%
% This also corrects any time offsets between the different data streams
% and handles differences between BB and NB data files.
%
% input  :  data      - LADCP data structure
%           params    - LADCP parameter structure
%           messages  - array of warnings
%           values    - LADCP values structure
%
% output :  data      - modified LADCP data structure
%           params    - modified LADCP parameter structure
%           messages  - array of warnings
%           values    - LADCP values structure
%
% version 0.4	last change 17.01.2013

% G.Krahmann, LDEO Nov 2004

% added possibility to use SADCP nav data     GK, Jul 2007    0.1-->0.2
% renamed besttlag to calc_adcp_ctd_lag       GK, 20.05.2011  0.2-->0.3
% new parameter params.forced_adcp_ctd_lag    GK, 17.01.2013  0.3-->0.4

% Merging LADCP data with other data is always tricky since the internal
% LADCP and the LADCP-startup-computer clocks tend to drift. If there is
% a GPS-based time serving computer on the ship, it is extremely useful to use it
% to set the time of the LADCP-startup-computer to the proper time before
% each cast. This is automatically done by the LDEO's LINUX LADCP programs,
% if a time server has been identified.
% If the computer is standalone, you will have to set the time yourself.
% Setting the time can be missed or a wrong time can be entered. 
% params.timoff can correct such time offsets


%
% general function info:
%
disp(' ')
disp('MERGEDATA:  merge ADCP, navigational, and CTD-time data')




%
% just for display purposes here
%
disp('    ADCP data coverage:')
disp(['      ',datestr(gregoria(data.time_jul(1)),31),'  to  ',...
             datestr(gregoria(data.time_jul(end)),31)])
if ~isempty(data.tim_sadcp)
  disp('    SADCP data coverage:')
  disp(['      ',datestr(gregoria(data.tim_sadcp(1)),31),'  to  ',...
               datestr(gregoria(data.tim_sadcp(end)),31)])
end


%
% first try to merge the navigational information
% this assumes that the LADCP time base is correct
%
% if you know that it is wrong, you have to correct that in
% cast_params.m with   params.timoff
%
data.slat = data.time_jul*nan;
data.slon = data.time_jul*nan;
values.nav_source = 0;
if ~isempty(data.nav_data)
  disp('    Navigation data coverage:')
  disp(['      ',datestr(gregoria(data.nav_time(1)),31),'  to  ',...
               datestr(gregoria(data.nav_time(end)),31)])
  dummy_time = data.nav_time;
  ind = find( diff(dummy_time)>0 );
  ind = [ind;length(dummy_time)];
  lon = data.nav_data(:,2);
  lat = data.nav_data(:,1);
  data.slon = interp1(dummy_time(ind),lon(ind),data.time_jul','linear',nan)';
  data.slat = interp1(dummy_time(ind),lat(ind),data.time_jul','linear',nan)';
  if any(isnan(data.slon))
    if ~all(isnan(data.slon))
      bad = find(isnan(data.slon));
      disp(['>   Found ',int2str(length(bad)),...
	' LADCP data outside navigational data coverage'])
      if length(bad)>50
        disp('>   There might be a problem with your navigational data extraction')
        disp('>     Will continue anyway.')
      end
      if ~isempty(bad)
        good = find(~isnan(data.slon));
        disp('>   Replacing missing values by nearest nav value')
        data.slon(bad) = interp1(good,data.slon(good),bad,'nearest','extrap');
        data.slat(bad) = interp1(good,data.slat(good),bad,'nearest','extrap');
      end
    else
      disp('>   Navigational data interpolation leads to only NaN.')
      disp('>     Check your navigation files, for time overlap with')
      disp('>     the ADCP data and for NaNs.')   
    end
  end
  values.nav_source = 1;
end
if all(isnan(data.slat)) 
  if ~isnan(prod(params.poss))

    if isempty( params.time_start )
      disp('>   No start time given. Using first record of ADCP')
      data.slat = linspace(params.poss(1)+params.poss(2)/60,...
    	params.pose(1)+params.pose(2)/60,length(data.time_jul));
      data.slon = linspace(params.poss(3)+params.poss(4)/60,...
    	params.pose(3)+params.pose(4)/60,length(data.time_jul));
    else
      p1 = params.poss(1)+params.poss(2)/60;
      p2 = params.pose(1)+params.pose(2)/60;
      t1 = julian( params.time_start );
      t2 = julian( params.time_end );
      data.slat = p1 + (p2-p1) * (data.time_jul-t1)/(t2-t1);
      p1 = params.poss(3)+params.poss(4)/60;
      p2 = params.pose(3)+params.pose(4)/60;
      data.slon = p1 + (p2-p1) * (data.time_jul-t1)/(t2-t1);
    end
    disp('    Found no navigational data but manual input')
    disp('      assuming no navigation files exist')
    disp('      will solve with linear movement between start and end')
    values.nav_source = 2;
  elseif isfield(data,'lat_sadcp')
    if ~isempty(data.lat_sadcp)
      ind = find( diff(data.tim_sadcp)~=0 );
      data.slon = interp1(data.tim_sadcp(ind+1),data.lon_sadcp(ind+1),data.time_jul);
      data.slat = interp1(data.tim_sadcp(ind+1),data.lat_sadcp(ind+1),data.time_jul);
      disp('    Found no other navigational data than SADCP, please check your nav data sources')
      values.nav_source = 4;
    else
      data.slat = data.time_jul*NaN;
      data.slon = data.slat;
      disp('    Found no navigational data and no manual input')
      disp('      please check files or enter positions in cast_params.m')
      disp('      attempting to solve without position information')
      values.nav_source = 3;
    end
  else
    data.slat = data.time_jul*NaN;
    data.slon = data.slat;
    disp('    Found no navigational data and no manual input')
    disp('      please check files or enter positions in cast_params.m')
    disp('      attempting to solve without position information')
    values.nav_source = 3;
  end
end


%
% extract a position used for general purposes 
% such as magnetic deviation and pressure-depth conversion
%
values.start_pos = [data.slat(1),data.slon(1)];
values.end_pos = [data.slat(end),data.slon(end)];
values.lat = nmean([values.start_pos(1),values.end_pos(1)]);
values.lon = nmean([values.start_pos(2),values.end_pos(2)]);
if isnan(values.lat)
  values.lat = params.lat_for_calc;
  values.lon = params.lon_for_calc;
  values.start_pos = [values.lat,values.lon];
  values.end_pos = [values.lat,values.lon];
  disp('> SERIOUS WARNING !!!!')
  disp('>   Got no usable position information, am using 0,0 for calculations')
  disp('>     If no other navigational data is available, consider setting')
  disp('>     params.lat_for_calc  and  params.lon_for_calc')
  disp('>     to the approximate numbers.')
end


%
% if we have CTD time data we will now shift it onto
% the same time base as the LADCP data
%
if ~isempty(data.ctdtime_data)

  disp('    CTD data coverage:')
  disp(['      ',datestr(gregoria(data.ctdtime_time(1)),31),'  to  ',...
               datestr(gregoria(data.ctdtime_time(end)),31)])
  disp(' ')

  % catch out of range CTD-time data
  if min(data.ctdtime_time)>max(data.time_jul) |...
    max(data.ctdtime_time)<min(data.time_jul)

    disp('>   CTD timeseries does not overlap! WRONG STATION????')
    disp('>     Will use vertical velocity to get depth')
    disp('>     Moving CTD-time data into different variables')
    data.ctdtime_time_bad = data.ctdtime_time;
    data.ctdtime_data_bad = data.ctdtime_data;
    data.ctdtime_time = [];
    data.ctdtime_data = [];

  else  
    
    % create raw w_CTD time series
    lat = values.lat*ones(length(data.ctdtime_data(:,1)),1);
    z = -sw_dpth(data.ctdtime_data(:,1),lat);
    wctd = -diff(z)./(diff(data.ctdtime_time)*86400);

    % catch a problem with spikey pressure data
    bad = find(abs(wctd)>(3*std(wctd)));
    if length(bad)>10
      disp('>   Found pressure spikes in CTD data')
      disp('>     Recommend stronger filtering in prepctdtime.m')
    end

    % create a w_CTD time series with roughly the same
    % time stepping as the LADCP
    % this needs to be done to avoid aliasing by a more highly
    % resolved pressure time series which we differentiate and then
    % compare to the w of the LADCP
    dtadcp = nmedian( diff(data.time_jul) );
    dtctd = meanmediannan( diff(data.ctdtime_time),...
      length(data.ctdtime_time)/5 );
    nshift = max([1,round(dtadcp/dtctd)])+1;
    nshift2 = fix(nshift/2);
%    disp(['    Using ',int2str(nshift2),' step time base for W_ctd'])
    wctd = z+nan;
    i2 = [nshift : length(z)-nshift2];
    i1 = [1:length(i2)];

    % calculate W from CTD with LADCP-similar time step
    wctd(i1+nshift2) = -(z(i2)-z(i1))./...
      ((data.ctdtime_time(i2)-data.ctdtime_time(i1))*24*3600);

    % now we can safely
    % interpolate W onto ADCP time series
%    ind = find(diff(data.ctdtime_time)>0);
%    data.ctdtime_time = data.ctdtime_time(ind);
%    data.ctdtime_data = data.ctdtime_data(ind,:);
%    wctd = wctd(ind);
    ctdtime_time = data.ctdtime_time;
    ctdtime_time(1) = 0;
    ctdtime_time(end) = 1e32;
    warning off
    data_int = interp1(ctdtime_time,data.ctdtime_data,data.time_jul','nearest');
    data.wctd = interp1(ctdtime_time,wctd,data.time_jul','nearest')';
    warning on
    lat = values.lat*ones(length(data_int(:,1)),1);
    data.ctd_z = -sw_dpth(data_int(:,1),lat)';
    disp(['    CTD max depth : ',int2str(-min(data.ctd_z))])
    
    % prepare for time lag check
    dt = mean(diff(data.time_jul));
    w = meanmediannan(data.rw,2);
    w = replace(w, sum(isfinite(data.rw))<4, nan);

    % check for timelag
    %    dtctd = nmedian( diff(data.ctdtime_time) )
    % this one was a bit problematic
    dtctd = meanmediannan( diff(data.ctdtime_time), fix(length(data.ctdtime_time)/4) );

    %
    % make up array to check for lag
    %
    ctdtw = [data.ctdtime_time, wctd];
    adcptw = [data.time_jul', w'];

    if isfield(params,'forced_adcp_ctd_lag')
      lag = params.forced_adcp_ctd_lag;
      lagdt = -lag*dtctd;
      co = nan;
      disp(['    Forcing W lag at: ',int2str(lag),' CTD scans ~',...
           int2str(lagdt*24*3600)]);
    else
      [lag,co] = calc_adcp_ctd_lag(ctdtw,adcptw,params.ctdmaxlag,params.ctdmaxlagnp);
      lagdt = -lag*dtctd;
      disp(['    Best W lag at: ',int2str(lag),' CTD scans ~',...
           int2str(lagdt*24*3600),' seconds  corr:',num2str(co)]);
    end
    
    %
    % reinterpolate w from the CTD onto the LADCP timing
    %
    warning off
    data_int = interp1(data.ctdtime_time-lagdt,data.ctdtime_data,...
      data.time_jul','linear');
    data.wctd = interp1(data.ctdtime_time-lagdt,wctd,data.time_jul','nearest')';
    warning on
    data.z = -sw_dpth(data_int(:,1),lat)';
    data.ctdtime_data = data_int;
    data.ctdtime_time = data.time_jul;


    %
    % plot two small sections for visual control
    % ahhh, now I understand this one
    % since the lowest point often includes a bottle stop
    % it makes no sense to plot the lag corrected data there
    % it will mostly show noise
    % thus two parts of the up and down cast were selected
    % and plotted. A problem arises sometimes when 
    % multiple profiles are in one LADCP file. Then one can
    % get all NaN in the plotted ranges.  GK
    %
    % will change this to a hopefully less error prone algorithm GK
    good_range = find(~isnan(data.wctd));
    good_range = good_range([1,end]);
    nsamples = 24;
    i1 = [1:length(w)];
%    ii = [-nsamples:nsamples] + fix(length(i1)*0.3);
%    ii = [ii, [-nsamples:nsamples]+fix(length(i1)*0.7)];
    
    ii = [-nsamples:nsamples] + fix( good_range(1)+0.3*diff(good_range) );
    ii = [ii, [-nsamples:nsamples] + fix( good_range(2)-0.3*diff(good_range) )];
    ii = ii( find( ii>0 & ii<length(i1) ) );
    
    figure(2)  
    
    subplot(212)
    plot(w(i1(ii)),'-b')
    hold on
    plot(data.wctd(i1(ii)),'-r')
    axis tight
    ax=axis;
    plot((2*nsamples+1)*[1 1],ax(3:4),'-k')
    text(nsamples,mean(w(i1(ii(2*nsamples+1:end)))),'down cast (sample at 30% of cast)','horizontalalignment','center')
    text(3*nsamples,mean(w(i1(ii(1:2*nsamples)))),'up cast (sample at 70% of cast)','horizontalalignment','center')
    title(['best lag W: ',int2str(lag),' scans ~ ',...
             int2str(lagdt*86400),' sec.  corr.:',num2str(co)]);
    ylabel('W used for lag correlation')
    xlabel('sample scans ADCP (b) CTD (r)')
    grid
    
    streamer([params.name,'   Figure 8']);
    hgsave('tmp/8')

    % in difference to the case when the ADCP is shifted onto the
    % CTD data, this case here is always assumed good
    disp('    Adjusting CTD time to ADCP time and shift depth record ')
%    data.ctdtime_time = data.ctdtime_time-lagdt;
% don't understand why this was here   GK
    if lagdt*86400>10
      warn=['    Shifted CTD timeseries by ',int2str(lagdt*86400),' seconds '];
      disp(warn)
      messages.warn(size(messages.warn,1)+1,1:length(warn))=warn;
    end
  end
end


