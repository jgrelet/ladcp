function prepctdtime(stn)
% function prepctdtime(stn)
%
% prepare CTD data against time for LADCP
%
% this file is cruise specific for NBP 04-08
% 
% to create a file for your own cruise, copy this file and
% modify where necessary

% copy CTD time data from public drives to raw CTD data
% on ladcp computer
eval(['!cp /nbp/science/NBP0408/ladcp/ctd/p408',int2str0(stn,3),...
  	'ladcp.cnv data/raw_ctd/p408',int2str0(stn,3),'.tim'])

% load this data and convert to standard format
[hdr,data] = hdrload(['data/raw_ctd/p408',int2str0(stn,3),'.tim']);
timctd = data(:,7) + julian([2004,1,0,0,0,0]);
data = data(:,[4,2,8]);

% catch spikes
good = find(data(:,3)>1);
data = data(good,:);
timctd = timctd(good);

% cut start and end of profiles to sea surface
% this is done by taking the maximum pressure
% finding the two values closest to half of this pressure
% on the up and the down casts and go towards the
% surface on up and down until one reaches 2dbar or
% the last value
[pmax,indmax] = maxnan(data(:,1));
[dummy,mid1] = min( abs(pmax/2-data(1:indmax,1)) );
[dummy,mid2] = min( abs(pmax/2-data(indmax+1:end,1)) );
mid2 = mid2+indmax-1;
inds = max( find( data(1:mid1,1)< 2 ) );
if isempty(inds)
  inds = 1;
end
mid2
inde = min( find( data(mid2:end,1)< 2 ) ) + mid2-1;
if isempty(inde)
  inde = size(data,1);
end
data = data(inds:inde,:);
timctd = timctd(inds:inde);

% interpolate to regular time series
min_t = min(timctd);
max_t = max(timctd);
delta_t = median(diff(timctd));
data = interp1q(timctd,data,[min_t:delta_t:max_t]');
timctd = [min_t:delta_t:max_t]';
disp(sprintf(' interpolated to %d CTD scans; delta_t = %.2f seconds',...
	length(timctd),median(diff(timctd))*24*3600));

if 0
    % check for spikes in pressure data
    dtadcp = mediannan( diff(data.time_jul) );
    dtctd = mediannan( diff(data.ctdtime_time),length(data.ctdtime_time)/5 );
    disp(['  median CTD time difference ',num2str(dtctd*24*3600),' s'])
    % fix jitter in CTD time base
    timctd2 = cumsum( round((diff(timctd))/dtctd)*dtctd );
    timctd2 = [timctd(1); timctd(1)+timctd2];
    ibad = [1:10];
    sbad = 0;
    nn = 0;
    while length(ibad)>9 & nn<10
      z = -sw_dpth(data(:,1),lat);
      wctd = [-diff(z)./(diff(timctd2)*24*3600)];
      ibad = find(abs(wctd)>(3*std(wctd)));
      ibad = [ibad; ibad+1];
      data.ctdtime_data(ibad,:) = [];
      data.ctdtime_time(ibad) = [];
      timctd2(ibad,:) = [];
      sbad = sbad + length(ibad);
      nn=nn+1;
    end
    disp(['  removed ',int2str(sbad),' pressure spikes'])
    z = -sw_dpth(data(:,1),lat);
end

% store data in the standard location
eval(['save data/ctd/ctdtime',int2str0(stn,3),' timctd data'])
