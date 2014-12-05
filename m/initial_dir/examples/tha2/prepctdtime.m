function prepctdtime(stn)
% function prepctdtime(stn)
%
% prepare CTD data against time for LADCP
%
% we need a vector 'timctd' with the time of the CTD in Julian days
% and an array 'data' with the 3 columns
% pressure in dbar    in situ temperature in degrees C    salinity in psu
%
% THIS FILE IS CRUISE SPECIFIC
%
% THALASSA 2
%
% to create a file for your own cruise, modify this file
%
% the data should typically be the data recorded during a CTD
% cast in 1 second steps
% it will be used to calculate the depth of the LADCP system

% G.Krahmann, IFM-GEOMAR, Aug 2005

% if you do no have CTD time data to be used in the
% LADCP processing, uncomment the next line

%return


% first copy CTD time data to raw CTD data directory
% data/raw_ctdtime
% this data could e.g. be coming from a mounted disk like in
% the example below
eval(['!cp /dosd/dvd/ladcp/th/data/raw_ctd/tha2_',int2str0(stn,3),...
  	'.prs data/raw_ctdtime'])

% load this data and convert to standard format
%
% in this example
% we skip the header of the file and extract the PTS columns
% into 'data' and the time vector into 'timctd'
%
% you might have to convert depth to pressure in dbar
% and/or conductivity to salinity
%
% and you will have to make sure that the time is stored in Julian days
[hdr,data] = hdrload(['data/raw_ctdtime/tha2_',int2str0(stn,3),'.prs']);
timctd = data(:,1) + julian([2005,1,0,0,0,0]);
data = data(:,[2,3,5]);

% the pressure data on one of our cruises had some spikes which
% could be removed by the following
% If your data quality is already good, you won't need the
% following lines. If it is bad you will need do create your
% own despiking.
good = find(data(:,3)>1);
data = data(good,:);
timctd = timctd(good);

% In our example the CTD cast recording began a bit before the
% actual down movement of the CTD. We want only the real
% down and uptrace of the cast.
%
% The extraction of this part is sometimes tricky. You might
% have to 'invent' your own methods to make sure that only the
% cast is extracted.
%
% Here we cut start and end of profiles to near sea surface.
% This is done by taking the maximum pressure
% finding the two values closest to half of this pressure
% on the up and the down casts and go towards the
% surface on up and down casts until one reaches either 2dbar or
% the last value
[pmax,indmax] = maxnan(data(:,1));
[dummy,mid1] = min( abs(pmax/2-data(1:indmax,1)) );
[dummy,mid2] = min( abs(pmax/2-data(indmax+1:end,1)) );
mid2 = mid2+indmax-1;
inds = max( find( data(1:mid1,1)< 2 ) );
if isempty(inds)
  inds = 1;
end
inde = min( find( data(mid2:end,1)< 2 ) ) + mid2-1;
if isempty(inde)
  inde = size(data,1);
end
data = data(inds:inde,:);
timctd = timctd(inds:inde);

% The following might not be necessary. But we needed it
% in some cases.
% Interpolate to a regular time stepping.
min_t = min(timctd);
max_t = max(timctd);
delta_t = median(diff(timctd));
data = interp1q(timctd,data,[min_t:delta_t:max_t]');
timctd = [min_t:delta_t:max_t]';
disp(sprintf('  interpolated to %d CTD scans; delta_t = %.2f seconds',...
	length(timctd),median(diff(timctd))*24*3600));

% store data in the standard location
eval(['save data/ctdtime/ctdtime',int2str0(stn,3),' timctd data'])
