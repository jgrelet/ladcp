function [values] = prepctdprof(params, values)
% function [values] = prepctdprof(stn,values)
%
% prepare CTD profile for LADCP
% we need an array 'ctdprof' containing the 3 columns
% pressure in dbar    in situ temperature in degrees C    salinity in psu
%
%
% THIS FILE IS CRUISE SPECIFIC
% 
% to create a file for your own cruise, modify this file
%
% the data should typically be a profile in 1dbar or 1m steps
% (a lower resolution of down to 10dbar or 10m might be sufficient)
% it will be used to calculate depth dependent sound speed corrections
%
% If such data is not available, a sound speed profile will be
% derived from the ADCP's temperature sensor, the integrated
% vertical velocity and a constant salinity.

% G.Krahmann, IFM-GEOMAR, Aug 2005

% if you do no have CTD profile data to be used in the 
% LADCP processing, uncomment the next two line, otherwise edit the following

% disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepctdprof.m !')
% pause
% return

% first copy CTD profile to the raw CTD data directory
% data/raw_ctd
% this data could e.g. be coming from a mounted disk like in
% the example below

% eval(['!copy z:\IFM_Leg4\CTD\for_use_uncalibrated\dATA4_',int2str0(stn,3),...
%   	'_1dbar.cnv data\raw_ctdprof'])
global pathFile;

fname = strcat(pathFile,'/data-processing/CTD/data/cnv/fr26',...
  params.ladcp_station_name, '.cnv');

fprintf('    PREPCTDPROF:');

if exist(fname,'file')
  fprintf('  read %s\n', fname);
  copyfile(fname,['data/raw_ctdprof/',...
    params.ladcp_station_name, '.ctd']);
  
  % load the data and convert to standard format
  %
  % in this example
  % we extract the PTS columns and get position and time data from the header
  %
  % you might have to convert depth to pressure in dbar
  % and/or conductivity to salinity
  [hdr,data] = read_sbe_cnv(['data/raw_ctdprof/',...
    params.ladcp_station_name, '.ctd']);
  
  ctdprof = [data.p,data.t_pri,data.s_pri];
  values.ctd_time = julian(hdr.nmea_utc);
  values.ctd_lat = hdr.nmea_lat;
  values.ctd_lon = hdr.nmea_lon;
  
  % the pressure data in the example had some spikes which
  % could be removed by the following
  % If your data quality is already good, you won't need the
  % following lines
  %good = find(ctdprof(:,3)>1);
  %ctdprof = ctdprof(good,:);
  
else
  error('ladcp:prepctdprof', 'file not exist: %s\nCheck the configuration of cruise_param.m file\n', fname);
end


% store data at the standard location
save6(['data/ctdprof/ctdprof',...
  params.ladcp_station_name], 'ctdprof')

% save filename
file = ['data/ctdprof/ctdprof', params.ladcp_station_name];
