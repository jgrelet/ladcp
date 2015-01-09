function [values] = prepctdtime(params, values)
% function [values] = prepctdtime(stn,values)
%
% prepare CTD data against time for LADCP
%
% we need a vector 'timctd' with the time of the CTD in Julian days
% and an array 'data' with the 3 columns
% pressure in dbar    in situ temperature in degrees C    salinity in psu
%
% THIS FILE IS CRUISE SPECIFIC
%
% to create a file for your own cruise, modify this file
%
% the data should typically be the data recorded during a CTD
% cast in about 1 second steps
% it will be used to calculate the depth of the LADCP system

% G.Krahmann, IFM-GEOMAR, Aug 2005

% if you do no have CTD time data to be used in the
% LADCP processing, uncomment the next two lines, otherwise edit the following


%
% first copy CTD time data to raw CTD data directory
% data/raw_ctdtime
% this data could e.g. be coming from a mounted disk like in
% the example below
%
% uncomment the following and MODIFY it,
% if you have the raw data stored elsewher
% and want to copy it to the ADCP processing
%
global pathFile

if ~exist('data/raw_ctdtime/', 'dir')
  mkdir('data/raw_ctdtime');
end

fname = strcat(pathFile,'/data-processing/CTD/data/ladcp/fr24',...
  params.ladcp_station_name, '_ladcp.cnv');

fprintf('    PREPCTDTIME:');

if exist(fname,'file')
  fprintf('  read %s\n', fname);
  copyfile(fname,['data/raw_ctdtime/',...
    params.ladcp_station_name, '.cnv']);
  
  
  % station three was an interrupted CTD file. Thus there is no
  % full CTD-time and NAV (from CTD) data available
  
  
  %
  % load this data and convert to standard format
  % we need the data:
  %
  % time in decimal julian days ( January 1, 2000 = 2451545 )
  % pressure in dbar
  % in situ temperature in degrees C
  % salinity in psu
  %
  % time is stored as acolumn vector in 'timctd'
  % the other variables as columns PTS in the array 'data'
  %
  %
  % in this example
  % we skip the header of the file and extract the PTS columns
  % into 'data' and the time vector into 'timctd'
  %
  % you might have to convert depth to pressure in dbar
  % and/or conductivity to salinity
  %
  % and you will have to make sure that the time is stored in Julian days
  %
  % in this example we add the julian day January 0 of the year of the
  % cast to the time stored in the file
  % THIS IS APPROPRIATE FOR SEABIRD CNV FILES that contain the
  % variable 'time in julian days'
  %
  [hdr,ctd] = read_sbe_cnv(['data/raw_ctdtime/',...
    params.ladcp_station_name, '.cnv']);
  timctd = julian(datevec(ctd.datenum));
  data = [ctd.p,ctd.t_pri,ctd.s_pri];
  
  % remove NaN values
  ind = find( isnan(timctd) | any(isnan(data),2) );
  timctd(ind) = [];
  data(ind,:) = [];
  
  %
  % the pressure data on one of our cruises had some spikes which
  % could be removed by the following
  % If your data quality is already good, you won't need the
  % following lines. If it is bad you will need do create your
  % own despiking.
  %
  %good = find(data(:,3)>1);
  %data = data(good,:);
  %timctd = timctd(good);
  
  %
  % In our example the CTD cast recording began a bit before the
  % actual down movement of the CTD. We want only the real
  % down and uptrace of the cast. A good part of this will be
  % done again in the main processing, but sometimes those
  % routines failed and it was simple enough to do here.
  %
  % The extraction of this part is sometimes tricky. You might
  % have to 'invent' your own methods to make sure that only the
  % cast is extracted.
  %
  % Here we cut start and end of profiles to when CTD is under water
  % This is done by taking the maximum pressure
  % finding the two values closest to half of this pressure
  % on the up and the down casts and go towards the
  % surface on up and down casts until salinity is lower than 10.0
  %
  if ~isfield(values,'ctd_time')
    values.ctd_time = nmedian(timctd);
  end
  if ~isfield(values,'start_time')
    values.start_time = timctd(1);
  end
  if ~isfield(values,'end_time')
    values.end_time   = timctd(end);
  end
  
  
  %
  % The following might not be necessary. But we needed it
  % in some cases.
  % Interpolate to a regular time stepping.
  %
  % uncomment the following only, if you experience problems with the
  % CTD interpolation in the merging part of the processing
  %
  %min_t = min(timctd);
  %max_t = max(timctd);
  %delta_t = median(diff(timctd));
  %data = interp1q(timctd,data,[min_t:delta_t:max_t]');
  %timctd = [min_t:delta_t:max_t]';
  %disp(sprintf('  interpolated to %d CTD scans; delta_t = %.2f seconds',...
  %	length(timctd),median(diff(timctd))*24*3600));
  %
  
else
  error('ladcp:prepctdtime', 'file not exist: %s\nCheck the configuration of prepctdtime.m and cruise_param.m file\n', fname);
end

%
% store data in the standard location
%
save6(['data/ctdtime/ctdtime',...
  params.ladcp_station_name], 'timctd', 'data')
