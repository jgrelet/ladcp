function [values] = prepctdprof(stn,values)
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
% LADCP processing, uncomment the next line

%disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepctdprof.m !')
%return

if stn==21
  return
end

% first copy CTD profile to the raw CTD data directory
% data/raw_ctd
% this data could e.g. be coming from a mounted disk like in
% the example below

eval(['!cp /mnt/merian_public/public_wiss/MSM-052/CTD/',...
	'uncalibrated/dmsm052_',int2str0(stn,3),...
  	'_1dbar.cnv data/raw_ctdprof'])

% load the data and convert to standard format
% 
% in this example 
% we extract the PTS columns and get position and time data from the header
%
% you might have to convert depth to pressure in dbar
% and/or conductivity to salinity
[hdr,data] = hdrload(['data/raw_ctdprof/dmsm052_',int2str0(stn,3),'_1dbar.cnv']);
%ctdprof = data(:,[5,6,11]); % primary sensors
ctdprof = data(:,[5,7,17]); % secondary sensors

% extract position and time
for n=1:size(hdr,1)
  ind = findstr('NMEA UTC (Time)',hdr(n,:));
  strs = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
  if ~isempty(ind)
    [dummy,dat] = strtok(hdr(n,:),'=');
    values.ctd_time(2) = strmatch(dat(3:5),strs);
    values.ctd_time(1) = sscanf(dat(10:13),'%d')';
    values.ctd_time(3) = sscanf(dat(7:8),'%d')';
    values.ctd_time(4) = sscanf(dat(16:17),'%d')';
    values.ctd_time(5) = sscanf(dat(19:20),'%d')';
    values.ctd_time(6) = sscanf(dat(22:23),'%d')';
  end  
  ind = findstr('NMEA Latitude',hdr(n,:));
  if ~isempty(ind)
    [dummy,dat] = strtok(hdr(n,:),'=');
    dummy = sscanf(dat(3:10),'%d %f')';
    if findstr(dat,'N')
      values.ctd_lat = dummy(1)+dummy(2)/60;
    else
      values.ctd_lat = -dummy(1)-dummy(2)/60;
    end
  end  
  ind = findstr('NMEA Longitude',hdr(n,:));
  if ~isempty(ind)
    [dummy,dat] = strtok(hdr(n,:),'=');
    dummy = sscanf(dat(3:11),'%d %f')';
    if findstr(dat,'E')
      values.ctd_lat = dummy(1)+dummy(2)/60;
    else   
      values.ctd_lat = -dummy(1)-dummy(2)/60;
    end
  end  
end
values.ctd_time = julian(values.ctd_time(1),values.ctd_time(2),...
                         values.ctd_time(3),values.ctd_time(4)+...
                         values.ctd_time(5)/60+values.ctd_time(6)/3600);

% the pressure data in the example had some spikes which
% could be removed by the following

% If your data quality is already good, you won't need the
% following lines
good = find(ctdprof(:,3)>1);
ctdprof = ctdprof(good,:);

% store data at the standard location
if mversion>=70
  save(['data/ctdprof/ctdprof',int2str0(stn,3)],'ctdprof','-v6')
else
  save(['data/ctdprof/ctdprof',int2str0(stn,3)],'ctdprof')
end

% save filename
file = ['data/ctdprof/ctdprof',int2str0(stn,3)];
