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
% LADCP processing, uncomment the next two line, otherwise edit the following

% disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepctdprof.m !')
% pause
% return

% first copy CTD profile to the raw CTD data directory
% data/raw_ctd
% this data could e.g. be coming from a mounted disk like in
% the example below

fname = ['m:/PANDORA/data-processing/CTD/PN',int2str0(stn,5),'.txt'];
if ~exist(fname,'file')
   fname = ['m:/PANDORA/data-processing/CTD/pn',int2str0(stn,5),'.txt'];
end
if ~exist(fname,'file')
   return
end
copyfile(fname,['data/raw_ctdprof/',int2str0(stn,5),'.ctd']);

% load the data and convert to standard format
% 
% in this example 
% we extract the PTS columns and get position and time data from the header
%
% you might have to convert depth to pressure in dbar
% and/or conductivity to salinity
fid = fopen(['data/raw_ctdprof/',int2str0(stn,5),'.ctd'],'r');
% date
for nl=1:3
    ligne = fgets(fid);
end
values.ctd_time = julian(datevec(sscanf(ligne,' DATE:  %20c')));
% position
ligne = fgets(fid);
pos_lat = sscanf(ligne,' LATITUDE: %d %f %c  LONGITUDE: %*d %*f %*c',[3 1]);
values.ctd_lat = pos_lat(1) + pos_lat(2)/60;
if strcmp(pos_lat(3),'N')
   values.ctd_lat = -pos_lat(1) + pos_lat(2)/60;
else 
   values.ctd_lat = -pos_lat(1) - pos_lat(2)/60;
end
pos_lon = sscanf(ligne,' LATITUDE: %*d %*f %*c  LONGITUDE: %d %f %c',[3 1]);
if strcmp(pos_lat(3),'E')
   values.ctd_lon = -pos_lon(1) + pos_lon(2)/60;
else
   values.ctd_lon = -pos_lon(1) - pos_lon(2)/60;
end
% heading
for nl=1:3
    fgets(fid);
end
% data
ctdprof = fscanf(fid,'%f %f %f %*f %*f %*d',[3 Inf])';
% close file
fclose(fid);

% remove NaN values
ind = find(any(ctdprof==-99,2));
ctdprof(ind,:) = [];

% store data at the standard location
save6(['data/ctdprof/ctdprof',int2str0(stn,5)],'ctdprof')

% save filename
file = ['data/ctdprof/ctdprof',int2str0(stn,5)];