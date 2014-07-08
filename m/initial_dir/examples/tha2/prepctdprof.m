function [file] = prepctdprof(stn)
% function [file] = prepctdprof(stn)
%
% prepare CTD profile for LADCP
% we need an array 'ctdprof' containing the 3 columns
% pressure in dbar    in situ temperature in degrees C    salinity in psu
%
% THIS FILE IS CRUISE SPECIFIC
%
% THALASSA 2
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

%return

% first copy CTD profile to the raw CTD data directory
% data/raw_ctd
% this data could e.g. be coming from a mounted disk like in
% the example below
eval(['!cp /dosd/dvd/ladcp/th/data/raw_ctd/tha2_',int2str0(stn,3),...
  	'.dat data/raw_ctdprof'])

% load the data and convert to standard format
% 
% in this example 
% we skip the header of the file and extract only the PTS columns
%
% you might have to convert depth to pressure in dbar
% and/or conductivity to salinity
[hdr,data] = hdrload(['data/raw_ctdprof/tha2_',int2str0(stn,3),'.dat']);
ctdprof = data(:,[1,2,4]);

% the pressure data in the example had some spikes which
% could be removed by the following
% If your data quality is already good, you won't need the
% following lines
good = find(ctdprof(:,3)>1);
ctdprof = ctdprof(good,:);

% store data at the standard location
eval(['save data/ctdprof/ctdprof',int2str0(stn,3),' ctdprof'])

% save filename
file = ['data/ctdprof/ctdprof',int2str0(stn,3)];
