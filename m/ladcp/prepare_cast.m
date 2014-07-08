function [values] = prepare_cast(stn)
% function [values] = prepare_cast(stn)
%
% prepare input data (CTD, LADCP, SADCP, NAV) so that it is digestible
% by the LADCP processing routines
%
% input  :	stn     station number
%           values  LADCP structure containing various values
%
% output :  values  LADCP structure containing various values

% this is to separate the processing from cruise dependent input data
% variations
%
% version 0.4	last change 07.2008

% G.Krahmann, LDEO

% modified the initial time and pos               GK, 06.03.2006  0.2-->0.3
% modified to save information of existing files  MV,    07.2008  0.3-->0.4


%
% general function info
%
disp(' ')
disp('PREPARE_CAST:  load raw data and store into MAT files')


values.initial = [];

% we need to prepare up to 5 different data sets
% CTD as a profile
% CTD as a timeseries
% NAV as a timeseries
% LADCP raw data in proper directories
% SADCP as one mat-file
%
% Unfortunately the formats and location in which this data is available
% varies from cruise to cruise. The poor user of these routines likely
% has to create a set of routines which bring the data into the desired standard
% format.
% Each data set has its own conversion routine. As there are typically
% some similarities between cruise setups such routines should be saved for 
% future (modified) use.
% Five loading routines must be exist and possibly be created in 'm/cruise'
% prepctdtime.m
% prepctdprof.m
% prepnav.m
% prepladcp.m
% prepsadcp.m
%
% Some examples can be found in the directory 'm/cruise/cruise_specific' .
% The examples will have filenames containing the cruise id. To use such a
% file copy it one directory up, remove the cruise id in the file name,
% and modify if necessary.

% In the following these routines are called
% and their output is specified

% PREPCTDPROF
% prepare CTD data as a profile
% this data is mainly used to calculate sound velocity profiles
% its data has to be in the form of a 3 column variable 'ctdprof'
% containing [pressure in dbar, in situ temperature in C, salinity in psu]
% 1 dbar steps is typically ok
% the data has to be stored in a mat-file named 'data/ctd/ctdprofSTN.mat'
% where STN is the 3-digit station number
%
% the structure 'values' will get position and time from the CTD
% data. This information will subsequently be used to cut the
% navigational and SADCP data (just to reduce the loading times and the
% necessary storage space)
if exist(['data/ctdprof/ctdprof',int2str0(stn,3),'.mat'],'file')
  disp(['    Found previously prepared CTD-PROFILE data.'])
else
  [values] = prepctdprof(stn,values);
end
if exist(['data/ctdprof/ctdprof',int2str0(stn,3),'.mat'],'file')
  values.ctdprofdata=1;
else
  values.ctdprofdata=0;
end

% prepare CTD data against time
% this data is used to determine the beginning and end of the LADCP cast
% its data has to be in the form of a 4 column variable 'ctdtime'
% containing [sec , pressure in dbar, in situ temperature in C, salinity in psu]
% 1 sec steps are necessary
% the data will be stored in a mat-file named 'data/ctd/ctdtimeSTN.mat'
% where STN is the 3-digit station number
if exist(['data/ctdtime/ctdtime',int2str0(stn,3),'.mat'],'file')
  disp(['    Found previously prepared CTD-TIME data.'])
else
  [values] = prepctdtime(stn,values);
end
if exist(['data/ctdtime/ctdtime',int2str0(stn,3),'.mat'],'file')
   values.ctdtimedata=1;
else
   values.ctdtimedata=0;
end

% define the start and end times which are to be used to
% cut navigational and SADCP data
%
% three ways are possible:
% 1     prescribed in set_cast_params
% 2     extracted from CTD profile (we use +/- 7 hours)
% 3     none (all navigational and SADCP data is stored, this might be slow)
if isfield(values,'start_time')     % this is the 'prescribed' case
  values.start_cut = values.start_time;
  values.end_cut = values.end_time;
elseif isfield(values,'ctd_time')   % this is the CTD info case
  values.start_cut = values.ctd_time-7/24;
  values.end_cut = values.ctd_time+7/24;
else
  values.start_cut = 2444240;       % 1980
  values.end_cut = 2488070;         % 2100
end  

% prepare navigational data against time
% this data is used to determine the location at the beginning and end of the LADCP cast
% its data has to be in the form of a 3 column variable 'nav'
% containing [sec , latitude in decimal degrees N, longitude in decimal degrees E]
% 1 sec steps are necessary
% the data will be stored in a mat-file named 'data/nav/navSTN.mat'
% where STN is the 3-digit station number
if exist(['data/nav/nav',int2str0(stn,3),'.mat'],'file')
  disp(['    Found previously prepared NAV data.'])
else
  prepnav(stn,values);
end
if exist(['data/nav/nav',int2str0(stn,3),'.mat'],'file')
   values.navdata=1;
else
   values.navdata=0;
end

% prepare LADCP data 
% this is only to put the raw data into the proper directories
% if you are using smart downloading routines this can be properly set up in there
% data files are supposed to be in
% 'data/raw_ladcp/STN'
% where STN is the 3-digit station number
prepladcp(stn);

% prepare SADCP data 
% this is to bring the SADCP data into the proper variables
% data file is supposed to be in
% 'data/sadcp/sadcp'
% the file has to contain
% tim_sadcp		vector containing the time in Julian days (get with julian.m)
% lat_sadcp		vector containing latitude in decimal degrees N
% lon_sadcp		vector containing longitude in decimal degrees E
% u_sadcp		matrix containing zonal velocities in m/s
% v_sadcp		matrix containing meridional velocities in m/s
% z_sadcp		vector containing the depths of the velocity data in m
if exist(['data/sadcp/sadcp',int2str0(stn,3),'.mat'],'file')
  disp('    Found previously prepared SADCP data.')
else
  prepsadcp(stn,values);
end
if exist(['data/sadcp/sadcp',int2str0(stn,3),'.mat'],'file')
  values.sadcpdata=1;
else  
  values.sadcpdata=0;
end

