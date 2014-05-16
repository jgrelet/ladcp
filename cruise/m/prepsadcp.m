function [] = prepsadcp(stn_str,values,pathFile)
% function [] = prepsadcp(stn,values)
%
% prepare Ship-ADCP data for LADCP processing
%
% we need the vectors 'tim_sadcp' , 'lon_sadcp' , 'lat_sadcp'
% and 'z_sadcp'
% and the arrays 'u_sadcp' and 'v_sadcp'
%
% THIS FILE IS CRUISE SPECIFIC
% 
% to create a file for your own cruise, modify this file
%
% the data should be the result of shipboard or later
% SADCP processing

% G.Krahmann, IFM-GEOMAR, Aug 2005


% if you do no have SADCP data to be used in the
% LADCP processing, uncomment the next two lines. Otherwise edit the following.

% disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepsadcp.m !')
% pause
% return

% first copy the SADCP files to the raw SADCP data directory
% data/raw_sadcp
% In our example
% this is data mounted via SMBMOUNT from //peale/adcp_home
%!copy z:\pos350\VMADCP\* data\raw_sadcp

% load this data and convert to standard format
%
% in this example we load the velocity and position/time files
% (the processing of the SADCP was also done in matlab which
% made the loading of files easy)
% and extract the necessary information
%
% again make sure that the time is in Julian days
% In the example the cruise was in 2004 and the processing
% stored only the day of the year not the actual year !!!
%  SADCP data file

% initialize data arrays
tim_sadcp = [];
lon_sadcp = [];
lat_sadcp = [];
z_sadcp = [];
u_sadcp = [];
v_sadcp = [];


flname = strcat(pathFile,'\data-processing\SADCP\CASCADE\ncc\PIRATA-FR24-ALL_osite.nc');

fprintf(1, '\nPREPSADCP:  load SADCP data\n');

% check that the SADCP file exists
if ~exist(char(flname),'file')
   sadcp = [];
   fprintf(1, 'cannot open SADCP file: %s\n', flname);
   return
end

% read the cascade file
nc = netcdf.open(flname,'NC_NOWRITE');
% time
varid =  netcdf.inqVarID(nc,'JULD_ADCP');
time = netcdf.getVar(nc,varid)';
tim_sadcp = julian(datevec(datenum(time + datenum(1950,1,1))));
% longitude
varid =  netcdf.inqVarID(nc,'LONGITUDE');
tab = netcdf.getVar(nc,varid);
% valid_max = netcdf.getAtt(nc,varid,'valid_max');
% valid_min = netcdf.getAtt(nc,varid,'valid_min');
% tab(tab>=valid_max | tab<=valid_min) = NaN;
lon_sadcp = tab;
% latitude
varid =  netcdf.inqVarID(nc,'LATITUDE');
tab = netcdf.getVar(nc,varid);
% valid_max = netcdf.getAtt(nc,varid,'valid_max');
% valid_min = netcdf.getAtt(nc,varid,'valid_min');
% tab(tab>=valid_max | tab<=valid_min) = NaN;
lat_sadcp = tab;
% depth
varid =  netcdf.inqVarID(nc,'DEPH');
z_sadcp = -netcdf.getVar(nc,varid);
% u
varid =  netcdf.inqVarID(nc,'UVEL_ADCP');
tab = netcdf.getVar(nc,varid);
valid_max = netcdf.getAtt(nc,varid,'valid_max');
valid_min = netcdf.getAtt(nc,varid,'valid_min');
tab(tab>=valid_max | tab<=valid_min) = NaN;
u_sadcp = tab;
% v
varid =  netcdf.inqVarID(nc,'VVEL_ADCP');
tab = netcdf.getVar(nc,varid);
valid_max = netcdf.getAtt(nc,varid,'valid_max');
valid_min = netcdf.getAtt(nc,varid,'valid_min');
tab(tab>=valid_max | tab<=valid_min) = NaN;
v_sadcp = tab;
% close file
netcdf.close(nc);


% restrict the data to the time of the cast
%good = find( tim_sadcp>values.start_time-0.1 & tim_sadcp<values.end_time+0.1);
good = find( tim_sadcp>values.start_cut-0.1 & tim_sadcp<values.end_cut+0.1 );
tim_sadcp = tim_sadcp(good);
lat_sadcp = lat_sadcp(good);
lon_sadcp = lon_sadcp(good);
u_sadcp = u_sadcp(:,good);
v_sadcp = v_sadcp(:,good);
z_sadcp = z_sadcp;

if ~isempty(tim_sadcp)
  fprintf(1, 'find sadcp data for %d ensemble\n', length(good));
end

% remove data below 1200 meters
if (0)
[zmax,imax] = max(abs(z_sadcp));
ind = find(z_sadcp<=-1200);
z_sadcp(ind) = [];
u_sadcp(ind,:) = [];
v_sadcp(ind,:) = [];
end

% store the data
save6(['data/sadcp/sadcp',stn_str],...
	'tim_sadcp','lon_sadcp','lat_sadcp','u_sadcp','v_sadcp','z_sadcp')
