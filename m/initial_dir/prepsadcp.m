function [] = prepsadcp(p,files,values)
% function [] = prepsadcp(p,files,values)
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

if isunix
    flname = '/m/CASSIOPEE\data-processing\SADCP\OS38\codas/cas38_ladcp.txt';
else
    flname = 'M:\CASSIOPEE\data-processing\SADCP\OS38\codas/cas38_ladcp.txt';  
end
fprintf('    PREPSADCP  :');

% check that the SADCP file exists
if ~exist(char(flname),'file') 
	error('ladcp:prepsadcp', 'Cannot read S-ADCP data from: %s\nCheck the file name and its location in prepadcp.m\n', flname);
   return
end

% read the ascii file
 fprintf('  read %s\n', flname);
fid = fopen(flname,'r');
data = fscanf(fid,'%f %f %f %f %f %f\n',[6 Inf])';
fclose(fid);

% replace missing value by NaN
data(data==999.999) = NaN;

% add reference year to time (+ 2min30)
data(:,1) = julian(datevec(datenum(p.correct_year,1,1)+data(:,1)));

% restrict the data to the time of the cast
igood = find( data(:,1)>values.start_cut-0.1 & data(:,1)<values.end_cut+0.1 );
data = data(igood,:);

% lon/lat
tim_sadcp = data(:,1);
lon_sadcp = data(:,2);
lat_sadcp = data(:,3);
NT = length(tim_sadcp);

% compute a unique vertical grid
z_sadcp = -unique(data(:,4));
NDEPTH = length(z_sadcp);

% load zonal/meridional velocities
u_sadcp = NaN(NDEPTH,NT);
v_sadcp = NaN(NDEPTH,NT);
for it=1:NT
  ind = find(data(:,1)==tim_sadcp(it));
  kk = [];
  for i=1:length(ind)
      kk(i) = find(z_sadcp==-data(ind(i),4));
  end
  u_sadcp(kk,it) = data(ind,5);
  v_sadcp(kk,it) = data(ind,6);
end

% convert velocity data in m/s
u_sadcp = 0.01*u_sadcp;
v_sadcp = 0.01*v_sadcp;

% remove data below 1200 meters
[zmax,imax] = max(abs(z_sadcp));
kk = find(z_sadcp<-1200);
z_sadcp(kk) = [];
u_sadcp(kk,:) = [];
v_sadcp(kk,:) = [];

% store the data
save6(files.sadcp,'tim_sadcp','lon_sadcp','lat_sadcp','u_sadcp','v_sadcp','z_sadcp')
