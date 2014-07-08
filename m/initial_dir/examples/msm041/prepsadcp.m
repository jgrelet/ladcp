function [] = prepsadcp(stn,values)
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
% LADCP processing, uncomment the next line

%disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepsadcp.m !')
%return

% first copy the SADCP files to the raw SADCP data directory
% data/raw_sadcp
% In our example
% this is data mounted via SMBMOUNT from //peale/adcp_home
%!cp /nbp/peale/www/data/* data/raw_sadcp

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
lon_sadcp = [];
lat_sadcp = [];
u_sadcp = [];
v_sadcp = [];
z_sadcp = [];
tim_sadcp = [];
for n=7:17
    load(['data/raw_sadcp/msm04_1100_000000_',int2str(n),'_hc'])
    lon_sadcp = [lon_sadcp,b.nav.txy1(2,:)];
    lat_sadcp = [lat_sadcp,b.nav.txy1(3,:)];
    tim_sadcp = [tim_sadcp,b.nav.txy1(1,:)];
    u_sadcp = [u_sadcp,squeeze(b.vel(:,1,:))];
    v_sadcp = [v_sadcp,squeeze(b.vel(:,2,:))];
end
z_sadcp = b.depth;
tim_sadcp = tim_sadcp + julian([2006 1 1 0 0 0]);


% restrict the data to the time of the cast
good = find( tim_sadcp>values.start_cut & tim_sadcp<values.end_cut);
tim_sadcp = tim_sadcp(good);
lat_sadcp = lat_sadcp(good);
lon_sadcp = lon_sadcp(good);
u_sadcp = u_sadcp(:,good);
v_sadcp = v_sadcp(:,good);
z_sadcp = z_sadcp;


% store the data
if str2num(version('-release'))>=14
  eval(['save data/sadcp/sadcp',int2str0(stn,3),...
	' tim_sadcp lon_sadcp lat_sadcp u_sadcp v_sadcp z_sadcp -v6'])
else
  eval(['save data/sadcp/sadcp',int2str0(stn,3),...
	' tim_sadcp lon_sadcp lat_sadcp u_sadcp v_sadcp z_sadcp'])
end
