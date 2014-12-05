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
% LADCP processing, uncomment the next two lines. Otherwise edit the following.

disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepsadcp.m !')
pause
return

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
fnames = dir('data/raw_sadcp/*.mat');
uv = [];
xyt = [];
for n=1:length(fnames)
    load(['data/raw_sadcp/',fnames(n).name])
    if n==1
        uv = b.vel;
        xyt = b.nav.txy1;
    else
        uv = cat(3,uv,b.vel);
        xyt = cat(2,xyt,b.nav.txy1);
    end
end
lon_sadcp = xyt(2,:);
lat_sadcp = xyt(3,:);
tim_sadcp = xyt(1,:) + julian([2007 1 1 0 0 0]);
u_sadcp = squeeze(uv(:,1,:));
v_sadcp = squeeze(uv(:,2,:));
z_sadcp = c.depth;
[tim_sadcp,ind] = sort(tim_sadcp);
lon_sadcp = lon_sadcp(ind);
lat_sadcp = lat_sadcp(ind);
u_sadcp = u_sadcp(:,ind);
v_sadcp = v_sadcp(:,ind);

% restrict the data to the time of the cast
%good = find( tim_sadcp>values.start_time-0.1 & tim_sadcp<values.end_time+0.1);
good = find( tim_sadcp>values.start_cut-0.1 & tim_sadcp<values.end_cut+0.1);
tim_sadcp = tim_sadcp(good);
lat_sadcp = lat_sadcp(good);
lon_sadcp = lon_sadcp(good);
u_sadcp = u_sadcp(:,good);
v_sadcp = v_sadcp(:,good);
z_sadcp = z_sadcp;

% delete deepest two bins
z_sadcp = z_sadcp(1:end-2);
u_sadcp = u_sadcp(1:end-2,:);
v_sadcp = v_sadcp(1:end-2,:);


% store the data
save6(['data/sadcp/sadcp',int2str0(stn,3)],...
	'tim_sadcp','lon_sadcp','lat_sadcp','u_sadcp','v_sadcp','z_sadcp')
