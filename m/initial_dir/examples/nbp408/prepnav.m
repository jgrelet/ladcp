function prepnav(stn,values)
% function prepnav(stn,values)
%
% prepare navigational data
%
% this file is cruise specific for NBP 04-08
% (CTD is used here as it contains navigational data)
% 
% to create a file for your own cruise, copy this file and
% modify where necessary

% copy CTD time data to nav data
eval(['!cp /nbp/science/NBP0408/ladcp/ctd/p408',int2str0(stn,3),...
  	'ladcp.cnv data/raw_ctd/p408',int2str0(stn,3),'.tim'])

% load this data and convert to standard format
[hdr,data] = hdrload(['data/raw_ctd/p408',int2str0(stn,3),'.tim']);
timnav = data(:,7) + julian([2004,1,0,0,0,0]);
data = data(:,[5,6]);

% crop data
good = find(timnav>=values.start_time & timnav<=values.end_time);
timnav = timnav(good);
data = data(good,:);

% store data in the standard location
eval(['save data/nav/nav',int2str0(stn,3),' timnav data'])
