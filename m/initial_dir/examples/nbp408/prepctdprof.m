function [file] = prepctdprof(stn)
% function [file] = prepctdprof(stn)
%
% prepare CTD profile for LADCP
%
% this file is cruise specific for NBP 04-08
% 
% to create a file for your own cruise, copy this file and
% modify where necessary

% copy CTD profile to raw CTD data
% this data is a profile in 1dbar steps
eval(['!cp /nbp/science/NBP0408/ladcp/processed/p408',int2str0(stn,3),...
  	'ladcp1dbar.asc data/raw_ctd/p408',int2str0(stn,3),'.prf'])

% load this data and convert to standard format
% skip the header and extract only the PTS columns
[hdr,data] = hdrload(['data/raw_ctd/p408',int2str0(stn,3),'.prf']);
ctdprof = data(:,[4,2,8]);

% catch spikes
good = find(ctdprof(:,3)>1);
ctdprof = ctdprof(good,:);

% store data at the standard location
eval(['save data/ctd/ctdprof',int2str0(stn,3),' ctdprof'])

% save filename
file = ['data/ctd/ctdprof',int2str0(stn,3)];
