function [file] = prepladcp(stn)
% function [file] = prepladcp(stn)
%
% prepare LADCP data
%
% this file is cruise specific for NBP 04-08
% nothing to do here, as everything is already handled by the downloading routines
% 
% to create a file for your own cruise, copy this file and
% modify where necessary

if ~exist(['data/raw_ladcp/',int2str0(stn,3)])
  eval(['!mkdir data/raw_ladcp/',int2str0(stn,3)])
  eval(['!cp -a /nbp/science/NBP0408/ladcp/raw_backup/',int2str0(stn,3),...
	'DN/* data/raw_ladcp/',int2str0(stn,3)])
  eval(['!cp -a /nbp/science/NBP0408/ladcp/raw_backup/',int2str0(stn,3),...
	'UP/* data/raw_ladcp/',int2str0(stn,3)])
end
% set file name
file = ['data/raw_ladcp/',int2str0(stn,3)];
