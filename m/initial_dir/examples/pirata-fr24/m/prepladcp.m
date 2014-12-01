function [file] = prepladcp(params, files)
% function [file] = prepladcp(stn)
%
% prepare LADCP data for LADCP processing
%
% we need the raw LADCP data to be in the correct place and
% have the correct names.
%
% THIS FILE IS CRUISE SPECIFIC
%
% to create a file for your own cruise, modify this file
%
% you will just need to copy and possibly rename the files
% In case of old BB and NB systems you might need to append
% the raw data files.
%
% the convention for filenames is
%
% xxxDN000.000  and  xxxUP000.000  	with xxx the 3-digit station number
%
% they need to be copied into one directory per station
% data/raw_ladcp/xxx		with xxx the 3-digit station number

% G.Krahmann, IFM-GEOMAR, Aug 2005
global pathFile;

% path pour windows
if ispc
  if ~exist(['data\raw_ladcp\', params.ladcp_station_name])
    eval(['!mkdir data\raw_ladcp\',params.ladcp_station_name])
  end
else
  if ~exist(['data/raw_ladcp/',params.ladcp_station_name])
    eval(['!mkdir data/raw_ladcp/',params.ladcp_station_name])
  end
end

% downward-looking L-ADCP file
fname = strcat( pathFile, '/data-processing/LADCP/data/FR24M',...
  params.ladcp_station_name, '.000');
if exist(fname,'file')
  % LADCP was inverted before station 28
  if params.ladcp_station <= 27
    copyfile(fname, files.ladcpup);
  else
    copyfile(fname, files.ladcpdo);
 end
end

% upward-looking L-ADCP file
fname = strcat(pathFile, '/data-processing/LADCP/data/FR24S',...
  params.ladcp_station_name, '.000');
if exist(fname,'file')
    % LADCP was inverted before station 28
    if params.ladcp_station <= 27
      copyfile(fname, files.ladcpdo);
    else
      copyfile(fname, files.ladcpup);
    end
end

% set file name
file = ['data/raw_ladcp/',params.ladcp_station_name];

