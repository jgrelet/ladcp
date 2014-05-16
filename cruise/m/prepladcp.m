function [file] = prepladcp(stn_str, pathFile)
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
% path pour windows
if ispc
  if ~exist(['data\raw_ladcp\',stn_str])
    eval(['!mkdir data\raw_ladcp\',stn_str])
  end
else
  if ~exist(['data/raw_ladcp/',stn_str])
    eval(['!mkdir data/raw_ladcp/',stn_str])
  end
end

% downward-looking L-ADCP file
fname = strcat(pathFile,'/data-processing/LADCP/data/FR24M',stn_str,'.000');
if exist(fname,'file')
  if str2double(stn_str) <= 27
    copyfile(fname,['data/raw_ladcp/',stn_str,'/',stn_str,'UP000.000']);
  else
    copyfile(fname,['data/raw_ladcp/',stn_str,'/',stn_str,'DN000.000']);
 end
end

% upward-looking L-ADCP file
fname = strcat(pathFile,'/data-processing/LADCP/data/FR24S',stn_str,'.000');
if exist(fname,'file')
    if str2double(stn_str) <= 27
      copyfile(fname,['data/raw_ladcp/',stn_str,'/',stn_str,'DN000.000']);
    else
      copyfile(fname,['data/raw_ladcp/',stn_str,'/',stn_str,'UP000.000']);
    end
end

% set file name
file = ['data/raw_ladcp/',stn_str];

