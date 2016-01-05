function misc_composefilenames(params, f)
% function [f] = misc_composefilenames(params,stn);
%
% compose the output filenames, this can't be done earlier
%
% input  :      params          - parameter structure
%               stn             - station number
%
% output :      f               - modified filename structure
%
% version 0.2  last change 08.11.2012

% GK, IFM-GEOMAR, Sep 2010

% moved stuff from default_params.m to here    GK, 08.11.2012  0.1-->0.2

% directory names
f.logs_dir        = 'logs';
f.plots_dir       = 'plots';
f.prof_dir        = 'profiles';
f.raw_dir         = 'data/raw_ladcp';
f.ctd_ts_dir      = 'data/ctdtime';
f.ctd_prof_dir    = 'data/ctdprof';
f.nav_dir         = 'data/nav';
f.sadcp_dir       = 'data/sadcp';

% file names
f.ladcpdo = strcat(f.raw_dir, filesep, params.ladcp_station_name, filesep, params.ladcp_station_name,'DN000.000');
f.ladcpup = strcat(f.raw_dir, filesep, params.ladcp_station_name, filesep, params.ladcp_station_name, 'UP000.000');

f.nav =     ['data/nav/nav',params.ladcp_station_name,'.mat'];
f.ctdprof = ['data/ctdprof/ctdprof',params.ladcp_station_name,'.mat'];
f.ctdtime = ['data/ctdtime/ctdtime',params.ladcp_station_name,'.mat'];
f.sadcp =   ['data/sadcp/sadcp',params.ladcp_station_name,'.mat'];

% file name for results (extensions will be added by software)
%  *.bot            bottom referenced ASCII data
%  *.lad            profile ASCII data
%  *.mat            MATLAB  format >> dr p ps f
%  *.cdf            NETCDF  (binary) LADCP data format 
%  *.log            ASCII log file of processing
%  *.txt            ASCII short log file
%  *.ps             post-script figure of result 

f.res =   strcat(f.prof_dir, filesep, params.name);
f.prof =  strcat(f.prof_dir, filesep, params.name);
f.plots = strcat(f.plots_dir,filesep, params.name);
f.log =   strcat(f.logs_dir, filesep, params.name);

if length(f.log) > 1                    % open log file
  logFile = strcat(f.log,'.log');
  if exist(logFile, 'file') == 2
    delete(logFile)
  end
  diary(logFile)
  diary on
end

