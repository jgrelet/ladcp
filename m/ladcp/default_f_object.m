classdef default_f_object < handle
  % class default_f_object contains file names
  %
  % Gerd Krahmann, Kiel July 2006
  % gkrahmann@ifm-geomar.de
  % Frederic Marin IRD - LEGOS - Noumea - june 2012
  % Frederic.Marin@ird.fr
  % Jacques Grelet R/V Le Suroit May 2014 - IRD US191
  % jacques.grelet@ird.fr
  %
  % $Id$
  %
  % file name for results (extensions will be added by software)
  %  *.bot            bottom referenced ASCII data
  %  *.lad            profile ASCII data
  %  *.mat            MATLAB  format >> dr p ps f
  %  *.cdf            NETCDF  (binary) LADCP data format
  %  *.log            ASCII log file of processing
  %  *.txt            ASCII short log file
  %  *.ps             post-script figure of result
  
  properties (Access = public)
    ladcpdo
    ladcpup
    res
    prof
    plots
    log
    nav
    ctdprof
    ctdtime
    sadcp
  end
  
  properties (Access = public, Hidden)
    logs_dir        = 'logs'
    plots_dir       = 'plots'
    prof_dir        = 'profiles'
    raw_dir         = 'data/raw_ladcp'
    ctd_ts_dir      = 'data/ctdtime'
    ctd_prof_dir    = 'data/ctdprof'
    nav_dir         = 'data/nav'
    sadcp_dir       = 'data/sadcp'
  end
  
  %% public methods
  methods
    % constructor
    % varargin must in the '%0xd' format

    function self = default_f_object(stn, ndigits)
      % pre initialization
      
      % if stn is numeric with number of valid digit
      if nargin == 2 && isnumeric(stn)
        stn_str = int2str0(stn, ndigits);
      end
      if nargin == 1 && isnumeric(stn)
        stn_str = int2str0(stn, 3);
      end
      
      % if stn is a string
      if ischar(stn)
        stn_str = stn;
      end

      self.ladcpdo = strcat(self.raw_dir, filesep, stn_str, filesep,...
        stn_str,'DN000.000');
      self.ladcpup = strcat(self.raw_dir, filesep, stn_str, filesep,...
        stn_str,'UP000.000');
      self.res = strcat(self.prof_dir, filesep, stn_str);
      self.prof = strcat(self.prof_dir, filesep, stn_str);
      self.plots = strcat(self.plots_dir, filesep, stn_str);
      self.log = strcat(self.logs_dir, filesep, stn_str);
      self.nav = strcat('data/nav/nav',stn_str,'.mat');
      self.ctdprof = strcat('data/ctdprof/ctdprof',stn_str,'.mat');
      self.ctdtime = strcat('data/ctdtime/ctdtime',stn_str,'.mat');
      self.sadcp = strcat('data/sadcp/sadcp',stn_str,'.mat');
    end
    

    
  end % end of public methods
  
end % end of class

