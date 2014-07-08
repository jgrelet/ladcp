function [f] = misc_composefilenames(f,params,stn);
% function [f] = misc_composefilenames(f,params,stn);
%
% compose the output filenames, this can't be done earlier
%
% input  :      f               - filename structure
%               params          - parameter structure
%               stn             - station number
%
% output :      f               - modified filename structure
%
% version 0.1  last change 01.09.2010

% GK, IFM-GEOMAR, Sep 2010

f.res = [f.prof_dir,'/',params.name];
f.prof = [f.prof_dir,'/',params.name];
f.plots = [f.plots_dir,'/',params.name];
f.log = [f.logs_dir,'/',params.name];

if length(f.log) > 1                    % open log file
  if exist([f.log,'.log'],'file')==2
    delete([f.log,'.log'])
  end
  diary([f.log,'.log'])
  diary on
end
