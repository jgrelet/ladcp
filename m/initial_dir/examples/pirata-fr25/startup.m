% LADCP processing 10.0

% change where is your code
if ispc
  addpath(genpath('M:\PIRATA-FR25\data-processing\LADCP\v10.16.2'));
elseif isunix
  addpath(genpath('/M/PIRATA-FR25\data-processing\LADCP\v10.16.2'));
end

%addpath ../m/ladcp;
%addpath ../m/sw;
%addpath ../m/netcdf;
%addpath ../m/netcdf/nctype;
%addpath ../m/netcdf/ncutility;
%addpath ../m/mexnc;
addpath m;

set(0,'defaultsurfaceedgecolor','none');
%set(0,'defaultfigurerenderer','zbuffer')
%set(0,'defaultfigurerenderer','opengl')
set(0,'defaultfigurerenderer','painter')
