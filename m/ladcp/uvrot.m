function [ur,vr]=uvrot(u,v,rot);
% function [ur,vr]=uvrot(u,v,rot);
%
% rotate velocities clockwise
%
% input  :  u   - zonal/x velocity matrix
%           v   - meridional/y velocity matrix
%           rot	- scalar/vector/matrix of angles in degrees
%                 has to be same size as u/v or
%                 a scalar or has to have one agreeing dimension with u/v
%
% output :  ur  - rotated zonal/x velocity matrix
%           vr  - rotated meridional/y velocity matrix
%
% version 0.2	last change  29.05.2011

% M.Visbeck

% added dimension matching                    GK, 29.05.2011  0.1-->0.2

if size(rot,1)==1
  rot = ones(size(u,1),1)*rot;
end
if size(rot,2)==1
  rot = rot*ones(1,size(u,2),1);
end
rot = -rot*pi/180;
cr = cos(rot);
sr = sin(rot);
ur = u.*cr-v.*sr;
vr = u.*sr+v.*cr;
