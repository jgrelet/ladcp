function [ur,vr]=uvrot(u,v,rot);
% function [ur,vr]=uvrot(u,v,rot);
%
% rotate velocities clockwise
%
% input  :	u		- zonal/x velocity matrix
% 		v		- meridional/y velocity matrix
%		rot		- matrix of angles in degrees
%				  has to be same size as u/v or
%				  a scalar
%
% output :	ur		- rotated zonal/x velocity matrix
%		vr		- rotated meridional/y velocity matrix
%
% version 0.1	last change ?

% M.Visbeck
rot=-rot*pi/180;
cr=cos(rot);
sr=sin(rot);
ur=u.*cr-v.*sr;
vr=u.*sr+v.*cr;
