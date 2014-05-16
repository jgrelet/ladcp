
function vel = sw_gvel(ga,r,lat)

% SW_CYCLOGEOVEL    Cyclogeostrophic velocity
%===================================================================
%
% USAGE:  vel = sw_cyclogeovel(ga,r,lat)
%
% DESCRIPTION:
%    Calculates cyclogeostrophic velocity given the geopotential anomaly,
%    radial position of each station and mean latitude of the section.
%
% INPUT:
%    ga   = geopotential anomoly relative to the sea surface.
%           dim(mxnstations)
%    r    = radial position of each station in m 
%    lat  = mean latitude
%
% OUTPUT:
%    vel  = cyclogeostrophic velocity RELATIVE to the sea surface.
%           dim(m,nstations-1)
%
% AUTHOR:    Gerd Krahmann 2008/04/06  (gkrahmann@ifm-geomar.de)
%            Phil Morgan   1992/03/26  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: S. Pond & G.Pickard  2nd Edition 1986
%            Introductory Dynamical Oceanogrpahy
%            Pergamon Press Sydney.  ISBN 0-08-028728-X
%            Equation 8.9A p73  Pond & Pickard
%
%==================================================================

% CALLER:   general purpose
% CALLEE:   sw_dist.m
%


DEG2RAD = pi/180;
RAD2DEG = 180/pi;
OMEGA   = 7.292e-5;  % Angular velocity of Earth  [radians/sec]

[m,n] = size(ga);
dr    = diff(r);
dr    = ones(m,1)*dr;
r     = ( r(1:end-1)+r(2:end) )/2;
r     = ones(m,1)*r;
f     = 2*OMEGA*sin( lat*DEG2RAD )*ones(m,n-1);
vel   = zeros(m,n-1);

for ii=1:5
disp(vel(100,3))
  lf    = (f+2*vel./r).*dr;
  vel   = -( ga(:,2:n)-ga(:,1:n-1) ) ./ lf;
end

return
%--------------------------------------------------------------------


