function [strlat,strlon] = str2pos(pos,pos2)
% function [strlat,strlon] = str2pos(pos)
%
% converts a position string into a decimal position
%
% Input: - pos                  : decimal position [lat,lon] in degrees
% 
% Output: - strlat              : position string of latitude ##o ##.###' N
%         - strlon              : position string of longitude ###o ##.###' E
%
% version 1.1.0         last change 17.09.2007

% Gerd Krahmann, IfM Kiel, Jun 1995
% added backwards compatibility of input args G.Krahmann    1.0.0-->1.0.1
% added degree signs G.K.                                   1.0.1-->1.0.2
% degree signs for different versions                       1.0.2-->1.1.0

if nargin==2
  pos = [pos,pos2];
end
latd = fix(pos(1));
latm = (pos(1)-latd)*60;
lond = fix(pos(2));
lonm = (pos(2)-lond)*60;
if latd<0 | latm<0
  ns = 'S';
  nsv = -1;
else
  ns = 'N';
  nsv = 1;
end
if lond<0 | lonm<0
  ew = 'W';
  ewv = -1;
else
  ew=  'E';
  ewv = 1;
end

% form latitude string
v = ver;
if str2num(v(1).Version(1))>5
  strlat = sprintf(['%2d^o %7.4f''',ns],abs([latd,latm]));
  strlon = sprintf(['%3d^o %7.4f''',ew],abs([lond,lonm]));
else
  strlat = sprintf(['%2d',char(176),' %7.4f''',ns],abs([latd,latm]));
  strlon = sprintf(['%3d',char(176),' %7.4f''',ew],abs([lond,lonm]));
end
