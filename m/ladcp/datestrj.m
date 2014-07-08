function a=datestrj(b)
%
% print julian date string
%
a=datestr(b-julian([0 1 0 0 0 0]));
