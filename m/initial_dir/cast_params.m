% Parameters setting files are called in this order
%
% default_params.m 
% cruise_params.m
% cast_params.m     <--- you are here
%
% this is the location to enter special settings which apply
% only to a single or a few casts. E.g. positions and times
%
% parameters to consider setting here are:
%
% p.poss and p.pose
% p.zpar
% p.time_start and p.time_end
%
% descriptions of the formats can be found in m/ladcp/default_params.m


% remove the following three lines after modifying the parameters
disp('edit  cruise_id/cast_params.m')
pause
return




if stn==1
%  p.time_start = [2006 05 16 11 24 00];
%  p.time_end = [2006 05 16 12 45 00];
%  p.zpar = [0, 1091,0];
%  p.poss = [10 32.436 -64 -33.302];
%  p.pose = [10 33.464 -64 -32.023];
elseif stn==2
elseif stn==3
end
