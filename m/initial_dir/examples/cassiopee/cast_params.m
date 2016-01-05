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

% Set list of bins to always remove from data.
% update for FR24, all profils
%p.edit_mask_up_bins = [1 3:4];

% Give bin number for the best W to compute depth of the ADCP
%	default uses bin 2-3 but be careful when up/down instruments
%	are used. The good bins are in the middle! 
p.trusted_i = 5:8;

% switch stn_str
%     case '001'
%         p.edit_mask_up_bins = [1 3:4];
%         p.trusted_i = 5:8;
%         p.down_sn = 12818;
%         p.up_sn = 12817;%
% end
