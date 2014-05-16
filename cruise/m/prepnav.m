function prepnav(stn,values)
% function prepnav(stn,values)
%
% prepare navigational data for LADCP
% we an array 'data' containing the 2 columns
% latitude in decimal degrees    longitude in decimal degrees
% and the vector 'timnav' containing the time of the navigational
% data in Julian days
%
% THIS FILE IS CRUISE SPECIFIC
%
% to create a file for your own cruise, modify this file
%
% The navigational data should be at a resolution of 1 per second.
% Lower resolution will lead to worse processing results.

% G.Krahmann, IFM-GEOMAR, Aug 2005

% if you do no have navigational data to be used in the
% LADCP processing, uncomment the next two lines. Otherwise edit the following.


% first copy navigational data to the raw NAV data directory
% data/raw_nav

% the navigation is read from the CTD file, if it exists
fname = ['m:/PANDORA/data-processing/CTD/data-ladcp/pn',int2str0(stn,5),'_ladcp.cnv'];
if exist(fname,'file')

   copyfile(fname,['data/raw_nav/',int2str0(stn,5),'.cnv']);
   [hdr,data] = read_sbe_cnv(['data/raw_nav/',int2str0(stn,5),'.cnv']);
   % time (in julian days)
   timnav = julian(datevec(data.datenum));
   % latitude/longitude array
   data = [data.latitude,data.longitude];

end

% remove NaN values
ind = find( isnan(timnav) | isnan(data(:,1)) | isnan(data(:,2)) );
timnav(ind) = [];
data(ind,:) = [];

% store data in the standard location
save6(['data/nav/nav',int2str0(stn,5)],'timnav','data')
