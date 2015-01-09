% Parameters setting files are called in this order
%
% default_params.m
% cruise_params.m   <--- you are here
% cast_params.m
%
% this is the location to enter special settings which apply
% to a whole cruise or to your special LADCP system setup


%
% set the  Cruise id
%
% this will appear on the top of all plots
% and in all file names
%
global pathFile;

if ispc
  thePath = 'g:\campagnes\';
elseif isunix
  thePath = '/M/';
end
p.cruise_id	= 'PANDORA';
p.name  = sprintf('pn%s', p.ladcp_station_name);
pathFile = strcat(thePath, p.cruise_id);


%
% some software does only record the day of the year
% to be able to process such data properly enter the
% year which will be used for data without year information
%
% if you are measuring over newyear, you will need to introduce an
% if-statement here
%
p.correct_year = 2012;

%
% If you want you can give the serial numbers of up and down instrument
%
% this is just used in one plot
p.down_sn = 133;
p.up_sn = 14131;


% distance between up and downlooker
% this will be added in rdiload to the distance of the
% first bin from the uplooking instrument
p.dist_up_down = 1.9;
%
% Output resolution and superensemble averaging depth
%
% 20 is good for standard full ocean depth
% smaller (10 or even 5) can be used for special shallow casts
%
% default is down-looker bin-length
%
%ps.dz	= 20;			% output depth resolution
p.avens	= 10;		% pre-average data


%
% Standard thresholds, beyond which data will be discarded
%
% elim : ADCP internal error velocity limit   0.5 is reasonable and default
% vlim : ADCP horizontal velocity limit       2.5 is reasonable and default
% wlim : ADCP vertical velocity bin limit     0.2 is reasonable and default
%
% (wlim is the deviation from the median of all bins in each ensemble)
%
%p.elim= 0.5;
%p.vlim= 2.5;
%p.wlim= 0.2;

% restrict time range to profile and disregard data close to surface
% p.cut = 0 dont restrict
% p.cut > 0 restrict time to adcp depth below a depth of p.cut
p.cut = 10;

% Write matlab file
%p.savemat = 0;

% Error using netcdflib
% The NetCDF library encountered an error during execution of 'defVar' function
% - 'String match to name in use (NC_ENAMEINUSE)'.
%
% Error in netcdf.defVar (line 38)
% varid = netcdflib('defVar', ncid, varname, xtype, dimids);
%
% Write netcdf file
% Error in ladcp2cdf (line 104)
%     varID = netcdf.defVar(nc,fnames{n},'float',zbot_dimID);
p.savecdf = 1;

%p.saveplot = [];

% In case the two instruments are running not synchronous, one
% is resampled onto the other. This is done by simply taking
% one instrument as the reference (default the downlooker) and
% each of its ensembles pick the closest in time of the other
% instrument. Depending on the ping rates and which instrument
% is pinging faster, this will result in whole ensembles being
% dropped or used multiple times.
%
% params.up2down==0 will not resample, unless different ping rates are detected
% params.up2down==1 will resample the uplooker onto the downlooker
% params.up2down==2 will resample the downlooker onto the uplooker
p.up2down = 1;

p.outlier = [];

% things to consider setting (mostly for experts)
% see default_params.m for descriptions

% BOTTOM TRACK
% 	The are several options to get bottom track data
% 
% mode = 1 :   use only RDI bottom track
%        2 :   use only own bottom track
%        3 :   use RDI, if existent, own else (default)
%        0 :   use not bottom track at all
%
% btrk_mode is the one you set to be used, btrk_used is the 
% the routine has used. They do not need to agree since the
% software can override your command, e.g. in case some
% necessary data is missing or all BTRK are bad.
p.btrk_mode = 0; %(and other related ones)

% p.rotup2down
% p.offsetup2down
% p.outlier
% p.pglim
% p.tiltmax
% p.tiltweight
% p.trusted_i
% p.barofac
% p.botfac
% p.smoofac
% p.smallfac

% Weight for SADCP data
% ps.sadcpfac=1 about equal weight for SADCP profile
ps.sadcpfac = 0;

% use SADCP data (1) or not (0)
% default is to use the data. This is just a simple switch
% to simplify testing of results
p.use_sadcp = 0;

% p.dragfac
% p.urange/zrange
% p.sadcp_dtok

% what is the format of the output plots
% can be multiple ones, separated by comma
% e.g.  params.print_formats = 'ps,jpg,png';
p.print_formats = 'png';
