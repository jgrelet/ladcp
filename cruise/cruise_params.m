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
if ispc
  thePath = 'M:\';
elseif isunix
  thePath = '/M/';
end
p.cruise_id	= 'PIRATA-FR24';
p.name  = ['fr24',stn_str];
global pathFile;
pathFile = strcat(thePath, p.cruise_id);


%
% some software does only record the day of the year
% to be able to process such data properly enter the
% year which will be used for data without year information
%
% if you are measuring over newyear, you will need to introduce an
% if-statement here
%
p.correct_year = 2014;

%
% If you want you can give the serial numbers of up and down instrument
%
% this is just used in one plot
p.up_sn   = 10307;
p.down_sn = 14131;

% distance between up and downlooker
% this will be added in rdiload to the distance of the
% first bin from the uplooking instrument
p.dist_up_down = 1.5;
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

p.cut = 10;

%p.savemat = 0;
% p.savecdf = 1;
% Error using netcdflib
% The NetCDF library encountered an error during execution of 'defVar' function
% - 'String match to name in use (NC_ENAMEINUSE)'.
%
% Error in netcdf.defVar (line 38)
% varid = netcdflib('defVar', ncid, varname, xtype, dimids);
%
% Error in ladcp2cdf (line 104)
%     varID = netcdf.defVar(nc,fnames{n},'float',zbot_dimID);
p.savecdf = 0;
%p.saveplot = [];

p.up2down = 1;

p.outlier = [];

% things to consider setting (mostly for experts)
% see default_params.m for descriptions
%
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
ps.sadcpfac = 0;
p.use_sadcp = 0;
% p.dragfac
% p.urange/zrange
% p.sadcp_dtok
