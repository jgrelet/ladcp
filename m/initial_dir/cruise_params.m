% Parameters setting files are called in this order
%
% default_params.m
% cruise_params.m   <--- you are here
% cast_params.m
%
% this is the location to enter special settings which apply
% to a whole cruise or to your special LADCP system setup


% remove the following three lines after modifying the parameters
disp('edit  cruise_id/cruise_params.m')
pause
return


%
% set the  Cruise id
%
% this will appear on the top of all plots
% and in all file names
%
p.cruise_id	= 'CARI3';
p.name  = ['car3',int2str0(stn,3)];


%
% some software does only record the day of the year
% to be able to process such data properly enter the
% year which will be used for data without year information
%
% if you are measuring over newyear, you will need to introduce an
% if-statement here
%
p.correct_year = 2006;


%
% If you want you can give the serial numbers of up and down instrument
%
% this is just used in one plot
%
p.down_sn = NaN;
p.up_sn = NaN;


%
% Output resolution and superensemble averaging depth
%
% 20 is good for standard full ocean depth
% smaller (10 or even 5) can be used for special shallow casts
%
% default is down-looker bin-length
% 
%ps.dz	= 20;			% output depth resolution
%p.avdz	= ps.dz;		% pre-average data


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


% things to consider setting (mostly for experts)
% see default_params.m for descriptions
%
% p.btrk_mode (and other related ones)
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
% ps.sadcpfac	
% p.dragfac
% p.urange/zrange
% p.sadcp_dtok
