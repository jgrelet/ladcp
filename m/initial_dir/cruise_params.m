% Parameters setting files are called in this order
%
% default_params.m
% cruise_params.m   <--- you are here
% cast_params.m
%
% this is the location to enter special settings which apply
% to a whole cruise or to your special LADCP system setup
%
% version 0.2  last change 20.05.2011

% changed handling of cruise_id                GK, 20.05.2011  0.1-->0.2


% remove the following three lines after modifying the parameters
disp('edit  cruise_id/cruise_params.m')
pause
return


%
% set the  Cruise id
%
% this will appear on the top of all plots
%
p.cruise_id	= 'CRUISE';


%
% file names will be composed by the directory name followed by '_NNN'
% where NNN is the station number
%
% here you could override this default behaviour
% you should however leave the station number part as in the example
%
%p.name  = ['cruise_id',int2str0(stn,3)];


%
% some software does only record the day of the year
% to be able to process such data properly enter the
% year which will be used for data without year information
%
% if you are measuring over newyear, you will need to introduce an
% if-statement here
%
p.correct_year = 2011;


%
% If you want you can give the serial numbers of up and down instrument
%
% this is just used in one plot
%
% these should starting at version 10.12 be extracted from the raw data
% files. Use these positions to override the extraction.
%
%p.down_sn = NaN;
%p.up_sn = NaN;


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
% p.edit_mask_dn_bins
% p.edit_mask_up_bins
