classdef default_p_object < dynamicprops
  % set default values for parameter in LADCP processing
  %
  %  structure p.??? contains parameter relevant to reading and preparing
  %             the data
  %
  % Gerd Krahmann, Kiel July 2006
  % gkrahmann@ifm-geomar.de
  % Frederic Marin IRD - LEGOS - Noumea - june 2012
  % Frederic.Marin@ird.fr
  % Jacques Grelet IRD US191 - Plouzane - May/Dec 2014
  % jacques.grelet@ird.fr
  %
  % $Id$
  
  properties (Access = public, Hidden, Dependent)

  end
  
  properties (Access = public)
    
    % software version display at startup
    software = 'IFM-GEOMAR LADCP software: Version 10.16.2'
    software_date = '07 February 2009 - 01 December 2014'
    whoami = 'unknown'
    
    % store station number
    ladcp_station
    
    % store station number in char
    ladcp_station_name
    
    % for cruise
    name = ' '
    cruise_id = []
    correct_year
    
    % preset start and end time vectors
    time_start = []
    time_end = []
    
    % restrict time range to profile and disregard data close to surface
    % p.cut = 0 dont restrict
    % p.cut > 0 restrict time to adcp depth below a depth of p.cut
    cut = 10
    
    % manually set POSITION of the start and end point
    % [degree lat, minute lat, degree lon, minute lon]
    %  i.e. [-59 -30.5697 -44 -22.4986]
    pose = [NaN NaN NaN NaN]
    poss = [NaN NaN NaN NaN]
    position_fixed = 0
    lat_for_calc = 0
    lon_for_calc = 0
    
    % navigation error in m
    %
    % This one is later used to determine the weight of the ship
    % movement constraint. It should be set to something like the
    % uncertainty of the position of the CTD when it is coming on
    % deck. I.e. the GPS error plus something accounting for the
    % position difference between GPS antenna and CTD (remember
    % the ship can rotate !).
    % 30 m is a reasonable number.
    nav_error = 30
    
    % SUPER ENSEMBLES
    % 	are calculated in prepinv.m to reduce the number of raw profiles
    % 	The ides is to obtain one average profile for each vertical dz=const
    % 	that the CTD traveled through. As a result a constant number of super
    % 	ensembles are obtained for the up and down cast.
    %   but a fixed number of ensembles can also be averaged
    
    
    %
    % p.avdz sets the depth interval between adjacent super-ensembles
    % default one bin length
    %
    % since the bin length is profile dependent we can not set it here
    %
    % if avdz is negative, -avdz bin lengths will be used for avdz
    % if avdz is positive, avdz will be used as is
    % the calculation takes place in  calc_ens_av.m
    %
    % p.avens will take precedence over p.avdz
    avdz = -1
    
    % p.avens overrides p.avdz and sets a fixed number of ensembles to average
    % default NAN means that it is not used
    avens = NaN
    
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
    btrk_mode = 3
    btrk_used = 0
    
    % p.btrk_ts is in dB to detect bottom above bin1 level (for own btm track)
    %
    % The following parameter (bottom-tracking target strength)
    % is quite iffy. Setting it to too small a value (e.g. 10, the old default)
    % makes instrument interference appear as a false bottom, which can
    % be a problem. For example, on CLIVAR P02 station 32 there was a
    % long stop 16m above the sea bed, which is too close for bottom
    % tracking.  True bottom detection is only possible during the approach
    % and the beginning of the upcast. These two short times are swamped
    % out by the interference-related false bottom detections. On the other
    % hand, when this value is set to too large a value, true seabed returns
    % are rejected. It would be fairly easy to set the correct value in
    % terms of instrument-returned target strength by plotting the target
    % strength with imagesc(d.wts). However, the values of this parameter
    % are not instrument-returned target strengths but db. The current value
    % is valid for the combo of downlooking BB150 / uplooking WH300 used on
    % the CLIVAR P02 cruise. It was derived by trial and error, involving
    % stations 2 and 32.
    btrk_ts = 30
    
    % p.btrk_below gives binoffset used below target strength maximum
    % to make bottom track velocity
    btrk_below = 0.5000
    
    % p.btrk_range gives minumum / maximum distance for bottom track
    % will be set to 2 / 10 binlength later
    btrk_range = [50 300]
    %
    % p.btrk_wstd gives maximum accepted wstd for super ensemble averages
    %
    % this one is usually calculated from the BTRK data but you can
    % override it
    btrk_wstd
    %btrk_wstd = 0.1;
    %
    % maximum allowed difference between reference layer W and W bottom track
    btrk_wlim = 0.0500
    
    % save bias and std of bottom track anomaly, used inside savearch
    btrk_u_bias
    btrk_u_std
    btrk_v_bias
    btrk_v_std
    btrk_w_bias
    btrk_w_std
    
    % force to recalculate bottom distance using target strength
    % this is turned off (0) by default and is being automatically
    % turned on (1), if the routine recognizes a large number of
    % bad (equal 0) RDI distances to the bottom
    bottomdist = 0
    
    % p.surfdist = 1 use surface reflections of up looking ADCP to get start
    % depth
    surfdist = 1
    
    % MAGNETIC deviation in degree
    magdev = 0
    
    %  time offset between upward-looking and downward-looking L-ADCPs
    %	>0 means down is in advance with respect to up
    %	<0 means up is in advance with respect to down
    %  this parameter may be necessary only when the 2 L-ADCP are not synchronized
    up2down_time_offset = 0
    
    % COMPASS manipulation
    % experts only
    % fix_compass:1 means hdg_offset gets added
    % fix_compass:2 means up looker gets down compass + hdg_offset
    % fix_compass:3 means down looker gets up compass + hdg_offset
    fix_compass = 0
    
    % In case the two instruments are running not synchronous, one
    % is resampled onto the other. This is done by simply taking
    % one instrument as the reference (default the downlooker) and
    % each of its ensembles pick the closest in time of the other
    % instrument. Depending on the ping rates and which instrument
    % is pinging faster, this will result in whole ensembles being
    % dropped or used multiple times.
    %
    % params.up2down==1 will resample the uplooker onto the downlooker
    up2down = 0
    
    % give compass offset in addition to declination (1) for down (2) for up
    %
    % p=setdefv(p,'hdg_offset',[0 0]);
    hdg_offset = [0 0]
    
    % COMPASS:
    % how to best adjust compass to best match
    % if 1 rotate up-looking and down-looking instrument to mean heading
    %    2 rotate up-looking and down-looking velocities to match up velocities
    %        (not really recommended)
    %    3 rotate up-looking velocities to down heading
    %        (use if suspect the up heading is bad
    %    4 rotate down-looking velocities to up heading
    %        (use if suspect the down heading is bad
    rotup2down = 1
    
    % Offset correction
    % if 1 remove velocity offset between up and down looking ADCP
    % this will correct errors due to tilt biases etc.
    offsetup2down = 1
    
    % DEPTH of the start, bottom and end of the profile
    % positive downwards
    zpar = [0 NaN 0]
    
    % maximum number of bins to be used
    % 0 : all will get used
    maxbinrange = 0
    
    % set ctdmaxlag.=100 to the maximum pings that the ADCP data can be shifted to
    % best match W calculated from CTD pressure time series (loadctd)
    % If you have good times set it to 10... if your time base is questionable
    % you can set it 100 or more
    ctdmaxlag = 150
    ctdmaxlagnp = 600
    
    % set forced_adcp_ctd_lag in case the built-in routine does not
    % properly recognize the lag
    % usually this is not set at all
    % forced_adcp_ctd_lag;
    
    % save individual target strength p.ts_save=[1 2 3 4]
    ts_save = 0
    
    % save individual correlation p.cm_save=[1 2 3 4]
    cm_save = 0
    
    % save individual percent good pings p.pg_save=[1 2 3 4]
    pg_save = 0
    
    % 1200kHz WH data is too high resolution to be merged with other data
    % it thus needs to be averaged before being used
    % this gives the number of values to be averaged
    nav_1200 = 4
    
    % a 1200kHz measures very close to the rosette
    % high error velocities result when measured in the eddy tail of the rosette
    % this parameter sets the threshold of values to be discarded
    error_limit_1200 = 20 % not used GK
    extra_blank_1200 = 5	% extra blank in meters
    
    %OUTLIER detection is called twice once to clean the raw data
    %	and a second time to clean the super ensembles
    %        [n1 n2 n3 ...] the length gives the number of scans and
    %	each value the maximum allowed departure from the mean in std
    %	applied for the u,v,w fields for each bin over blocks
    %   of p.outlier_n profiles
    %
    % 2: very strong  3: medium  4:only largest outliers
    outlier = [4 3]
    
    % default for p.outlier_n is number of profiles in 5 minutes
    % p=setdefv(p,'outlier_n',100);
    % minimum std for horizontal velocities of super ensemble
    % p=setdefv(p,'superens_std_min',0.01);
    
    %SPIKES
    % 	maximum value for abs(V-error) velocity
    elim = 0.5000
    % 	maximum value for horizontal velocity
    vlim = 2.5000
    % 	minimum value for %-good
    pglim = 0
    %	maximum value for W difference between the mean W and actual
    %        W(z) for each profile.
    wlim = 0.2000
    
    % TILT  flag data with large tilt or tilt differences as bad
    % [22  (max tilt allowed)
    %  4 (maximum tilt difference between pings allowed)]
    % WH systems have reported decent profiles with up to 35 deg tilt ...
    tiltmax = [22 4]
    
    % TILT  reduce weight for large tilts
    % calculated after an obscure formula in prepinv.m GK
    tilt_weight = 10
    
    % fix TIME of the ADCP in days
    % positive numbers shift the ADCP's time to a later time
    % (i.e. timoff is added to the julian time of the ADCP)
    % Please note that here is only one offset possible.
    % This is the offset of the MASTER system.
    % Should the LADCP slave system have a different offset, it
    % will be handled by the data shifting of the slave.
    timoff = 0
    
    % fix time of uplooking ADCP relative to downlooking ADCP
    % e.g. when the clocks were never properly set
    %
    % Usually that is taken care of by the automatic shifting.
    % But this will fail when the uplooking system apparently start
    % before the downlooking one, but in reality starts later.
    %
    % In that case you can shift the uplooking data to later
    % times by giving a positive number. Units are days.
    % There will be a paused control plot to compare the vertical
    % velicities.
    %
    timoff_uplooker = 0;
    
    % by up to how many ensembles shall the slave be shifted
    % against the master so that the vertical velocities match
    % usually 20 is enough. But if the times of master and
    % slave were not properly set a much larger number might
    % be necessary
    maxlag = 20
    
    
    % there is still an old slower up/down lag routine implemented
    % it is turnned off by default and the newer one is used
    %
    % it will automatically turn on, if the bestlag found has a correlation
    % of less than 0.9
    %
    bestlag_testing_on = 0;
    
    % apply tilt correction
    % tiltcor(1)=down-pitch bias
    % tiltcor(2)=down-rol bias
    % tiltcor(3)=up-pitch bias
    % tiltcor(4)=up-rol bias
    tiltcor = 0
    
    % Give bin number for the best W to compute depth of the ADCP
    %	default uses bin 2-3 but be careful when up/down instruments
    %	are used. The good bins are in the middle!
    trusted_i = [2 3 4 5]
    
    % SET ambiguity velocity used [m/s]
    % not used for anything, but stored in the archival data !? GK
    ambiguity = 2.5000
    
    % Give single ping accuracy;
    % another strange one, set to NaN, multiplies something and
    % that's it  GK
    single_ping_accuracy = NaN
    
    % clear LADCPs pressure sensor
    % sometimes has strange data and w integration might be preferred
    %
    % 1 clears the pressure records
    clear_ladcp_pressure = 0
    weight_ladcp_pressure = 0.1
    
    % Write matlab file
    savemat = 0
    
    % Write netcdf file
    savecdf = 1
    
    % Save Plots
    % Save figure numbers to ps file
    %    1 : Summary Plot
    %    2 : Engineering Data
    %    3 : Data Quality
    %    4 : Depth
    %    5 : Heading Corrections
    %    6 : Up/Down Differences
    %    7 : CTD Position
    %    8 : Shear
    %    9 : SADCP U, V
    %   10 : U, V Offsets, Tilt Error
    %   11 : Processing Warnings
    %   12 : Inversion Constraints
    %   13 : Bottom Track detail
    %   14 : Target Strength
    %   15 : Correlation
    %   16 : Weights
    
    % FM - IRD - June 2011
    % Add figure 15 to the list of the figures to be plotted
    saveplot = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]
    
    % The following parameter (slack time for including SADCP data outside
    % LADCP-cast time interval) is by default set to 5 minutes. On CLIVAR_P02
    % this is too much, since we ran an efficient operation. Sometimes, the
    % ship was getting underway less than 5 minutes after the end of the cast.
    % This led to SADCP data outliers, which increases the standard deviation,
    % which makes the inversion reject all SADCP data (low weight).
    sadcp_dtok = 0 % no time slack for SADCP near ends of stn
    
    % use SADCP data (1) or not (0)
    % default is to use the data. This is just a simple switch
    % to simplify testing of results
    use_sadcp = 1
    
    % set the ranges for the main velocity plot
    % [nan,nan,nan,0] will use standard variable axis ranges
    % with this setting you can set it to something fixed for easier
    % comparison
    % units are [cm/s cm/s m m] (depth is negative)
    plot_range = [NaN NaN NaN 0]
    
    % set the maximum distance for the bottom track velocity plot
    btrk_plot_range = 400
    
    % PARAMETERS FOR EDIT_DATA
    % Set list of bins to always remove from data.
    edit_mask_dn_bins = []
    edit_mask_up_bins = []
    
    % set list of bins to fully remove from data
    % this will be done already at the loading stage and should very rarely be
    % necessary
    edit_hardremove_mask_dn_bins = [];
    edit_hardremove_mask_up_bins = [];
    
    % Set to 1 to remove side-lobe contaminated data near seabed and
    % surface.
    edit_sidelobes = 1
    
    % Set to finite value to implement time-domain spike filter on the data;
    % this removes interference from other acoustic instruments but,
    % more importantly, can get rid of PPI when staggered pings
    % are used.
    %
    % Spike filtering is done using 2nd-difference
    % peak detection in time. This parameter gives the maximum
    % target-strength 2nd derivative that's allowed. Set to larger
    % values to weaken the filtering. (Check figure 14 to see if
    % filter is too strong or too weak.)
    %
    % has been normalized to handle different instruments
    % old values from pre-10 versions will not work !!!  GK
    edit_spike_filter_max_curv = 0.1000
    
    % Set to 1 to remove data contaminated by previous-ping interference.
    % NB: using the spike filter seems to work more robustly, as long
    %     as staggered pings are used.
    edit_PPI = 0
    
    % PPI layer thickness in meters; the value is taken directly from Eric
    % Firing's default (2*clip_margin = 180m).
    edit_PPI_layer_thickness = 180
    
    % max distance from seabed at which PPI should be removed. This is
    % an observed parameter and depends on the clarity of the water.
    % Check Figure 14 to see whether this should be changed.
    edit_PPI_max_hab = 1000
    
    % set this vector to enable skipping of ensembles; skipping vector
    % is wrapped around, i.e. [1 0] skips all odd ensembles, [0 1 0] skips
    % ensembles 2 5 8 11.... This filter is useful to process the casts
    % with only half the data to see whether the two halves agree, which
    % probably means that the cast can be trusted. Note that if staggered
    % ping setup is used to avoid PPI, the skipping vector should leave
    % adjacent ensembles intact, i.e. use something like [1 1 0 0] and
    % [0 0 1 1].
    edit_skip_ensembles = []
    
    % a detection alogrithm for asynchronous ping interference between
    % master and slave has been developed. This is by default off, as
    % we assume that the system is run synchronous
    detect_asynchronous = 0
    
    % bottom depth derived from data
    % this one is used later for editing data
    %
    % if  0    : will be calculated, if data is available
    %     nan  : will not be used
    zbottom = 0
    zbottomerror = 0
    
    % from experience it appears as if the most distant data containing
    % bin is very often not good. This can be seen in Plot 3 where the
    % fringes of the data in the leftmost subplot are colored.
    % With this parameter this data is always discarded.
    edit_mask_last_bin = 0
    
    % misc other
    % LADCP cast number
    ladcp_cast = 1
    warnp = []
    
    % serial numbers of instruments
    % this is just used in one plot
    % overwrite them in cruise_params.m
    down_sn = NaN
    up_sn = NaN
    
    % distance between up and downlooker
    % this will be added in rdiload to the distance of the
    % first bin from the uplooking instrument
    dist_up_down = 0
    
    % what is the format of the output plots
    % can be multiple ones, separated by comma
    % e.g.  params.print_formats = 'ps,jpg,png';
    print_formats = 'ps'
    
    % discard all values in bins higher than the first one that is below
    % the threshold
    % the two values are for down and up looking instruments
    % if there is only a single value, it will be applied to both
    minimum_correlation_threshold = [0,0]
    
    % multiply the weight of up and/or down looker by a factor
    down_up_weight_factors = [1,1]
    
    % properties used inside processing
    up_range
    dn_range
    bins_u
    bins_d
    all_trusted_i
    ladcpr_CTD_depth_std
    dt_profile
    up_dn_comp_off
    up_dn_pit_rol_comp_off
    up_dn_rol_off
    up_dn_pit_off
    hbot_0
    
    % need iside process_cast
    warnings = []
    
    % the properties are dynamically setting inside processing step
    % using updated setdefv methods
    %
    %     outlier_n
    %     btrk_range
    %
    %     ts_att_dn
    %     ts_att_up
    
  end
  
  properties (Access = public, Hidden)
    logs_dir        = 'logs'
    plots_dir       = 'plots'
    prof_dir        = 'profiles'
    raw_dir         = 'data/raw_ladcp'
    ctd_ts_dir      = 'data/ctdtime'
    ctd_prof_dir    = 'data/ctdprof'
    nav_dir         = 'data/nav'
    sadcp_dir       = 'data/sadcp'
  end
  
  properties (Access = public, Dependent = true)
  end
  
  
  %% public methods
  methods
    % constructor
    function self = default_p_object(stn, ndigits)
      
      % pre initialization
      
      % if stn is numeric with number of valid digit
      if nargin == 2 && isnumeric(stn)
        self.ladcp_station = stn;
        self.ladcp_station_name = int2str0(stn, ndigits);
      end
      if nargin == 1 && isnumeric(stn)
        self.ladcp_station = stn;
        self.ladcp_station_name = int2str0(stn, 3);
      end      
      
      % if stn is a string
      if ischar(stn)
        self.ladcp_station_name = stn;
        self.ladcp_station = str2double(stn);
      end
      
      % getting information from the host and user
      [~,self.whoami] = system('whoami');
      
      % show version software
      fprintf(1, '%s\n%s\n', self.software, self.software_date);
      
    end % end of constructor
    
    % diplay object if public properties are set to hidden
    %     function display(self)
    %       mco = metaclass(self);
    %       plist = mco.PropertyList;
    %       for i = 1 : length(plist)
    %         fprintf(1, '%s: ', plist(i).Name);
    %         disp(self.(plist(i).Name));
    %       end
    %     end
    
    % getters and setters
    % usefull for type ckecking
    % this is just an example
    function value = get.software(self)
      value = self.software;
    end
    
    function set.software(self, value)
      self.software = value;
    end
    
    % get station number and name (in proper format) since v10.16
    function value = get.ladcp_station_name(self)
      value = self.ladcp_station_name;
    end
    
    function value = get.ladcp_station(self)
      value = self.ladcp_station;
    end
    
  end % end of public methods
  
end % end of class

