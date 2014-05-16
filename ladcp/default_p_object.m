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
  % Jacques Grelet R/V Le Suroit May 2014 - IRD US191
  % jacques.grelet@ird.fr
  %
  % $Id$
  
  properties (Access = public)
    software = 'IFM-GEOMAR LADCP software: Version 10.8: 07 February 2009 '
    whoami = 'unknown'
    name = ' '
    
    % preset start and end time vectors
    time_start = []
    time_end = []
    
    % store station number
    ladcp_station = ''
    
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
    %  params.btrk_range = [50 300];
    %
    % p.btrk_wstd gives maximum accepted wstd for super ensemble averages
    %
    % this one is usually calculated from the BTRK data but you can
    % override it
    %
    % p.btrk_wstd = 0.1;
    %
    % maximum allowed difference between reference layer W and W bottom track
    btrk_wlim = 0.0500
    
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
    
    % by up to how many ensembles shall the slave be shifted
    % against the master so that the vertical velocities match
    % usually 20 is enough. But if the times of master and
    % slave were not properly set a much larger number might
    % be necessary
    maxlag = 20
    
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
    
    % FM - IRD - June 2011
    % Add figure 15 to the list of the figures to be plotted
    saveplot = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
    
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
    
    % PARAMETERS FOR EDIT_DATA
    % Set list of bins to always remove from data.
    edit_mask_dn_bins = []
    edit_mask_up_bins = []
    
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
    down_sn = 0
    up_sn = 0
    
    % distance between up and downlooker
    % this will be added in rdiload to the distance of the
    % first bin from the uplooking instrument
    dist_up_down = 0
    
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
    function self = default_p_object()
      self.whoami = whoami; %#ok<CPROP>
      disp(self.software);          % show version
    end
    
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
    % usefull for realtime type ckecking
    function value = get.software(self)
      value = self.software;
    end
    
    function value = get.whoami(self)
      value = self.whoami;
    end
    function value = get.name(self)
      value = self.name;
    end
    function value = get.time_start(self)
      value = self.time_start;
    end
    function value = get.time_end(self)
      value = self.time_end;
    end
    function value = get.ladcp_station(self)
      value = self.ladcp_station;
    end
    function value = get.cut(self)
      value = self.cut;
    end
    function value = get.pose(self)
      value = self.pose;
    end
    function value = get.poss(self)
      value = self.poss;
    end
    function value = get.position_fixed(self)
      value = self.position_fixed;
    end
    function value = get.lat_for_calc(self)
      value = self.lat_for_calc;
    end
    function value = get.lon_for_calc(self)
      value = self.lon_for_calc;
    end
    function value = get.nav_error(self)
      value = self.nav_error;
    end
    function value = get.avdz(self)
      value = self.avdz;
    end
    function value = get.avens(self)
      value = self.avens;
    end
    function value = get.btrk_mode(self)
      value = self.btrk_mode;
    end
    function value = get.btrk_used(self)
      value = self.btrk_used;
    end
    function value = get.btrk_ts(self)
      value = self.btrk_ts;
    end
    function value = get.btrk_below(self)
      value = self.btrk_below;
    end
    function value = get.btrk_wlim(self)
      value = self.btrk_wlim;
    end
    function value = get.bottomdist(self)
      value = self.bottomdist;
    end
    function value = get.surfdist(self)
      value = self.surfdist;
    end
    function value = get.magdev(self)
      value = self.magdev;
    end
    function value = get.up2down_time_offset(self)
      value = self.up2down_time_offset;
    end
    function value = get.fix_compass(self)
      value = self.fix_compass;
    end
    function value = get.up2down(self)
      value = self.up2down;
    end
    function value = get.rotup2down(self)
      value = self.rotup2down;
    end
    function value = get.offsetup2down(self)
      value = self.offsetup2down;
    end
    function value = get.zpar(self)
      value = self.zpar;
    end
    function value = get.maxbinrange(self)
      value = self.maxbinrange;
    end
    function value = get.ctdmaxlag(self)
      value = self.ctdmaxlag;
    end
    function value = get.ctdmaxlagnp(self)
      value = self.ctdmaxlagnp;
    end
    function value = get.ts_save(self)
      value = self.ts_save;
    end
    function value = get.cm_save(self)
      value = self.cm_save;
    end
    function value = get.pg_save(self)
      value = self.pg_save;
    end
    function value = get.nav_1200(self)
      value = self.nav_1200;
    end
    function value = get.error_limit_1200(self)
      value = self.error_limit_1200;
    end
    function value = get.extra_blank_1200(self)
      value = self.extra_blank_1200;
    end
    function value = get.outlier(self)
      value = self.outlier;
    end
    function value = get.elim(self)
      value = self.elim;
    end
    function value = get.vlim(self)
      value = self.vlim;
    end
    function value = get.pglim(self)
      value = self.pglim;
    end
    function value = get.wlim(self)
      value = self.wlim;
    end
    function value = get.tiltmax(self)
      value = self.tiltmax;
    end
    function value = get.tilt_weight(self)
      value = self.tilt_weight;
    end
    function value = get.timoff(self)
      value = self.timoff;
    end
    function value = get.maxlag(self)
      value = self.maxlag;
    end
    function value = get.tiltcor(self)
      value = self.tiltcor;
    end
    function value = get.trusted_i(self)
      value = self.trusted_i;
    end
    function value = get.ambiguity(self)
      value = self.ambiguity;
    end
    function value = get.single_ping_accuracy(self)
      value = self.single_ping_accuracy;
    end
    function value = get.clear_ladcp_pressure(self)
      value = self.clear_ladcp_pressure;
    end
    function value = get.weight_ladcp_pressure(self)
      value = self.weight_ladcp_pressure;
    end
    function value = get.savemat(self)
      value = self.savemat;
    end
    function value = get.savecdf(self)
      value = self.savecdf;
    end
    function value = get.saveplot(self)
      value = self.saveplot;
    end
    function value = get.sadcp_dtok(self)
      value = self.sadcp_dtok;
    end
    function value = get.use_sadcp(self)
      value = self.use_sadcp;
    end
    function value = get.plot_range(self)
      value = self.plot_range;
    end
    function value = get.edit_mask_dn_bins(self)
      value = self.edit_mask_dn_bins;
    end
    function value = get.edit_mask_up_bins(self)
      value = self.edit_mask_up_bins;
    end
    function value = get.edit_sidelobes(self)
      value = self.edit_sidelobes;
    end
    function value = get.edit_spike_filter_max_curv(self)
      value = self.edit_spike_filter_max_curv;
    end
    function value = get.edit_PPI(self)
      value = self.edit_PPI;
    end
    function value = get.edit_PPI_layer_thickness(self)
      value = self.edit_PPI_layer_thickness;
    end
    function value = get.edit_PPI_max_hab(self)
      value = self.edit_PPI_max_hab;
    end
    function value = get.edit_skip_ensembles(self)
      value = self.edit_skip_ensembles;
    end
    function value = get.detect_asynchronous(self)
      value = self.detect_asynchronous;
    end
    function value = get.zbottom(self)
      value = self.zbottom;
    end
    function value = get.zbottomerror(self)
      value = self.zbottomerror;
    end
    function value = get.edit_mask_last_bin(self)
      value = self.edit_mask_last_bin;
    end
    function value = get.ladcp_cast(self)
      value = self.ladcp_cast;
    end
    function value = get.warnp(self)
      value = self.warnp;
    end
    function value = get.down_sn(self)
      value = self.down_sn;
    end
    function value = get.up_sn(self)
      value = self.up_sn;
    end
    function value = get.dist_up_down(self)
      value = self.dist_up_down;
    end
    
    function set.whoami(self, value)
      self.whoami = value;
    end
    function set.name(self, value)
      self.name = value;
    end
    function set.time_start(self, value)
      self.time_start = value;
    end
    function set.time_end(self, value)
      self.time_end = value;
    end
    function set.ladcp_station(self, value)
      self.ladcp_station = value;
    end
    function set.cut(self, value)
      self.cut = value;
    end
    function set.pose(self, value)
      self.pose = value;
    end
    function set.poss(self, value)
      self.poss = value;
    end
    function set.position_fixed(self, value)
      self.position_fixed = value;
    end
    function set.lat_for_calc(self, value)
      self.lat_for_calc = value;
    end
    function set.lon_for_calc(self, value)
      self.lon_for_calc = value;
    end
    function set.nav_error(self, value)
      self.nav_error = value;
    end
    function set.avdz(self, value)
      self.avdz = value;
    end
    function set.avens(self, value)
      self.avens = value;
    end
    function set.btrk_mode(self, value)
      self.btrk_mode = value;
    end
    function set.btrk_used(self, value)
      self.btrk_used = value;
    end
    function set.btrk_ts(self, value)
      self.btrk_ts = value;
    end
    function set.btrk_below(self, value)
      self.btrk_below = value;
    end
    function set.btrk_wlim(self, value)
      self.btrk_wlim = value;
    end
    function set.bottomdist(self, value)
      self.bottomdist = value;
    end
    function set.surfdist(self, value)
      self.surfdist = value;
    end
    function set.magdev(self, value)
      self.magdev = value;
    end
    function set.up2down_time_offset(self, value)
      self.up2down_time_offset = value;
    end
    function set.fix_compass(self, value)
      self.fix_compass = value;
    end
    function set.up2down(self, value)
      self.up2down = value;
    end
    function set.rotup2down(self, value)
      self.rotup2down = value;
    end
    function set.offsetup2down(self, value)
      self.offsetup2down = value;
    end
    function set.zpar(self, value)
      self.zpar = value;
    end
    function set.maxbinrange(self, value)
      self.maxbinrange = value;
    end
    function set.ctdmaxlag(self, value)
      self.ctdmaxlag = value;
    end
    function set.ctdmaxlagnp(self, value)
      self.ctdmaxlagnp = value;
    end
    function set.ts_save(self, value)
      self.ts_save = value;
    end
    function set.cm_save(self, value)
      self.cm_save = value;
    end
    function set.pg_save(self, value)
      self.pg_save = value;
    end
    function set.nav_1200(self, value)
      self.nav_1200 = value;
    end
    function set.error_limit_1200(self, value)
      self.error_limit_1200 = value;
    end
    function set.extra_blank_1200(self, value)
      self.extra_blank_1200 = value;
    end
    function set.outlier(self, value)
      self.outlier = value;
    end
    function set.elim(self, value)
      self.elim = value;
    end
    function set.vlim(self, value)
      self.vlim = value;
    end
    function set.pglim(self, value)
      self.pglim = value;
    end
    function set.wlim(self, value)
      self.wlim = value;
    end
    function set.tiltmax(self, value)
      self.tiltmax = value;
    end
    function set.tilt_weight(self, value)
      self.tilt_weight = value;
    end
    function set.timoff(self, value)
      self.timoff = value;
    end
    function set.maxlag(self, value)
      self.maxlag = value;
    end
    function set.tiltcor(self, value)
      self.tiltcor = value;
    end
    function set.trusted_i(self, value)
      self.trusted_i = value;
    end
    function set.ambiguity(self, value)
      self.ambiguity = value;
    end
    function set.single_ping_accuracy(self, value)
      self.single_ping_accuracy = value;
    end
    function set.clear_ladcp_pressure(self, value)
      self.clear_ladcp_pressure = value;
    end
    function set.weight_ladcp_pressure(self, value)
      self.weight_ladcp_pressure = value;
    end
    function set.savemat(self, value)
      self.savemat = value;
    end
    function set.savecdf(self, value)
      self.savecdf = value;
    end
    function set.saveplot(self, value)
      self.saveplot = value;
    end
    function set.sadcp_dtok(self, value)
      self.sadcp_dtok = value;
    end
    function set.use_sadcp(self, value)
      self.use_sadcp = value;
    end
    function set.plot_range(self, value)
      self.plot_range = value;
    end
    function set.edit_mask_dn_bins(self, value)
      self.edit_mask_dn_bins = value;
    end
    function set.edit_mask_up_bins(self, value)
      self.edit_mask_up_bins = value;
    end
    function set.edit_sidelobes(self, value)
      self.edit_sidelobes = value;
    end
    function set.edit_spike_filter_max_curv(self, value)
      self.edit_spike_filter_max_curv = value;
    end
    function set.edit_PPI(self, value)
      self.edit_PPI = value;
    end
    function set.edit_PPI_layer_thickness(self, value)
      self.edit_PPI_layer_thickness = value;
    end
    function set.edit_PPI_max_hab(self, value)
      self.edit_PPI_max_hab = value;
    end
    function set.edit_skip_ensembles(self, value)
      self.edit_skip_ensembles = value;
    end
    function set.detect_asynchronous(self, value)
      self.detect_asynchronous = value;
    end
    function set.zbottom(self, value)
      self.zbottom = value;
    end
    function set.zbottomerror(self, value)
      self.zbottomerror = value;
    end
    function set.edit_mask_last_bin(self, value)
      self.edit_mask_last_bin = value;
    end
    function set.ladcp_cast(self, value)
      self.ladcp_cast = value;
    end
    function set.warnp(self, value)
      self.warnp = value;
    end
    function set.down_sn(self, value)
      self.down_sn = value;
    end
    %     function set.up_sn(self, value)
    %       self.up_sn = value;
    %     end
    function set.dist_up_down(self, value)
      self.dist_up_down = value;
    end
    
  end % end of public methods
  
end % end of class

