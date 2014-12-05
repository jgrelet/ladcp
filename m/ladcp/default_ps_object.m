classdef default_ps_object < dynamicprops
  % class default_ps_object contains file names
  %
  % Gerd Krahmann, Kiel July 2006
  % gkrahmann@ifm-geomar.de
  % Frederic Marin IRD - LEGOS - Noumea - June 2011
  % Frederic.Marin@ird.fr
  % Jacques Grelet IRD US191 - Plouzane - May/Dec 2014
  % jacques.grelet@ird.fr
  %
  % $Id$
  %
  
  
  properties (Access = public)
    % Process data using shear based method
    % compute shear based solution
    % ps.shear=2  ; use super ensemble
    % ps.shear=1  ; use raw data
    shear = 1
    
    
    % decide how to weigh data
    % 1 : use super ensemble std
    % 0 : use correlation based field
    std_weight = 1
    
    % Weight for the barotropic constraint
    barofac = 1
    
    % Weight for the bottom track constraint
    botfac = 1
    
    % Process up and down cast seperately
    down_up = 1
    
    % Depth resolution for final profile
    %	default one bin length
    % ps=setdefv(ps,'dz',medianan(abs(diff(di.izm(:,1)))));
    
    % Smoothing of the final profile
    smoofac = 0
    
    % Request that shears are small  (experts only)
    smallfac = [1 0]
    
    % weight bottom track data with distance of bottom
    %  use Gaussian with offset (btrk_weight_nblen(1) * bin)
    %  and width (btrk_weight_nblen(2) * bin)
    %  one might set this to [15 5] to reduce the weight of close by bottom track data
    btrk_weight_nblen = [0 0]
    
    % Weight for SADCP data
    % ps.sadcpfac=1 about equal weight for SDACP profile
    sadcpfac = 1
    
    % average over data within how many standard deviations
    shear_stdf = 2
    
    % the minimum weight a bin must have to be accepted for shear
    shear_weightmin = 0.1
    
    % super ensemble velocity error
    % try to use the scatter in W to get an idea of the "noise"
    % in the velocity measurement
    %This is a bit of code used in GETINV.m
    % nmax=min(length(di.izd),7);
    % sw=nstd(di.rw(di.izd(1:nmax),:)); ii=find(sw>0);
    % sw=medianan(sw(ii))/tan(p.Beam_angle*pi/180);
    % ps=setdefv(ps,'velerr',max([sw,0.02]));
    %
    % velerr = 0.02
    
    
    % How to solve the inverse
    %     ps.solve = 0  Cholseky transform
    %              = 1  Moore Penrose Inverse give error for solution
    solve = 1
    
    
    % Threshold for minimum weight, data with smaller weights
    %  	will be ignored
    weightmin = 0.05
    
    % Change the weights by
    %	weight=weight^ps.weightpower
    weightpower = 1
    
    % Remove 1% of outliers (values in a superensemble deviating most from
    % the inversion result) after solving the inversion
    % ps.outlier defines how many times 1% is removed and the
    % inversion is recalculated.
    outlier = 1
    
    
    % Weight for the cable drag constraint
    % only for experts
    dragfac = 0
    drag_tilt_vel = 0.5
    
    % average ctdvel back in time to take the inertia of wire into account
    % smooth over max of 20 minutes for depth ~ 2000m
    drag_lag_tim = 20
    drag_lag_depth = 2000
    
    % getinv
    up_dn_looker
    ucorr
    vcorr
    
    
  end
  
  properties (Access = public, Hidden)
  end
  
  properties (Access = public, Dependent = true)
  end
  
  
  %% public methods
  methods
    % constructor
    function self = default_ps_object
      
    end
    
  end % end of public methods
  
end % end of class

