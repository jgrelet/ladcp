function [da]=savearch(values,dr,d,p,ps,f)
% function [da]=savearch(values,dr,d,p,ps,f)
%
% store LADCP result LADCP archive format
%

% Attributes For Structure dr
att.name.long_name      	  = 'Cast ID';
att.date.long_name      	  = 'Date';
att.date.units          	  = 'Y M D H M S';
att.lat.long_name       	  = 'Latitude';
att.lat.units           	  = 'Degree North';
att.lon.long_name       	  = 'Longitude';
att.lon.units           	  = 'Degree East';
att.zbot.long_name      	  = 'Bottom Referenced Profile Depth';
att.zbot.units          	  = 'm';
att.ubot.long_name      	  = 'Bottom Referenced Profile U';
att.ubot.units          	  = 'm/s';
att.vbot.long_name      	  = 'Bottom Referenced Profile V';
att.vbot.units          	  = 'm/s';
att.uerrbot.long_name   	  = 'Bottom Referenced Profile Velocity Error';
att.uerrbot.units       	  = 'm/s';
att.z_sadcp.long_name   	  = 'SADCP Profile Depth';
att.z_sadcp.units   		  = 'm';
att.u_sadcp.long_name   	  = 'SADCP Profile U';
att.u_sadcp.units   		  = 'm/s';
att.v_sadcp.long_name   	  = 'SADCP Profile V';
att.v_sadcp.units   		  = 'm/s';
att.uerr_sadcp.long_name   	  = 'SADCP Profile Velocity Error';
att.uerr_sadcp.units   		  = 'm/s';
att.z.long_name                   = 'Depth';
att.z.units                       = 'Meters';
att.u.long_name                   = 'U';
att.u.units                       = 'm/s';
att.v.long_name                   = 'V';
att.v.units                       = 'm/s';
att.uerr.long_name                = 'Velocity Error';
att.uerr.units                    = 'm/s';
att.nvel.long_name                = 'LADCP number of ensembles per bin';
att.ubar.long_name                = 'LADCP U Barotropic';
att.ubar.units              	  = 'm/s';
att.vbar.long_name                = 'LADCP V Barotropic';
att.vbar.units              	  = 'm/s';
att.tim.long_name 		  = 'Station Time Series';
att.tim.units 			  = 'Julian Days';
att.tim_hour.long_name 		  = 'Station Time Series';
att.tim_hour.units 		  = 'Hour of Day';
att.shiplon.long_name 		  = 'Longitude';
att.shiplon.units     		  = 'Degree East';
att.shiplat.long_name 		  = 'Latitude';
att.shiplat.units     		  = 'Degree North';
att.xship.long_name   		  = 'Ship Position E';
att.xship.units       		  = 'Meters East';
att.yship.long_name   		  = 'Ship Position N';
att.yship.units       		  = 'Meters North';
att.uship.long_name   		  = 'Ship Velocity U';
att.uship.units       		  = 'm/s';
att.vship.long_name   		  = 'Ship Velocity V';
att.vship.units       		  = 'm/s';
att.zctd.long_name    		  = 'Depth of CTD';
att.zctd.units        		  = 'm';
att.wctd.long_name    		  = 'CTD Velocity W';
att.wctd.units        		  = 'm/s';
att.uctd.long_name    		  = 'CTD Velocity U';
att.uctd.units        		  = 'm/s';
att.vctd.long_name    		  = 'CTD Velocity V';
att.vctd.units        		  = 'm/s';
att.uctderr.long_name 		  = 'CTD Velocity Error';
att.uctderr.units     		  = 'm/s';
att.xctd.long_name    		  = 'CTD Position Relative to Start E';
att.xctd.units        		  = 'm';
att.yctd.long_name    		  = 'CTD Position Relative to Start N';
att.yctd.units        		  = 'm';
att.range.long_name      	  = 'ADCP total range of data';
att.range.units          	  = 'm';
att.range_do.long_name   	  = 'ADCP down looking range of data';
att.range_do.units       	  = 'm';
att.range_up.long_name   	  = 'ADCP up looking range of data';
att.range_up.units       	  = 'm';
att.ts.long_name         	  = 'ADCP echo amplitude profile bin 1';
att.ts.units                      = 'dB';
att.ts_out.long_name              = ['ADCP echo amplitude ' ...
		                          'profile last down bin'];
att.ts_out.units                  = 'dB';
att.u_do.long_name       	  = 'LADCP down only profile U';
att.u_do.units           	  = 'm/s';
att.v_do.long_name       	  = 'LADCP down only profile V';
att.v_do.units           	  = 'm/s';
att.u_up.long_name       	  = 'LADCP up only profile U';
att.u_up.units           	  = 'm/s';
att.v_up.long_name       	  = 'LADCP up only profile V';
att.v_up.units           	  = 'm/s';
att.p.long_name          	  = 'Pressure';
att.p.units              	  = 'dBar';
att.ensemble_vel_err.long_name    = 'ADCP ensemble velocity error';
att.ensemble_vel_err.units        = 'm/s';
att.u_shear_method.long_name      = 'LADCP shear method profile U';
att.u_shear_method.units 	  = 'm/s';
att.v_shear_method.long_name      = 'LADCP shear method profile V';
att.v_shear_method.units 	  = 'm/s';
att.ctd_t.long_name      	  = 'CTD profile temperature';
att.ctd_t.units          	  = 'Degree C';
att.ctd_s.long_name      	  = 'CTD profile salinity';
att.ctd_s.units          	  = 'psu';
att.ctd_ss.long_name     	  = 'CTD profile sound speed';
att.ctd_ss.units         	  = 'm/s';
att.ctd_N2.long_name     	  = 'CTD profile stability';
att.ctd_N2.units         	  = '1/s^2';




g=gregoria(d.time_jul(1));
p=setdefv(p,'ref_year',g(1));
year0=julian([p.ref_year,0,0,0,0,0]);


 da.GEN_Velocity_Units                = 'm/s';
 da.GEN_LADCP_station                 = p.ladcp_station;
 da.GEN_LADCP_cast                    = p.ladcp_cast;
 da.GEN_Profile_start_decimal_day     = d.time_jul(1)-year0;

[m,ii]=min(d.z);
 da.GEN_Profile_bottom_decimal_day    = d.time_jul(ii)-year0; 
 da.GEN_Profile_end_decimal_day       = d.time_jul(end)-year0; 

 da.GEN_Profile_start_longitude       = p.poss(3)+p.poss(4)/60; 
 da.GEN_Profile_end_longitude         = p.pose(3)+p.pose(4)/60; 
 da.GEN_Profile_start_latitude        = p.poss(1)+p.poss(2)/60; 
 da.GEN_Profile_end_latitude          = p.pose(1)+p.pose(2)/60; 

 da.GEN_Ocean_depth_m                 = round(p.zbottom);
 da.GEN_Profile_max_depth_m           = round(values.maxdepth);
 da.GEN_Magnetic_deviation_deg        = values.magdev;

 da.BAR_ref_U                     = dr.ubar; 
 da.BAR_ref_V                     = dr.vbar;
 da.BAR_ref_error                 = 2*p.nav_error/p.dt_profile; 
 da.BAR_tide_U                    = NaN;
 da.BAR_tide_V                    = NaN;
 da.INPUT_SADCP_profile_avail           = isfield(dr,'u_sadcp'); 
 da.INPUT_Pegasus_profile_avail         = 0; 
 da.INPUT_Bottom_track_profile_avail    = (p.btrk_used>0 & isfield(dr,'zbot'));
 da.INPUT_Nav_time_series_avail         = length(f.nav)>1 ; 
 da.INPUT_CTD_time_series_avail         = length(f.ctdtime)>1 ; 

 da.BAR_ref_descr                 = [];
 if ps.barofac>0, da.BAR_ref_descr= [da.BAR_ref_descr,'[Navigation]']; end
 if ps.botfac>0 & isfield(dr,'zbot'), 
        da.BAR_ref_descr= [da.BAR_ref_descr,'[Bottom-track]']; 
 end
 if ps.dragfac>0, da.BAR_ref_descr= [da.BAR_ref_descr,'[CTD-drag]']; end
 if ps.sadcpfac>0, da.BAR_ref_descr= [da.BAR_ref_descr,'[SADCP]']; end
 if length(da.BAR_ref_descr) < 1; da.BAR_ref_descr='[NA]'; end
 da.BAR_tide_model=                 '[NA]';

% down instrument information

 da.LADCP_dn_hard_type=['[RDI-',int2str(d.down.Frequency),'BB]'];
 if round(d.down.Frequency)==300
  da.LADCP_dn_hard_type=['[RDI-',int2str(d.down.Frequency),'WH]'];
 end
 if isfield(d.down,'NarrowBand')
  da.LADCP_dn_hard_type=['[RDI-',int2str(d.down.Frequency),'NB]'];
 end
 da.LADCP_dn_hard_freq_kHz            = d.down.Frequency;
 da.LADCP_dn_hard_SN                  = p.down_sn;
 da.LADCP_dn_hard_FV                  = d.down.Firm_Version;
 da.LADCP_dn_hard_TNO                 = '[convex4]';
 da.LADCP_dn_hard_beam_ang_deg        = d.down.Beam_angle;
 da.LADCP_dn_hard_comp_type           = '[RDI]';
 da.LADCP_dn_hard_general_comments    = '  ';   

 da.LADCP_dn_conf_blank_intvl_m       = d.down.Blank/100;
 da.LADCP_dn_conf_bin_len_m           = d.down.Cell_length/100; 
 da.LADCP_dn_conf_pulse_len_m         = d.down.Pulse_length/100;
 da.LADCP_dn_conf_number_bins         = length(d.izd); 
 da.LADCP_dn_conf_ping_stagr          = '[NA]';
 da.LADCP_dn_conf_ping_trns_intvl_sec = d.down.Time_Pings;
 da.LADCP_dn_conf_number_pings        = d.down.Pings_per_Ensemble;
 da.LADCP_dn_conf_vel_ambiguity       = p.ambiguity;
 da.LADCP_dn_conf_single_ping_acc     = d.down.Single_Ping_Err;
 da.LADCP_dn_xmit_cur            = values.xmc(1);
 da.LADCP_dn_xmit_vol            = values.xmv(1);
 da.LADCP_dn_xmit_pings          = values.nping_total(1);
 da.LADCP_dn_beam_range          = p.dn_range;
 if d.down.Coordinates==3
  da.LADCP_dn_conf_coord_system='[earth]';
 elseif d.down.Coordinates==1
  da.LADCP_dn_conf_coord_system='[beam]';
 else
  da.LADCP_dn_conf_coord_system='[unknown]';
 end
 da.LADCP_dn_conf_bottom_trkr=p.btrk_used;
 if isfinite(p.zbottom) & p.btrk_used>0 & d.down.Up==0
  da.LADCP_dn_btrk_u_bias = p.btrk_u_bias;
  da.LADCP_dn_btrk_v_bias = p.btrk_v_bias;
  da.LADCP_dn_btrk_u_std =  p.btrk_u_std;
  da.LADCP_dn_btrk_v_std =  p.btrk_v_std;
 end
 da.LADCP_dn_conf_general_comments='   ';

if isfield(d,'up')
% up instrument information

 da.LADCP_up_hard_type=['[RDI-',int2str(d.up.Frequency),'BB]'];
 if round(d.up.Frequency)==300
  da.LADCP_up_hard_type=['[RDI-',int2str(d.up.Frequency),'WH]'];
 end
 da.LADCP_up_hard_freq_kHz            = d.up.Frequency;
 da.LADCP_up_hard_SN                  = p.up_sn;
 da.LADCP_up_hard_FV                  = d.up.Firm_Version;
 da.LADCP_up_hard_TNO                 = '[convex4]';
 da.LADCP_up_hard_beam_ang_deg        = d.up.Beam_angle;
 da.LADCP_up_hard_comp_type           = '[RDI]';
 da.LADCP_up_hard_general_comments    ='  ';   

 da.LADCP_up_conf_blank_intvl_m       = d.up.Blank/100;
 da.LADCP_up_conf_bin_len_m           = d.up.Cell_length/100; 
 da.LADCP_up_conf_pulse_len_m         = d.up.Pulse_length/100;
 da.LADCP_up_conf_number_bins         = length(d.izd); 
 da.LADCP_up_conf_ping_stagr          = '[NA]';
 da.LADCP_up_conf_ping_trns_intvl_sec = d.up.Time_Pings;
 da.LADCP_up_conf_number_pings        = d.up.Pings_per_Ensemble;
 da.LADCP_up_conf_vel_ambiguity       = p.ambiguity;
 da.LADCP_up_conf_single_ping_acc     = d.up.Single_Ping_Err;
 da.LADCP_up_xmit_cur                 = values.xmc(2);
 da.LADCP_up_xmit_vol                 = values.xmv(2);
 da.LADCP_up_xmit_pings               = values.nping_total(2);
 da.LADCP_up_beam_range               = p.up_range;

 if p.rotup2down==2
  da.LADCP_up_compass='[velocity-match]';
 elseif p.rotup2down==1
  da.LADCP_up_compass='[down-compass]';
 else
  da.LADCP_up_compass='[up-compass]';
 end
 if d.up.Coordinates==3
  da.LADCP_up_conf_coord_system='[earth]';
 elseif d.up.Coordinates==1
  da.LADCP_up_conf_coord_system='[beam]';
 else
  da.LADCP_up_conf_coord_system='[unknown]';
 end
 da.LADCP_up_conf_general_comments='   ';
end

 da.GEN_LADCP_ensemble_time_mean_sec=mean(diff(d.time_jul*24*3600));
 da.GEN_LADCP_ensemble_time_std_sec=std(diff(d.time_jul*24*3600));
 da.GEN_conf_general_comments ='  ';

 da.GEN_Matlab_version=version;
 da.GEN_Processing_personnel= p.whoami;
 da.GEN_Processing_date=date;
 da.GEN_Proc_methodology= '[inverse]';
 da.GEN_Software_orig= p.software;
 da.GEN_Sound_sp_calc = values.GEN_Sound_sp_calc;

 % da.Depth_source:  (choose one) <w> <w&Pmax> <measured P (CTD)>
 %                <measured P (other)> <NA> <unconfirmed>
 da.GEN_Depth_source=  '[w]'; 
 if values.ladcpdepth==2
  da.GEN_Depth_source=  '[w&surface&bottom]';
 end
 if isfinite(p.zpar(2))
  da.GEN_Depth_source(end)=[];
  da.GEN_Depth_source=  [da.GEN_Depth_source,'&Pmax]']; 
 end
 if ~isempty(d.ctdtime_data)
  da.GEN_Depth_source=  '[measured P (CTD)]';
 end

jok = cumprod(size(find(~isnan(d.rw))));
j = cumprod(size(find(isnan(d.re) & ~isnan(d.rw))));
 da.GEN_Percent_3beam=round(j*100/jok); 
 da.GEN_Editing_parm_descr=p.outlier;
 da.GEN_Inverse_weight_bottom=ps.botfac;
 da.GEN_Inverse_weight_navigation=ps.barofac;
 da.GEN_Inverse_weight_smooth=ps.smoofac;
 da.GEN_Proc_general_comments='  ';


if exist([f.res,'.log'],'file')~=0
 diary off
 id=fopen([f.res,'.log']);
 da.LOG_Inverse_log=setstr(fread(id))';
 fclose(id);
end

% save to matlab file 
if p.savemat
  disp(['    save ',[f.res,'.ladcp.mat'],' da dr att'])
  save6([f.res,'.ladcp.mat'],'da','dr','att')
end

% save to netcdf file if you have NETCDF libary
if exist('netcdf')~=0 & p.savecdf
 dr.tim=dr.tim-year0;
 disp(['    saving level-II results in netcdf-file: ',f.res,'.nc'])
 ladcp2cdf([f.res,'.nc'],dr,da,p,ps,f,att)
else
 disp('>   Please install NETCDF libary to get the netcdf LADCP archive output')
end 
