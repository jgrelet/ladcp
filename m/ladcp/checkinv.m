function p = checkinv(dr,de,der,p,ps,values)
% function p = checkinv(dr,de,der,p,ps,values)
% check inversion for consistency 
%
% version 0.2	last change 02.06.2007

% Martin Visbeck 2004
% G.Krahmann, 2005

% changed bar plot, problems with 2006b		GK, Jun 2007	?->0.2


%
% general function info
%
disp(' ')
disp('CHECKINV:  displaying check values ')


uerrm = median(dr.uerr);
unoise = nmean(nstd([der.ru_err ; der.rv_err])) ;

disp('    The following values are given in [m/s] ')
disp(['    Velocity profile error: ',num3str(uerrm,6,3),...
       '  should be about noise: ',num3str(unoise,6,3)])

% check bottom track
if isfield(dr,'zbot') 
  if isfield(de,'bvels') 

    % std of missmatch between bottom track and U_ctd
    ubso = ([real(de.bvel)+dr.uctd]);
    vbso = ([imag(de.bvel)+dr.vctd]);
    uvbso = nmean(sqrt(ubso.^2+vbso.^2));

    % std of bottom track 
    uvbsi = nmean(de.bvels);
  else
    disp('    Do not know how uncertain the bottom track is')
    uvbso = nan;
    uvbsi = nan;
  end

  disp(['    Check bottom track rms: ',num3str(uvbso,6,3),...
       '  should be smaller than ',num3str(uvbsi,6,3),' / ',num3str(ps.botfac,6,3)])

end

% check SADCP
if isfield(dr,'z_sadcp')

  % std of missmatch between SADCP and LADCP velocity profile
  us = interp1(dr.z,dr.u,dr.z_sadcp);
  vs = interp1(dr.z,dr.v,dr.z_sadcp);
  usso = ([us-dr.u_sadcp]);
  vsso = ([vs-dr.v_sadcp]);
  uvsso = nmean(sqrt(usso.^2+vsso.^2));
 
  % given std of SADCP profile
  uvssi = nmedian(dr.uerr_sadcp);

  disp(['    Check SADCP        rms: ',num3str(uvsso,6,3),...
       '  should be smaller than ',num3str(uvssi,6,3),' / ',num3str(ps.sadcpfac,6,3)])

end

% check GPS
if isfield(dr,'uship')

  % difference of mean ships drift from LADCP and GPS 
  dtiv = gradient(dr.tim);
  us = sum(dtiv.*dr.uship)/sum(dtiv);
  vs = sum(dtiv.*dr.vship)/sum(dtiv);
  uvgso = sqrt((us-values.uship).^2 + (vs-values.vship).^2);
 
  % computed uncertainty of GPS derived ships speed
  uvgsi = ps.barvelerr;

  disp(['    GPS-LADCP ship spd diff:',num3str(uvgso,6,3),...
       '  should be smaller than ',num3str(uvgsi,6,3),' / ',num3str(ps.barofac,6,3)])
end


%
% plot inversion solution weights
%
if isfield(de,'type_constraints')
 
  ic = find(sum([de.ctd_constraints,de.ocean_constraints]')~=0);  

  figure(2)
  clf
  orient tall

  ic2 = 3*length(ic);
  col = jet(ic2);
  ic3 = [1:length(ic)]*3;
  ic3 = [ic3(1:2:end),fliplr(ic3(2:2:end))];
  colormap(col(ic3,:));
 
  subplot(211)
  simple_bar(dr.z',de.ocean_constraints(ic,:),de.type_constraints(ic,:)) 
  ylabel('sum of weights')
  title('ocean velocity constraints')
  xlabel('depth [m]')

  subplot(212)
  simple_bar([1:size(de.ctd_constraints(ic,:),2)],de.ctd_constraints(ic,:),...
	de.type_constraints(ic,:)) 
  title('CTD velocity constraints')
  ylabel('sum of weights')
  xlabel('super ensemble')

  streamer([p.name,'  Figure 12']);
  hgsave('tmp/12')
end
