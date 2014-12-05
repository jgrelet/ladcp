function l = geterr(dr,d,p,iplot)
% function l = geterr(dr,d,p,iplot)
%
% returns predictions of U_ocean and
% U_ctd on the raw data grid
% 
% CTD velocity
%
% version 0.3  last change 10.11.2012

% replaced finite with isfinite                   GK, 08.11.2009     -->0.2
% fixed problem with uctd extrapolation           GK, 10.11.2012  0.2-->0.3

if nargin<4
  iplot=1; 
end


%
% make sure that we do no interpolate beyond the values we have
%
tim = dr.tim;
uctd1 = dr.uctd;
vctd1 = dr.vctd;
if min(d.time_jul)<tim(1) | max(d.time_jul)>tim(end)
  disp('>   uctd timing interpolation problem')
  disp('>   using constant uctd to extrapolate')
  if min(d.time_jul)<tim(1)
    tim = tim([1,1:end]);
    tim(1) = min(d.time_jul);
    uctd1 = uctd1([1,1:end]);
    vctd1 = vctd1([1,1:end]);
  end
  if max(d.time_jul)<tim(end)
    tim = tim([1:end,end]);
    tim(end) = max(d.time_jul);
    uctd1 = uctd1([1:end,end]);
    vctd1 = vctd1([1:end,end]);
  end
end
uctd = -interp1q(tim',uctd1',d.time_jul');
vctd = -interp1q(tim',vctd1',d.time_jul');


[ib,it] = size(d.ru);

wm = meanmediannan(d.rw,3);
wz = gradient(-d.z,d.time_jul*86400);
l.ru_ctd = meshgrid(uctd,[1:ib])+d.weight*0;
l.rv_ctd = meshgrid(vctd,[1:ib])+d.weight*0;
l.rw_ctd = meshgrid(wm,[1:ib])+d.weight*0;
l.rw_ctd_z = meshgrid(wz,[1:ib])+d.weight*0;
if isfield(d,'wctd')
  l.rw_ctd_p = meshgrid(d.wctd,[1:ib])+d.weight*0;
end

% OCEAN velocity

z = -d.izm+d.ru*0;
dz = diff(d.izm(:,1))';

ii = find(z>=min(dr.z) & z<=max(dr.z));

uoce = interp1q(dr.z,dr.u,z(ii));
voce = interp1q(dr.z,dr.v,z(ii));

[prof,bin] = meshgrid([1:it],[1:ib]);

l.ru_oce = full(sparse(bin(ii),prof(ii),uoce));
l.rv_oce = full(sparse(bin(ii),prof(ii),voce));
l.ru_oce(ib,it) = NaN;
l.rv_oce(ib,it) = NaN;
l.ru_oce = l.ru_oce+d.weight*0;
l.rv_oce = l.rv_oce+d.weight*0;
ii = find(l.ru_oce==0 & l.rv_oce==0);
l.ru_oce(ii) = NaN;
l.rv_oce(ii) = NaN;

% ocean velocity as a function of depth and time

					% ib is number of bins
					% it is number of times (super ensembles)
itm = meshgrid([1:it],[1:ib]);		% each of ib rows of itm contains 1:it

					% d.izm contains for each time (colums),
					% list of absolute depths for each bin
dzdo = mean(abs(diff(d.izm(d.izd,1))));	% dzdo contains sound-speed corrected
					% mean bin length of downlooker at surface
					% NB: at depth, bins are smaller, because
					%     of increased soundspeed!

if length(d.izu)>1			% uplooker bin length
  dzup = mean(abs(diff(d.izm(d.izu,1))));
else
  dzup = dzdo;
end
dz = min([dzdo dzup]);			% dz is min bin length near surface
iz = -(d.izm/dz);			% iz is d.izm with depth coordinate given
					% as distance from surface, measured in 
					% near-surface bin lengths 

					% d.ru contains super-ensemble velocities
					% dr.z contains output depth grid
ii = find(isfinite(d.ru) & iz>0 & iz<max(dr.z)/dz);
					% ii contains indices (valid for d.ru,
					% d.izm, iz, ...) with valid velocities,
					% inside the output depth grid

ij = find( iz>0 & iz<max(dr.z)/dz);	% ij contains same as ii but also for
					% invalid velocities

if abs(dzup-dzdo)>dzup*0.1
  disp(['    Sorry dz not constant, looping over ',...
	int2str(length(ii)),' elements'])
  for j=1:length(ii)
    iiz = ceil(iz(ii(j)));
    iit = itm(ii(j));
    l.u_oce(iiz,iit) = d.ru(ii(j))-l.ru_ctd(ii(j));
    l.v_oce(iiz,iit) = d.rv(ii(j))-l.rv_ctd(ii(j));
    l.w_oce(iiz,iit) = d.rw(ii(j))-l.rw_ctd(ii(j));
    l.w_oce_z(iiz,iit) = d.rw(ii(j))-l.rw_ctd_z(ii(j));
    if isfield(l,'rw_ctd_p')
      l.w_oce_p(iiz,iit) = d.rw(ii(j))-l.rw_ctd_p(ii(j));
    end
    if isfield(d,'tg')
      l.tg_oce(iiz,iit) = d.tg(ii(j));
    end

    l.u_ocean(iiz,iit) = l.ru_oce(ii(j));
    l.v_ocean(iiz,iit) = l.rv_oce(ii(j));

    l.u_adcp(iiz,iit) = d.ru(ii(j));
    l.v_adcp(iiz,iit) = d.rv(ii(j));
  end
else % uplooker and downlooker bin sizes are equal
  l.u_oce = full(sparse(ceil(iz(ii)),itm(ii),d.ru(ii)-l.ru_ctd(ii)));
  l.v_oce = full(sparse(ceil(iz(ii)),itm(ii),d.rv(ii)-l.rv_ctd(ii)));
  l.w_oce = full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd(ii)));
  l.w_oce_z = full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd_z(ii)));
  if isfield(l,'rw_ctd_p')
    l.w_oce_p = full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd_p(ii)));
  end
  if isfield(d,'tg')
    l.tg_oce = full(sparse(ceil(iz(ij)),itm(ij),d.tg(ij)));
  end

  l.u_ocean = full(sparse(ceil(iz(ii)),itm(ii),l.ru_oce(ii)));
  l.v_ocean = full(sparse(ceil(iz(ii)),itm(ii),l.rv_oce(ii)));

  l.u_adcp = full(sparse(ceil(iz(ii)),itm(ii),d.ru(ii)));
  l.v_adcp = full(sparse(ceil(iz(ii)),itm(ii),d.rv(ii)));
end

ik = find(l.u_oce==0 & l.v_oce==0);
l.u_oce(ik) = NaN;
l.v_oce(ik) = NaN;
l.w_oce(ik) = NaN;
l.w_oce_z(ik) = NaN;
if isfield(l,'rw_ctd_p')
  l.w_oce_p(ik) = NaN;
end
l.u_adcp(ik) = NaN;
l.v_adcp(ik) = NaN;
if isfield(d,'tg')
  ik = find(l.tg_oce==0);
  l.tg_oce(ik) = NaN;
end

[lz,lt] = size(l.u_oce);
l.itv = [1:lt];

l.z_oce = ([1:lz]-.5)*dz;
l.u_oce_m = nmean(l.u_oce');
l.v_oce_m = nmean(l.v_oce');

l.u_oce_s = nstd(l.u_oce');
l.v_oce_s = nstd(l.v_oce');

l.ru_err = d.ru-l.ru_oce-l.ru_ctd;
l.rv_err = d.rv-l.rv_oce-l.rv_ctd;

l.izm = d.izm;

[lz,lt] = size(l.ru_err);
l.itv2 = [1:lt];

if iplot

  % blank out shallow/deep estimates
  ii = find(iz<0 | iz>max(dr.z)/dz);
  d.ru(ii) = nan;
  d.rv(ii) = nan;


  figure(2)
  clf
  orient tall
  colormap([[1 1 1]; jet(128)])

  subplot(231)
  ib = [1:size(l.ru_err,1)];
  ib = ib-length(d.izu);
  pcolorn(l.itv2,-ib,l.ru_err) 
  fac = nmean(l.u_oce_s);
  fac=  max([fac,1e-2]);
  caxis([-3 3]*fac)
  colorbar
  title(['U-err std:',num2str(nmean(nstd(l.ru_err')))])
  ylabel('BINS')
   
  subplot(232)
  pcolorn(l.itv,-l.z_oce,l.u_oce)
  if isfield(dr,'zbot')
    hold on
    plot(d.z-d.hbot,'.k')
    ax = axis;
    ax(3) = nmin([d.z-d.hbot,ax(3)]);
    axis(ax);
  end
  colorbar
  ca = caxis;
  ax = axis;
  title('U_{ocean}')

  subplot(233)
  iz = d.izu;
  plot(d.ru(iz,:)-l.ru_ctd(iz,:),l.izm(iz,:),'k.','markersize',3)
  hold on
  iz = d.izd;
  plot(d.ru(iz,:)-l.ru_ctd(iz,:),l.izm(iz,:),'b.','markersize',3)
  plot(dr.u,-dr.z,'-r')
  plot([l.u_oce_m'+l.u_oce_s', l.u_oce_m'-l.u_oce_s'],-l.z_oce','-g')
  if isfield(dr,'u_shear_method')& length(dr.z)==length(dr.u_shear_method)
    plot(dr.u_shear_method+nmean(dr.u),-dr.z,'-k')
  end
  ax = [ca,ax(3:4)];
  axis(ax)
  grid
  title('U (.k) up-looker (.b) down-looker')

  subplot(234)
  pcolorn(l.itv2,-ib,l.rv_err) 
  fac = nmean(l.v_oce_s);
  fac = max([fac,1e-2]);
  caxis([-3 3]*fac)
  colorbar
  title(['V-err std:',num2str(nmean(nstd(l.rv_err')))])
  ylabel('BINS')
  xlabel('super ensemble')

  subplot(235)
  pcolorn(l.itv,-l.z_oce,l.v_oce) 
  hold on
  if isfield(dr,'zbot')
    plot(d.z-d.hbot,'.k')
    ax = axis;
    ax(3) = nmin([d.z-d.hbot,ax(3)]);
    axis(ax);
  end
  colorbar
  ca = caxis;
  ax = axis;
  title('V_{ocean}')
  xlabel('super ensemble')

  subplot(236)
  iz = d.izu;
  plot(d.rv(iz,:)-l.rv_ctd(iz,:),l.izm(iz,:),'k.','markersize',3)
  hold on
  iz = d.izd;
  plot(d.rv(iz,:)-l.rv_ctd(iz,:),l.izm(iz,:),'b.','markersize',3)
  plot(dr.v,-dr.z,'-r')
  plot([l.v_oce_m'+l.v_oce_s', l.v_oce_m'-l.v_oce_s'],-l.z_oce','-g')
  if isfield(dr,'v_shear_method')& length(dr.z)==length(dr.v_shear_method)
    plot(dr.v_shear_method+nmean(dr.v),-dr.z,'-k')
  end
  ax = [ca,ax(3:4)];
  axis(ax)
  grid
  title('V (-k)shear (-r)inverse (-g)+/-std')
  xlabel('velocity [m/s]')
  
  streamer([dr.name,'  Figure 3']);
  hgsave('tmp/3')

  % reset colormap
  figure(2)
  colormap(jet(128))
 
end
