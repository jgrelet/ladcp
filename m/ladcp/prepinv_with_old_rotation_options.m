function [p,d,messages]=prepinv(messages,d,p,dr,values)
% function [p,d,messages]=prepinv(messages,d,p,dr,values)
%
% prepare for inverse solver
%
% - average velocities within depth range
%
% version 0.3	last change 28.05.2011

% 
%  Martin Visbeck, LDEO, 6/10/99
%  Gerd Krahmann, IFM-GEOMAR, September 2007

% fix problem with trusted_i > # of bins	GK, Sep 2007	0.1-->0.2
% bug in tilt_weight application          GK, 28.05.2011  0.2-->0.3


%
% general function info
%
disp(' ')
disp('PREPINV:  prepare data for inversion')


%
% remove bottom track values outside a distance-window from the bottom
%
if isfinite(p.zbottom)

  ind = find(abs(d.z+p.zbottom)>max(p.btrk_range) & isfinite(d.bvel(:,1)') );
  if length(ind)>0
    disp(['    Discarded ',int2str(length(ind)),...
	' bottom tracks velocities because of height above bottom > ',...
	int2str(max(p.btrk_range))])
    d.bvel(ind,:) = NaN;
    d.bvels(ind,:) = NaN;
  end

  ind = find(abs(d.z+p.zbottom)<min(p.btrk_range) & isfinite(d.bvel(:,1)') );
  if length(ind)>0
    disp(['    Discarded ',int2str(length(ind)),...
	' bottom tracks velocities because of height above bottom < ',...
	int2str(min(p.btrk_range))])
    d.bvel(ind,:) = NaN;
    d.bvels(ind,:) = NaN;
  end

  ind = find(abs(d.z+p.zbottom-d.hbot)> 100 );
  if length(ind)>0
    disp(['    Discarded ',int2str(length(ind)),...
	    ' bottom distances because of depth difference > 100'])
    d.hbot(ind) = NaN;
  end

else

  d.bvel(:,:) = NaN;
  d.bvels(:,:) = NaN;

end


%
% reduce weight for large tilts
%
if p.tilt_weight>0 & ~isfield(d,'tilt_weight')
  disp(['    Reducing weights for tilts larger than ',...
      num2str(p.tilt_weight),' degree'])
  fac = 1.-tanh(d.tilt/p.tilt_weight)/2;
  oldweight = d.weight;
%  d.weight = d.weight.*meshgrid(fac,d.weight(:,1));  <-- this appears
%  buggy , GK 28.05.2011
  d.weight = d.weight.*(ones(size(d.weight,1),1)*fac);
  d.tilt_weight = p.tilt_weight;
end


if 1
%%%%%%%%%%%%%%%%%%%% rotation

% prepare for heading averaging
u1d = exp(-1i*(d.hdg(1,:))*pi/180); 
if values.up==1

  u1u = exp(-1i*(d.hdg(2,:))*pi/180); 

  % get mean heading offset from COMPASS comparison
  hdg_offset = compoff(u1d,u1u);

  % hdg_offset=angle(u1d/u1u)*180/pi;
  p.up_dn_comp_off = hdg_offset;
  disp(['    Mean heading offset based on compass is ',num2str(hdg_offset),' deg'])

  % check tilt sensors
  % rotate by compass offset
  diary off
  hdg_offset2 = fminsearch('checktilt',0,[],...
    [d.rol(2,:);d.pit(2,:);d.rol(1,:);d.pit(1,:)]);
  diary on
  disp(['    Mean heading offset based on tilt is ',num2str(hdg_offset2),' deg'])
  p.up_dn_pit_rol_comp_off = hdg_offset2;
  [d.rol(3,:),d.pit(3,:)] = uvrot(d.rol(2,:),d.pit(2,:),-hdg_offset2);
  p.up_dn_rol_off = mean(d.rol(1,:)-d.rol(3,:));
  p.up_dn_pit_off = mean(d.pit(1,:)-d.pit(3,:));

  % plot compass and tilt meter differences
  figure(2)
  clf
  orient tall

  % plot compass comparison between up and down instrument
  da = -angle(u1d)+angle(u1u*exp(1i*hdg_offset*pi/180));
  ii = find(da>pi);
  da(ii) = da(ii)-2*pi;
  ii = find(da<-pi);
  da(ii) = da(ii)+2*pi;
  subplot(311)
  plot(d.hdg(1,:),da*180/pi,'.')
  ylabel('Heading difference up-down')
  title([' Heading offset : ',num2str(hdg_offset)])
  xlabel('Heading down')
  grid
  dhmax = sort(abs(da*180/pi));
  dhmax = dhmax(fix(end*0.95));
  if dhmax > 15 & ~isfield(d,'hrot')
    warn = ['>   Large up-down compass difference: ',num2str(dhmax)];
    messages.warnp = strvcat(messages.warnp,warn);
    disp(warn)
  end
  dhmax = max([10 dhmax*1.3]); 
  axis([0 360 [-1 1]*dhmax])

  % plot tilt meter difference
  tiltmax = p.tiltmax(1);
  subplot(312)
  pa = d.pit(1,:)-d.pit(3,:)-p.up_dn_pit_off;
  plot(d.pit(1,:),pa,'.')
  ylabel('Pitch difference up-down')
  xlabel('Pitch of down instrument')
  title([' Pitch offset : ',num2str(p.up_dn_pit_off)])
  grid;
  axis([[-1 1]*tiltmax [-1 1]*5])
 
  subplot(313)
  ra = d.rol(1,:)-d.rol(3,:)-p.up_dn_rol_off;
  plot(d.rol(1,:),ra,'.')
  ylabel('Roll difference up-down')
  xlabel('Roll of down instrument')
  title(['Roll offset : ',num2str(p.up_dn_rol_off)])
  grid;
  axis([[-1 1]*tiltmax [-1 1]*5])

%  streamer([p.name,' Figure 6']);
  suplabel([p.name,' Figure 6'],'t');
  hgsave('tmp/6')




  % offset upward looking ADCP to downward looking ADCP
%  if (p.offsetup2down~=0 & length(d.zd)~=length(d.ru(:,1)) & exist('dr')==1)
  if (p.offsetup2down~=0 & length(d.zd)~=length(d.ru(:,1)) & ~isempty(dr))
    if p.rotup2down==2
      % will not rotate to match velocities and correct offset
      p.rotup2down = 1;
    end
    disp('    Removing first guess ocean velocity from raw data')
    % OCEAN velocity
    [ib,it] = size(d.ru);
    z = -d.izm+d.ru*0;
    ii = find(z>=min(dr.z) & z<=max(dr.z));
    l.uoce(:,1) = interp1q(dr.z,dr.u,z(ii));
    l.uoce(:,2) = interp1q(dr.z,dr.v,z(ii));
    [prof,bin] = meshgrid(1:it,1:ib);
    l.ru = full(sparse(bin(ii),prof(ii),l.uoce(:,1)));
    l.rv = full(sparse(bin(ii),prof(ii),l.uoce(:,2)));
    l.ru(ib,it) = NaN;
    l.rv(ib,it) = NaN;
    ii = find(~(z>=min(dr.z) & z<=max(dr.z)));
    l.ru(ii) = NaN;
    l.rv(ii) = NaN;

    uu = meanmediannan(d.ru(d.izu,:)+1i*d.rv(d.izu,:)...
              -l.ru(d.izu,:)-1i*l.rv(d.izu,:)+d.weight(d.izu,:)*0,2);
    ud = meanmediannan(d.ru(d.izd,:)+1i*d.rv(d.izd,:)...
             -l.ru(d.izd,:)-1i*l.rv(d.izd,:)+d.weight(d.izd,:)*0,2);
    clear l

    ii = find(~isfinite(uu+ud));
    uu(ii) = 0;
    ud(ii) = 0;
    uoff = (ud-uu)*p.offsetup2down; 
    uoffm(d.izu,:) = meshgrid(uoff/2,d.izu);
    uoffm(d.izd,:) = meshgrid(-uoff/2,d.izd);
    d.ru = d.ru+real(uoffm); 
    d.rv = d.rv+imag(uoffm); 
    if isfield(d,'bvel')
      d.bvel(:,1) = d.bvel(:,1)+real(-uoff/2)';
      d.bvel(:,2) = d.bvel(:,2)+imag(-uoff/2)';
    end
    if isfield(d,'uoff')
      d.uoff = d.uoff+uoff;
    else
      d.uoff = uoff;
    end
    % estimate tilt error that could explain difference
    % beam velocity ~ W
    wbeam = (abs(d.wm)*cos(d.down.Beam_angle*pi/180));
    ii = find(wbeam<0.4);
    wbeam(ii) = NaN;
    % uerror for 1 degree tilt error
    du = wbeam/sin((d.down.Beam_angle+0.5)*pi/180)-...
      	wbeam/sin((d.down.Beam_angle-0.5)*pi/180);
    % linear projection to estimate tilt error
    d.tilterr = abs(d.uoff)./abs(du);

    figure(2)
    clf
    orient tall
    subplot(311)
    plot(real(d.uoff))
    uoffav = boxav(d.uoff',20);
    hold on
    plot([0.5:length(uoffav)]*20,real(uoffav),'-r')
    grid
    axis tight
    ax = axis;
    ax(3:4) = [-1 1]*max(nmedian(abs(real(uoffav)))*6,0.1);
    axis(ax)

    title('U offset [m/s]')
    subplot(312)
    plot(imag(d.uoff))
    hold on
    plot([0.5:length(uoffav)]*20,imag(uoffav),'-r')
    grid
    axis tight
    ax = axis;
    ax(3:4) = [-1 1]*max(nmedian(abs(imag(uoffav)))*6,0.1);
    axis(ax)
    title('V offset [m/s]')

    subplot(313)
    plot(d.tilterr)
    hold on
    tilterrav = boxav(d.tilterr',20);
    plot([0.5:length(tilterrav)]*20,tilterrav,'-r')
    axis tight
    ax=axis;
    ax(3:4) = [0 1]*nmedian(d.tilterr)*6;
    if isfinite(sum(ax))
      ax(4)=max(3,ax(4));
      axis(ax)
    end
    grid
    title('Tilt error [degree] consistent with offset')
 
    streamer([p.name,' Figure 10']);
    hgsave('tmp/10')

    disp('    Adjusted for velocity offset in up and down looking ADCP')

  end

  % check if rotated already

  % rotate upward looking ADCP to downward looking ADCP
  if (p.rotup2down~=0 & length(d.zd)~=length(d.ru(:,1)))
    da = angle(u1d)-angle(u1u*exp(-1i*hdg_offset*pi/180));
    ii = find(da>pi); 
    da(ii)=da(ii)-2*pi; 
    ii=find(da<-pi); 
    da(ii)=da(ii)+2*pi;
    d.diff_hdg = da;

    u1uc = exp(-1i*(d.hdg(2,:)-hdg_offset)*pi/180);
    hrotcomp = angle(u1uc./u1d)*180/pi;
    d.rot_comp = hrotcomp;

    % get mean heading offset from velocities
    if ~isempty(dr)
      disp('    Removing first guess ocean velocity from raw data')
      % OCEAN velocity
      [ib,it] = size(d.ru);
      z = -d.izm+d.ru*0;
      ii = find(z>=min(dr.z) & z<=max(dr.z));
      l.uoce(:,1) = interp1q(dr.z,dr.u,z(ii));
      l.uoce(:,2) = interp1q(dr.z,dr.v,z(ii));
      [prof,bin] = meshgrid(1:it,1:ib);
      l.ru = full(sparse(bin(ii),prof(ii),l.uoce(:,1)));
      l.rv = full(sparse(bin(ii),prof(ii),l.uoce(:,2)));
      l.ru(ib,it) = NaN;
      l.rv(ib,it) = NaN;
      ii = find(~(z>=min(dr.z) & z<=max(dr.z)));
      l.ru(ii) = NaN;
      l.rv(ii) = NaN;

      uu = meanmediannan(d.ru(d.izu,:)+1i*d.rv(d.izu,:)...
             -l.ru(d.izu,:)-1i*l.rv(d.izu,:)+d.weight(d.izu,:)*0,2);
      ud = meanmediannan(d.ru(d.izd,:)+1i*d.rv(d.izd,:)...
             -l.ru(d.izd,:)-1i*l.rv(d.izd,:)+d.weight(d.izd,:)*0,2);
      clear l

    else
      iz = p.trusted_i;
      good = find(iz<=length(d.izu));
      uu = meanmediannan(d.ru(d.izu(iz(good)),:)+...
        1i*d.rv(d.izu(iz(good)),:)+...
        d.weight(d.izu(iz(good)),:)*0,1);
      good = find(iz<=length(d.izd));
      ud = meanmediannan(d.ru(d.izd(iz(good)),:)+...
        1i*d.rv(d.izd(iz(good)),:)+...
        d.weight(d.izd(iz(good)),:)*0,1);
    end
    %  try to take speed into account
    hrotvel = angle(uu./ud)*180/pi;
    d.rot_vel = hrotvel;
  %  keyboard

    %  decide which to use 
    if p.rotup2down==1
      hrot = hrotcomp;
      disp('    Rot up2down is using average up/down compass')
    elseif p.rotup2down==2
      hrot = hrotvel;
      disp('>   Rot up2down is using velocties')
    elseif p.rotup2down==3
      hrot = hrotcomp;
      disp('>   Rot up2down is using down compass only')
    elseif p.rotup2down==4
      hrot = hrotcomp;
      disp('>   Rot up2down is using up compass only')
    else
      hrot = d.ru(1,:)*0;
      disp('>   Not sure what you want. Applying no compass rotation at all.')
    end  

    % plot what I am about to do
    figure(2)
    clf
    orient tall
    plot([hrotcomp-90;hrotvel-180;hrot]')
    hold on
    ax = axis;
    plot(ax(1:2),[0 0],'-k',ax(1:2),[0 0]-90,'-k',ax(1:2),[0 0]-180,'-k') 
    text(1,20,'rotation used ')
    text(1,-60,'rotation compass')
    text(1,-150,'rotation velocity')
    ax(3) = -250;
    axis(ax)
    ylabel(' heading correction for uplooking')
    streamer([p.name,' Figure 5']);
    hgsave('tmp/5')

    if isfield(d,'hrot')
      disp('    Rotated earlier, using only difference ')
      oldrot = d.hrot;
    else
      oldrot = 0;
    end

    % apply rotation
    hrotm = meshgrid(hrot-oldrot,d.izd);
    % dont rotate NaNs
    ii = find(~isfinite(hrotm));
    hrotm(ii) = nmean(hrot-oldrot);
    % rotate both by half difference
    if p.rotup2down==4
      [ru,rv] = uvrot(d.ru(d.izd,:),d.rv(d.izd,:),-hrotm);
      d.ru(d.izd,:) = ru;
      d.rv(d.izd,:) = rv;
      [bu,bv] = uvrot(d.bvel(:,1),d.bvel(:,2),-(hrot-oldrot)');
      d.bvel(:,1) = bu;
      d.bvel(:,2) = bv;
    elseif p.rotup2down~=3
      [ru,rv] = uvrot(d.ru(d.izd,:),d.rv(d.izd,:),-hrotm/2);
      d.ru(d.izd,:) = ru;
      d.rv(d.izd,:) = rv;
      [bu,bv] = uvrot(d.bvel(:,1),d.bvel(:,2),-(hrot-oldrot)'/2);
      d.bvel(:,1) = bu;
      d.bvel(:,2) = bv;
    end

    hrotm = meshgrid(hrot-oldrot,d.izu);
    figure(4)
    clf
    subplot(2,1,1)
    plot(hrot)
    subplot(2,1,2)
    plot(oldrot)
    % dont rotate NaNs
    ii = find(~isfinite(hrotm));
    hrotm(ii) = nmean(hrot);
    % rotate both by half difference
    if p.rotup2down==4
      [ru,rv] = uvrot(d.ru(d.izu,:),d.rv(d.izu,:),hrotm);
      d.ru(d.izu,:) = ru;
      d.rv(d.izu,:) = rv;
    elseif p.rotup2down~=3
      [ru,rv] = uvrot(d.ru(d.izu,:),d.rv(d.izu,:),hrotm/2);
      figure(5)
      plot(hrotm(1,:)/2)
      d.ru(d.izu,:) = ru;
      d.rv(d.izu,:) = rv;
    end

    d.hrot = hrot;

  end
end
end
%%%%%%%%%%%% rotation

%--------------------------------------------------
function  hdg_offset = compoff(u1,u2)
% compute mean compass offset 

h1 = -angle(u1)*180/pi;
h1 = h1+(h1<0)*360;
nhead = 36;
dhead = 360/2/nhead;
h0 = linspace(5,355,nhead);

for n=1:nhead
  ind = find(abs(h1-h0(n))<=dhead);
  if length(ind)>1
    u1a(n) = mean(u1(ind));
    u2a(n) = mean(u2(ind));
  else
    u1a(n) = NaN;
    u2a(n) = NaN;
  end
end

ii = find(isfinite(u1a+u2a));
if length(ii)>0
  hdg_offset = angle(u1a(ii)/u2a(ii))*180/pi;
else
  hdg_offset = 0;
end

return


%-----------------------------------------------------
function y = boxav(x,nav)
% boxaverage and subsample

[ld,lv] = size(x);
in = fix(ld/nav);

for n=1:lv
  xm = reshape(x(1:(in*nav),n),nav,in);
  y(:,n) = nmean(xm)';
end
