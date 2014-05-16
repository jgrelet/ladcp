function plotraw(d,params,values)
% function plotraw(d,params,values)
%
% plot some results
%
% version 0.2	last change 21.09.2007

% M.Visbeck, G.Krahmann, IFM-GEOMAR

% wrong bins picked looking for broken beam	GK, Sep 2007	0.1-->0.2

%
% open and clear figure
%
figure(2)
clf
orient tall


%
% do not plot all points since that takes too long and
% creates too large postscript files
% thus reduce the used indices to maximum of 200
%
if length(d.time_jul)>200
  ind = fix(linspace(2,length(d.time_jul)-2,200));
else
  ind = [1:length(d.time_jul)];
end


%
% top subplot with vertical velocities for all bins
%
sp = subplot(411);
z = -d.izm(:,1)+d.z(1);
zz = [];
rw = [];
dz = [];


%
% find three beam solutions
% a profile is considered 3-beam if half of the good data are 3-beam
%
% they will be specially marked in this plot
%
iz = fliplr(d.izu);
n3bu = 0;
n3bd = 0;
if length(iz)>1
  dzu = nmedian(diff(d.zu));
  zz = [zz;-z(iz)];
  dz = [dz,iz*0+dzu];
  rw = [rw;d.rw(iz,ind)+d.weight(iz,ind)*0];
  % check for 3-beam solutions
  iw = sum(~isnan(d.rw(iz,ind)));
  ie = sum(~isnan(d.re(iz,ind)));
  i3bu = find(iw>(2*ie));
  n3bu = length(i3bu)/length(ind)*100;
end
iz = d.izd;
if length(iz)>1
  dzd = nmedian(diff(d.zd));
  zz = [zz;-z(iz)];
  dz = [dz,iz*0+dzd];
  rw = [rw;d.rw(iz,ind)+d.weight(iz,ind)*0];
  % check for 3-beam solutions
  iw = sum(~isnan(d.rw(iz,ind)));
  ie = sum(~isnan(d.re(iz,ind)));
  i3bd = find(iw>(2*ie));
  n3bd = length(i3bd)/length(ind)*100;
end


%
% contour results
%
col = jet(64);
col = ([[1 1 1]; col]);
colormap(col)
gcolor(ind,zz,nmedian(diff(ind))*ones(1,length(ind)),dz,rw);
hold on
ax = axis;
plot(ax(1:2),[0 0],'-k')
axis(ax) 
colorbar('horiz','peer',sp)


%
% mark 3-beam solutions
%
if n3bu>10
  l3b = NaN*ind;
  l3b(i3bu) = ax(3);
  plot(ind,l3b,'-r','linewidth',8)
  text(mean(ax(1:2)),ax(4),['found ',int2str(n3bu),...
	'% profiles 3 beam solutions'],...
	'VerticalAlignment','top','HorizontalAlignment','center','fontsize',14)
end
if n3bd>10
  l3b = NaN*ind;
  l3b(i3bd) = ax(4);
  plot(ind,l3b,'-r','linewidth',8)
  text(mean(ax(1:2)),ax(3),['found ',int2str(n3bd),...
	'% profiles 3 beam solutions'],...
 	'VerticalAlignment','bottom','HorizontalAlignment','center',...
	'fontsize',14)
end

ylabel('range [m]')
xlabel('ensemble')
if isfield(params,'name')
  streamer([params.name,' Figure 2']);
end
title(' W as function of bindepth and time')


%
% subplot showing beam performance (target strength/echo amplitude)
%
if isfield(d,'tsd_m')
  subplot(427)
  plot(d.tsd_m(1:length(d.izd),:),-d.zd)
  axis tight
  ax = axis;
  if isfield(d,'tsu_m')
    hold on
    plot(d.tsu_m(1:length(d.izu),:),d.zu)
    axis tight
    ax = axis;
    plot(ax(1:2),ax(1:2)*0,'-k')
  end

  t = d.tsd_m;
  checkbeam(t,ax,1)

  if isfield(d,'tsu_m')
    t = d.tsu_m;
    checkbeam(t,ax,0)
  end
  axis(ax)
  ylabel('distance [m]')
  xlabel('echo amplitude [dB]')
  title('Beam Performance')
end


%
% subplot showing range and correlation
%
if isfield(d,'cmd_m')
  subplot(428)
  plot(d.cmd_m(1:length(d.izd),:),-d.zd)
  axis tight
  ax = axis;
  [dum,dum,dum,x,y] = makebars(-d.zd,sum(isfinite(d.weight(d.izd,:))'));
  hold on
  fill(-y/max(y)*10,x,'r')

  if isfield(d,'cmu_m')
    hold on
    plot(d.cmu_m(1:length(d.izu),:),d.zu)
    plot(ax(1:2),ax(1:2)*0,'-k')
    [dum,dum,dum,x,y] = makebars(d.zu,sum(isfinite(d.weight(d.izu,:))'));
    fill(-y/max(y)*10,x,'g')
    axis tight
    ax = axis;
    for n=1:4
      text((0.12*n+0.27)*ax(2),ax(4),int2str(params.up_range(n)),...
	'VerticalAlignment','top')
    end
    if ~isnan(params.up_sn)
        text(0.01*ax(2),ax(4),['#',int2str(params.up_sn),...
            ' range:'],'VerticalAlignment','top')
    end
  end
  ax(1) = -13; 

  for n=1:4
    text((0.12*n+0.27)*ax(2),ax(3),int2str(params.dn_range(n)),...
	'VerticalAlignment','bottom')
  end
  if ~isnan(params.down_sn)
      text(0.01*ax(2),ax(3),['#',int2str(params.down_sn),...
         ' range:'],'VerticalAlignment','bottom')
  end
  axis(ax)
  ylabel('distance [m]')
  xlabel('correlation ')
  title('Range of good data')
  
end


%
% subplot showing depth of the package
%
if isfield(d,'z')
  subplot(813)
  plot(d.z)
  ylabel('depth')
  ax = axis;
  ax(3) = min(d.z*1.05);
  ax(4) = 0;
  ax(2) = length(d.z);
  axis(ax)
end


%
% subplot showing the tilt of the package
%
if isfield(d,'tilt')
  subplot(814)
  plot(d.tilt(1,:))
  hold on
  ax = axis;
  patch(ax([1,2,2,1]),[20,20,30,30],[1,0.75,0.75],'edgecolor','none')
  patch(ax([1,2,2,1]),[30,30,40,40],[1,0.5,0.5],'edgecolor','none')
  plot(d.tilt(1,:))
  ylabel('tilt [deg]')
  ax(2) = length(d.z);
  ax(4) = 30;
  axis(ax)
  set(gca,'yaxislocation','right')
end


%
% subplot showing the heading of the package
%
if isfield(d,'hdg')
  subplot(815)
  plot(d.hdg(1,:))
  ylabel('heading [deg]')
  ax = axis;
  ax(4) = 360;
  ax(2) =length(d.z);
  axis(ax )
  set(gca,'YTick',[0 90 180 270 360])
end


%
% subplot showing the transmit voltage of the instruments
%
if isfield(d,'xmv')
  subplot(816)
  plot(d.xmv(1,:))
  text(length(d.xmv)/2,mean(d.xmv(1,:)),[' mean: ',...
	num2str(fix(values.xmv(1)*10)/10)])
  ylabel('X-mit volt down')
  xlabel('ensemble')
  ax = axis;
  ax(2) = length(d.z);
  axis(ax)
  set(gca,'yaxislocation','right')
end

hgsave('tmp/2')






%=============================================
function checkbeam(t,ax,do)
% check beam performance
bl = size(t,1);			% this seems to have been a bug
				% the 1 was a 2 and thus picked
				% the wrong dimension
iend = fix(bl/2):bl;

tax = mean(ax(1:2));
if do
  tay = ax(3);
  tflag = 'bottom';
else
  tay = ax(4);
  tflag = 'top';
end

for i=1:4;
  % first correct for source level
  t(:,i) = t(:,i)-mean(t(:,i));
  % find noise level
  tn(i) = meanmediannan(t(iend,i),2);
  s2n(i) = mean(t(1:2,i))-tn(i);
  if s2n(i)==0
    s2n(i) = 1e100;
  end
end

ifail = s2n<max(s2n)*0.5;
ibad = ~ifail & s2n<max(s2n)*0.65;
iweak = ~ifail & ~ibad & s2n<max(s2n)*0.8;
for i=1:4;
  text(ax(1)+0.2*i*diff(ax(1:2)),tay,[int2str(s2n(i)./max(s2n)*100),'%'],...
 	'VerticalAlignment',tflag)
end

if sum(ifail)>0
  it = find(ifail==1);
  text(tax,tay*0.5,[' beam ',int2str(it),' broken'],...
 	'VerticalAlignment',tflag,'HorizontalAlignment','center','fontsize',18)
end

if sum(ibad)>0
  it = find(ibad==1);
  text(tax,tay*0.65,[' beam ',int2str(it),' bad'],...
 	'VerticalAlignment',tflag,'HorizontalAlignment','center','fontsize',15)
end

if sum(iweak)>0
  it = find(iweak==1);
  text(tax,tay*0.8,[' beam ',int2str(it),' weak'],...
 	'VerticalAlignment',tflag,'HorizontalAlignment','center','fontsize',12)
end
