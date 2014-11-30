function plot_result(dr,d,p,ps,values)
% function plot_result(dr,d,p,ps,values)
%                              
% - plot final velocity profile
%
% version 0.6	last change 13.07.2012

% Martin Visbeck and Gerd Krahmann, LDEO, April-2000

% TS axis bug fix, minor stuff            GK, Sep 2007    0.1-->0.2
% Bottom track axis,                      MV, Jul 2008    0.2-->0.3
% replaced finite with isfinite           GK, 08.11.2009  0.3-->0.4
% variable btrk plot range                GK, 12.05.2011  0.4-->0.5
% display instrument serial number        GK, 13.07.2012  0.5-->0.6

if isfield(dr,'range_do');
  zpmax = values.maxdepth + nmax([0;dr.range_do]);
else
  zpmax = values.maxdepth + nmax([0;dr.range]);
end

if isfield(p,'zbottom')
  if isfinite(p.zbottom)
    zpmax = p.zbottom;
  end
end

if isnan(p.plot_range(3))
  p.plot_range(3) = -zpmax;
end


ua = dr.u;
va = dr.v;
z = dr.z;


% 
% set range for the velocity plots
%
if isnan(p.plot_range(1))
  p.plot_range(1:2) = [-1,1]*abs(nmax([nmax(abs(ua)), nmax(abs(va))]))*120;
end


%
% plot final absolute velocity
%
sfigure(2);
clf
axes('position',[0.1 0.23 0.4 0.7])
plot(ua*100,-z,'-r','linewidth',2.5)
grid
hold on
plot(va*100,-z,'--g','linewidth',2.5)
ct = [' U(-) V(--)'];
if isfield(dr,'uerr')
  plot([va-dr.uerr]*100,-z,':g','linewidth',1.8)
  plot([va+dr.uerr]*100,-z,':g','linewidth',1.8)
  plot([ua-dr.uerr]*100,-z,':r','linewidth',1.8)
  plot([ua+dr.uerr]*100,-z,':r','linewidth',1.8)
end
if isfield(dr,'u_do')
  plot((dr.u_do+dr.ubar)*100,-z,'-r','linewidth',0.5)
  plot((dr.u_up+dr.ubar)*100,-z,'-r','linewidth',0.5)
  plot((dr.v_do+dr.vbar)*100,-z,'--g','linewidth',0.5)
  plot((dr.v_up+dr.vbar)*100,-z,'--g','linewidth',0.5)
  iz =[3:5:length(dr.u_do)];
  plot((dr.u_do(iz)+dr.ubar)*100,-z(iz),'.b','markersize',6)
  plot((dr.v_do(iz)+dr.vbar)*100,-z(iz),'.b','markersize',6)
  ct = [ct,'; blue dots down cast'];
  if nstd(dr.u_do)+nstd(dr.v_do)==0
    keyboard
  end
end
if isfield(dr,'u_shear_method')
  plot((dr.u_shear_method+dr.ubar)*100,-dr.z,'-r','linewidth',0.9)
  plot((dr.v_shear_method+dr.vbar)*100,-dr.z,'--g','linewidth',0.9)
  iz = [2:4:length(dr.u_shear_method)];
  plot((dr.u_shear_method(iz)+dr.ubar)*100,-dr.z(iz),'.r','markersize',7)
  plot((dr.v_shear_method(iz)+dr.vbar)*100,-dr.z(iz),'.g','markersize',7)
  ct=[ct,'; dotted shear'];
end
if isfield(dr,'u_sadcp')
  plot(dr.u_sadcp*100,-dr.z_sadcp,'-r','linewidth',1.3)
  plot(dr.v_sadcp*100,-dr.z_sadcp,'-g','linewidth',1.3)
  plot(dr.u_sadcp*100,-dr.z_sadcp,'pr','linewidth',0.9)
  plot(dr.v_sadcp*100,-dr.z_sadcp,'pg','linewidth',0.9)
  ct = [ct,'; pentagon SADCP'];
end
if (p.btrk_used>0 & isfield(dr,'zbot') )
  plot(dr.ubot*100,-dr.zbot,'r-','linewidth',1.3)
  plot(dr.vbot*100,-dr.zbot,'g--','linewidth',1.3)
  plot(dr.ubot*100,-dr.zbot,'r^','linewidth',1.3)
  plot(dr.vbot*100,-dr.zbot,'g^','linewidth',1.3)
end


% 
% set the axis limits for the plots
%
ax = p.plot_range;
if isfield(dr,'onlyshear')
  if dr.onlyshear
    text(ax(1),ax(3)*0.95,'  SHEAR SOLUTION ONLY','color','b','fontsize',15)
  end
end
if ps.up_dn_looker==2
  text(ax(1),ax(3)*0.90,'  DOWN LOOKER ONLY','color','b','fontsize',13)
end
if ps.up_dn_looker==3
  text(ax(1),ax(3)*0.90,'    UP LOOKER ONLY','color','b','fontsize',13)
end
plot([0 0],ax(3:4),'-k')
axis(ax)

title(ct)
ylabel('depth [m]')

% plot bottom track data

if (p.btrk_used>0 & isfield(dr,'zbot') )
  plot(p.plot_range(1:2),[1 1]*abs(p.plot_range(3)),'-k','linewidth',2)
  set(gca,'XTickLabel',[]);
  zbprofr=-dr.zbot+p.zbottom;
  axes('position',[0.1 0.1 0.4 0.12])
  iz = find(zbprofr>10 & zbprofr <p.btrk_plot_range);
  plot(dr.ubot(iz)*100,zbprofr(iz),'r-','linewidth',2.5)
  grid
  hold on
  plot(dr.vbot(iz)*100,zbprofr(iz),'g--','linewidth',2.5)
  if isfield(dr,'uerrbot')
    plot([dr.ubot(iz)+dr.uerrbot(iz)]*100,zbprofr(iz),':r','linewidth',1.8)
    plot([dr.ubot(iz)-dr.uerrbot(iz)]*100,zbprofr(iz),':r','linewidth',1.8)
    plot([dr.vbot(iz)+dr.uerrbot(iz)]*100,zbprofr(iz),':g','linewidth',1.8)
    plot([dr.vbot(iz)-dr.uerrbot(iz)]*100,zbprofr(iz),':g','linewidth',1.8)
  end
  plot([0 0],[0 -nmax(-[p.btrk_plot_range, p.zbottom])],'-k')
  xlabel('velocity [cm/s]')
  ylabel('above bottom [m]')
  ax2 = [p.plot_range(1:2) 0 min([p.btrk_plot_range, p.zbottom])];
  axis(ax2);
  if p.btrk_used>1 
    text(0.9*ax2(1),ax2(4)*0.92,'post processed bottom track')
  else
    text(0.9*ax2(1),ax2(4)*0.92,'RDI bottom track')
  end
end


%
% write some text information
%
subplot(322)

iy =1.2;
ix = - 0.1;
idy = 1/8;

iy = iy-idy;
[slat,slon] = pos2str(values.start_pos);
text(ix,iy,[' Start:'])
text(ix+0.2,iy,[slat])
text(ix+0.7,iy,[slon])

date1 = gregoria(values.start_time);  
ds = datenum(date1(1),date1(2),date1(3),date1(4),date1(5),date1(6));
iy = iy-idy;
text(ix+0.2,iy,[datestr(ds,0)])

[slat,slon] = pos2str(values.end_pos);
iy = iy-idy;
text(ix,iy,[' End:'])
text(ix+0.2,iy,[slat])
text(ix+0.7,iy,[slon])

date2 = gregoria(values.end_time);    
ds = datenum(date2(1),date2(2),date2(3),date2(4),date2(5),date2(6));
iy = iy-idy;
text(ix+0.2,iy,[datestr(ds,0)])

iy = iy-idy;
text(ix,iy,sprintf('u-mean: %3.0f cm/s    v-mean %3.0f cm/s',...
    nmean(ua)*100,nmean(va)*100))

iy = iy-idy;
if length(d)>0
 if length(d.zu>2), 
   zu = d.zu; 
  else, 
   zu = [0 0]; 
  end
  text(ix,iy,['binsize do: ',num2str(diff(d.zd([1,2]))),...
	' m  binsize up:  ', num2str(diff(zu([1,2]))),' m'])
end

if isfield(values,'inst_serial')
  if length(values.inst_serial)==1
    tstr = ['S/N down : ',int2str(values.inst_serial)];
  elseif length(values.inst_serial)==2
    tstr = ['S/N down : ',int2str(values.inst_serial(1)),...
        '     S/N up : ',int2str(values.inst_serial(2))];
  else
    tstr = '';
  end
  iy = iy-idy;
  text(ix,iy,tstr);
end

idy = 1/10;
iy = iy-idy;
te = ['mag. deviation ',num3str(values.magdev,4,1),'^o'];
text(ix,iy,te)

iy = iy-idy;
text(ix,iy,['wdiff: ',num2str(p.wlim),...
'  pglim: ',num2str(p.pglim),'  elim ',num2str(p.elim)])

iy = iy-idy;
dum=' ';
if ps.smoofac~=0
  dum=[dum,'smo:',num3str(ps.smoofac,3,2)];
end
if ps.dragfac~=0
  dum=[dum,' dra:',num3str(ps.dragfac,3,2)]; 
end
if sum(ps.smallfac(:,2))>0
  dum=[dum,' smal:',int2str(ps.smallfac(1,1)),'-',int2str(ps.smallfac(end,1))]; 
end
if ps.barofac~=0
  dum=[dum,' bar:',num3str(ps.barofac,3,1)];
end
if ps.botfac~=0
  dum=[dum,' bot:',num3str(ps.botfac,3,1)];
end
if ps.sadcpfac~=0
  dum=[dum,' sad:',num3str(ps.sadcpfac,3,1)]; 
end
text(ix,iy,dum)

iy = iy-idy;
text(ix,iy,['weightmin ',num3str(ps.weightmin,3,1),...
	'  weightpower: ',num3str(ps.weightpower,3,1)])


iy = iy-idy;
ctext = ['max depth: ',int2str(values.maxdepth),' m'];
if p.zbottom>-Inf, 
  ctext=[ctext,'   bottom: ',int2str(p.zbottom),' m']; 
end
text(ix,iy,ctext)

axis off


%
% plot profiles of target strength
%
axes('position',[0.58 0.39 0.1 .25])
plot(dr.ts,-dr.z/1000,'b-','linewidth',1.5)
hold on
plot(dr.ts_out,(-max(d.zd)-dr.z)/1000,'k-')
axi = axis;
if isfinite(nmax(dr.ts))
  ax(1) = -0.9*nmax(-dr.ts_out); 
  ax(2) = nmax(dr.ts)*1.1; 
end
ax(3:4) = ax(3:4)/1000;
if isnan(ax(1))
  ax(1) = axi(1);
end
axis(ax)
ylabel('depth [km]')
xlabel('target strength [dB]','color','b')
set(gca,'fontsize',10);


%
% plot profiles of super ensemble data range
%
axes('position',[0.70 0.39 0.1 .25])

plot(sum(dr.range,2),-dr.z/1000,'r-','linewidth',1.5)
if isfield(dr,'range_up')
  hold on
  plot(dr.range_up,-dr.z/1000,'-k')
  plot(dr.range_do,-dr.z/1000,'-k')
  iz = round(linspace(2,length(dr.z)-2,10));
  plot(dr.range_do(iz),-dr.z(iz)/1000,'.b')
end
if isfinite(nmax(sum(dr.range,2)))
  ax(1)=0; 
  ax(2)=1.1*nmax(sum(dr.range,2)); 
end
axis(ax)
title('range [m]','color','r')
set(gca,'YtickLabel',[])
set(gca,'fontsize',10)


if isfield(dr,'uerr')==1

  % plot profiles of velocity error
  axes('position',[0.82 0.39 0.1 .25])

  ue = dr.uerr;
  plot(ue,-dr.z/1000,'k-','linewidth',1.5)
  hold on
  ax(1) = 0; 
  ax(2) = min([0.001,3*nmedian(ue)]);
  axis(ax)
  xlabel('vel error (-k) [m/s]')
  set(gca,'fontsize',10)
  set(gca,'YtickLabel',[])

end


%
% add single ping accuracy
%
if isfield(dr,'ensemble_vel_err')==1
  plot(dr.ensemble_vel_err,-dr.z/1000,'-b')
  ax(2) = max([2.5*nmedian(ue),2.5*nmedian(dr.ensemble_vel_err)]);
  axis(ax)
  title('single ping (-b)','color','b')
end


%
% plot position of CTD
%
subplot(326)
xctd = dr.xctd;
yctd = dr.yctd;
ii = fix(linspace(1,length(xctd),10));
[m,ib] = min(dr.zctd);
plot(xctd,yctd,'linewidth',2)
hold on
plot(xctd(ii),yctd(ii),'r.','markersize',10)
plot(xctd(ib),yctd(ib),'g+','markersize',9)
if isfield(dr,'xship')
  plot(dr.xship,dr.yship,'-g',dr.xship(ii),dr.yship(ii),'k.','markersize',10)
  plot([xctd(ii);dr.xship(ii)],[yctd(ii); dr.yship(ii)],'-y','linewidth',0.5)
  xlabel('CTD-position (blue) and ship (green) east-west [m]')
else
  xlabel('CTD-position east-west [m]')
end
uship_a = values.uship+sqrt(-1)*values.vship;
if ~(abs(uship_a)==0 & values.lat==0 & values.lon==0)
  plot(values.xdisp,values.ydisp,'gp')
  text(values.xdisp,values.ydisp,' GPS-end')
end
text(xctd(ib),yctd(ib),' bottom')
text(xctd(end),yctd(end),' end')
axis equal
text(xctd(1),yctd(1),' start')
ylabel('north-south [m]')
grid
set(gca,'fontsize',10)

suplabel(['Station : ',p.name,'  Figure 1'],'t');
suplabel(p.software,'x');
axis off


orient tall
hgsave('tmp/1');
