function p=checkbtrk(d,di,de,dr,p);
% function p=checkbtrk(d,di,de,dr,p);

%
% general function info
%
disp(' ')
disp(['CHECKBTRK:  check bottom track against U_ctd solution '])


%
% check whether any bottom was found or used at all
%
if ~isfinite(p.zbottom) 
  return 
end


%
% store the data in different variables
%
% now, why is that necessary   GK
%
if ~isfield(d,'bvel_own') 
  d.bvel_own = d.bvel*nan; 
  d.hbot_own = d.hbot*nan; 
end
if ~isfield(d,'bvel_rdi')
  d.bvel_rdi = d.bvel*nan; 
  d.hbot_rdi = d.hbot*nan; 
end


%
% extract the indices of ensembles where the bottom is close
%
ii = find(abs(d.z+p.zbottom)<max(p.btrk_range) & ...
          abs(d.z+p.zbottom)>min(p.btrk_range));
 
if isfield(d,'bvel_rdi')
  bu_rdi = d.bvel_rdi(ii,1);
  bv_rdi = d.bvel_rdi(ii,2);
  bw_rdi = d.bvel_rdi(ii,3);
else
  bu_rdi = d.bvel(ii,1)*nan;
  bv_rdi = d.bvel(ii,2)*nan;
  bw_rdi = d.bvel(ii,3)*nan;
end

if isfield(d,'hbot_rdi')
  bh_rdi = d.hbot_rdi(ii);
else
  bh_rdi = d.hbot(ii)*nan;
end

bu_own = d.bvel_own(ii,1);
bv_own = d.bvel_own(ii,2);
bw_own = d.bvel_own(ii,3);
bh_own = d.hbot_own(ii);


%u_o = interp1(dr.z,dr.u,-d.z(ii)');
%v_o = interp1(dr.z,dr.v,-d.z(ii)');
u_c = -interp1(dr.tim',dr.uctd',d.time_jul(ii)');
v_c = -interp1(dr.tim',dr.vctd',d.time_jul(ii)');
  
rw = d.rw(d.izd,ii);
ij = find((d.izm(d.izd,ii)+p.zbottom)<d.zd(1));
rw(ij) = nan;
if isfield(d,'wctd')
  w_c = d.wctd(ii)';
  disp('    Using CTD pressure derived W')
else
  w_c = meanmediannan(rw,3)';
  disp('    Using reference layer W')
end
w_c2 = w_c-[0;diff(w_c)/2 ];

bu_used = d.bvel(ii,1);
bv_used = d.bvel(ii,2);
bw_used = d.bvel(ii,3);
bh_used = d.hbot(ii);

%blen = d.down.Cell_length/100;
blen = d.blen;				% this seems to be better for NB
 
% extract raw velocity data near bottom
rube = d.ru(d.izd,ii);
rvbe = d.rv(d.izd,ii);
rwbe = d.rw(d.izd,ii);

% the same, but only velocities that passed editing
inan = find(~isfinite(d.weight));
ru = d.ru; 
ru(inan) = NaN; 
ru = ru(d.izd,ii);
rv = d.rv; 
rv(inan) = NaN; 
rv = rv(d.izd,ii);
rw = d.rw; 
rw(inan) = NaN; 
rw = rw(d.izd,ii);

% correct velocities with solved U_ctd
rua = ru-repmat(u_c',[length(d.izd) 1]); 
ruabe = rube-repmat(u_c',[length(d.izd) 1]);
rva = rv-repmat(v_c',[length(d.izd) 1]); 
rvabe = rvbe-repmat(v_c',[length(d.izd) 1]);
rwa = rw-repmat(w_c',[length(d.izd) 1]); 
rwabe = rwbe-repmat(w_c',[length(d.izd) 1]);
rwan = abs(rw)./abs(repmat(w_c',[length(d.izd) 1]));
rwanbe = abs(rwbe)./abs(repmat(w_c',[length(d.izd) 1]));



bua_rdi = bu_rdi-u_c;
bva_rdi = bv_rdi-v_c;
%
% try to use w inbetween pings for real bottom track
%
bwa1_rdi = bw_rdi-w_c;
bwa2_rdi = bw_rdi-w_c2;
if nstd(bwa1_rdi)>nstd(bwa2_rdi)
  bwa_rdi = bwa2_rdi;
  disp('    Using inbetween w-ref')
else
  bwa_rdi = bwa1_rdi;
end

bua_own = bu_own-u_c;
bva_own = bv_own-v_c;
bwa_own = bw_own-w_c;

bua_used = bu_used-u_c;
bva_used = bv_used-v_c;
bwa_used = bw_used-w_c;


%
% save bias and std of bottom track anomaly
%
p.btrk_u_bias = meanmediannan(bua_used,6);
p.btrk_u_std = nstd(bua_used);
p.btrk_v_bias = meanmediannan(bva_used,6);
p.btrk_v_std = nstd(bva_used);
p.btrk_w_bias = meanmediannan(bwa_used,6);
p.btrk_w_std = nstd(bwa_used);

disp(['    Profiles within give acceptable range: ',int2str(length(ii))])
disp(['    U bias :',num3str(p.btrk_u_bias,6,3),' [m/s]  std: ',...
       num3str(p.btrk_u_std,5,3),' [m/s]'])
disp(['    V bias :',num3str(p.btrk_v_bias,6,3),' [m/s]  std: ',...
       num3str(p.btrk_v_std,5,3),' [m/s]'])
disp(['    W bias :',num3str(p.btrk_w_bias,6,3),' [m/s]  std: ',...
       num3str(p.btrk_w_std,5,3),' [m/s]'])

if abs(p.btrk_u_bias)>0.1
  warn = ['>   Large U bottom track bias ',num2str(p.btrk_u_bias)];
  disp(warn)
  p.warnp(size(p.warnp,1)+1,1:length(warn)) = warn;
end

if abs(p.btrk_v_bias)>0.1
  warn = ['>   Large V bottom track bias ',num2str(p.btrk_v_bias)];
  disp(warn)
  p.warnp(size(p.warnp,1)+1,1:length(warn)) = warn;
end


%
% correct nominal depth with detected bottom depth
%
izm = d.izm(d.izd,ii);
zcor = (-d.hbot(ii)+d.z(ii));
zcor_std = nstd(zcor);
p.btrk_rough = zcor_std;
izm = izm-repmat(zcor,[length(d.izd) 1]);


%
% compute bias W slope
%
dz = abs(diff(d.zd(1:2)));
wbias = [];
zbias = [];
nbias = [];
for i=1:6
  dzi = dz*i/2;
  ij = find(abs(izm+dzi)<(dz/4));
  nok = sum(isfinite(rwanbe(ij)));
  if nok>5
    wbias = [wbias, 1-nmedian(rwanbe(ij))];
    zbias = [zbias, dzi];
    nbias = [nbias, nok];
  end
end
if length(wbias)>1
  p.wbslope = polyfit(zbias,wbias,1);
else
  wbias = nan;
  zbias = nan;
  p.wbslope = [nan nan];
end

disp(['    W slope fact :',num3str(p.wbslope(1),6,4),...
	' [1/m] lower W below bottom '])


%
% compute 'wave' action
%
p.btrk_wdiff = nmedian(abs(diff(w_c)));
disp(['    W diff :',num3str(p.btrk_wdiff,6,4),...
	' [m/s] ping to ping w difference '])



disp(['    H std :',num3str(zcor_std,6,1),...
	' [m]  large means bottom is rough/sloped'])

p.btrk_tilt_mean = mean(d.tilt(ii));
p.btrk_tilt_std = std(d.tilt(ii));
disp(['    Tilt mean :',num3str(p.btrk_tilt_mean,3,1),' +/- ',...
           num3str(p.btrk_tilt_std,3,1),' [^o]  '])


% select bins used for own bottom track
iz = find(abs(izm+blen*p.btrk_below)<blen);


% plot results
tx2 = 0.25;
figure(2)
clf
  
subplot(221); hold on;
plot(ruabe,izm,'.k');			% all velocities
plot(rua,izm,'.g'); 			% edited velocities
plot(rua(iz),izm(iz),'.r');		% used-for-BT velocities
axis tight;
ax = [-0.8,0.8,-160,50];
vh = floor(ax(1)*20)/20+0.0125 : 0.025: ceil(ax(2)*20)/20;
plot(ax(1:2),[0 0],'-')
plot([0 0],ax(3:4),'k-','linewidth',1.5)
plot([0 0]+p.btrk_u_bias,ax(3:4),'r-','linewidth',1.5)
text(tx2,10,' Bottom')

if sum(isfinite(bua_rdi))>0
  [y,x] = hist(bua_rdi,vh);
  y([1,end]) = 0;
  ys = sum(y);
  [x,y] = stairs(x,y);
  fill(x,y/max(y)*30-50,'-g')
  plot([0 0]+meanmediannan(bua_rdi,6),ax(3:4),'-g')
  text(tx2,-40,['RDI n: ',int2str(ys)])
  text(0.1,-10,[' bias ',num3str(meanmediannan(bua_rdi,6),6,3)],'color','g')
end

[y,x] = hist(bua_own,vh);
y([1,end]) = 0;
ys = sum(y);
[x,y] = stairs(x,y);
fill(x,y/max(y)*30-100,'-b')
plot([0 0]+meanmediannan(bua_own,6),ax(3:4),'-b')
text(tx2,-90,[' own n: ',int2str(ys)])
text(0.1,-60,[' bias ',num3str(meanmediannan(bua_own,6),6,3)],'color','b')

[y,x] = hist(real(de.bvel-de.uctd(:,1).'),vh);
y([1,end]) = 0;
ys = sum(y);
[x,y] = stairs(x,y);
fill(x,y/max(y)*30-150,'-r')
text(tx2,-140,['Sup. ens. n: ',int2str(ys)])
text(0.1,-110,[' bias ',num3str(p.btrk_u_bias,6,3)],'color','r')
title([' U bot-tr mean U_{ctd} : ',num3str(-nmean(real(de.uctd(:,1))),4,2),...
        ' [m/s] added'])
xlabel(' U [m/s]')
ylabel('depth [m]')
axis(ax)
text(-0.7,-154,['Bottom roughness: ',int2str(p.btrk_rough),' [m]'])

subplot(222); hold on;
plot(rvabe,izm,'.k'); 			% all velocities
plot(rva,izm,'.g'); 			% edited velocities
plot(rva(iz),izm(iz),'.r')		% used-for-BT velocities
axis tight
ax = [-0.8,0.8,-160,50];
plot(ax(1:2),[0 0],'-')
plot([0 0],ax(3:4),'k-')
plot([0 0]+p.btrk_v_bias,ax(3:4),'r-','linewidth',1.5)

text(tx2,10,' Bottom')
if sum(isfinite(bva_rdi))>0
  [y,x] = hist(bva_rdi,vh);
  y([1,end]) = 0;
  ys = sum(y);
  [x,y] = stairs(x,y);
  fill(x,y/max(y)*30-50,'-g')
  plot([0 0]+meanmediannan(bva_rdi,6),ax(3:4),'-g')
  text(tx2,-40,[' RDI n: ',int2str(ys)])
  text(0.1,-10,[' bias ',num3str(meanmediannan(bva_rdi,6),6,3)],'color','g')
end

[y,x] = hist(bva_own,vh);
y([1,end]) = 0;
ys = sum(y);
[x,y] = stairs(x,y);
fill(x,y/max(y)*30-100,'-b')
plot([0 0]+meanmediannan(bva_own,6),ax(3:4),'-b')
text(tx2,-90,[' own n: ',int2str(ys)])
text(0.1,-60,[' bias ',num3str(meanmediannan(bva_own,6),6,3)],'color','b')

[y,x] = hist(imag(de.bvel-de.uctd(:,1).'),vh);
y([1,end]) = 0;
ys = sum(y);
[x,y] = stairs(x,y);
fill(x,y/max(y)*30-150,'-r')
text(tx2,-140,['Sup. ens. n: ',int2str(ys)])
text(0.1,-110,[' bias ',num3str(p.btrk_v_bias,6,3)],'color','r')
title([' V bot-tr mean V_{ctd} : ',num3str(-nmean(imag(de.uctd(:,1))),4,2),...
        ' [m/s] added'])
xlabel(' V [m/s]')
ylabel('depth [m]')
text(-0.7,-154,['Tilt : ',num3str(p.btrk_tilt_mean,3,1),' +/- ',...
                 num3str(p.btrk_tilt_std,3,1),' [^o]'])
axis(ax)

tx2 = 0.05;

subplot(223); hold on;
plot(rwabe,izm,'.k');				% all velocities
plot(rwa,izm,'.g');				% edited velocities
plot(rwa(iz),izm(iz),'.r');			% used-for-BT velocities
axis tight
ax=axis;
ax(1:2)=[-1 1]*0.4;
vh=floor(ax(1)*20)/20+0.005 : 0.01: ceil(ax(2)*20)/20;
ax(3)=-160;
ax(4)=50;
plot(ax(1:2),[0 0],'-')
plot([0 0],ax(3:4),'k-')
plot([0 0]+p.btrk_w_bias,ax(3:4),'r-','linewidth',1.5)

text(tx2,10,' Bottom')
if sum(isfinite(bwa_rdi))>0
  [y,x] = hist(bwa_rdi,vh);
  y([1,end]) = 0;
  ys = sum(y);
  [x,y]=stairs(x,y);
  fill(x,y/max(y)*30-50,'-g')
  plot([0 0]+meanmediannan(bwa_rdi,6),ax(3:4),'-g')
  text(tx2,-40,[' RDI n: ',int2str(ys)])
  text(0.1,-10,[' bias ',num3str(meanmediannan(bva_rdi,6),6,3)],'color','g')
end

[y,x] = hist(bwa_own,vh);
y([1,end]) = 0;
ys = sum(y);
[x,y] = stairs(x,y);
fill(x,y/max(y)*30-100,'-b')
plot([0 0]+meanmediannan(bwa_own,6),ax(3:4),'-b')
text(tx2,-90,[' own n: ',int2str(ys)])
text(0.1,-60,[' bias ',num3str(meanmediannan(bwa_own,6),6,3)],'color','b')

[y,x] = hist(di.bvel(3,:)-nmedian(di.rw),vh);
y([1,end]) = 0;
ys = sum(y);
[x,y] = stairs(x,y);
fill(x,y/max(y)*30-150,'-r')
text(tx2,-140,['Sup. ens. n: ',int2str(ys)])
text(0.1,-110,[' bias ',num3str(p.btrk_w_bias,6,3)],'color','r')
text(-0.3,-154,['W diff: ',num3str(p.btrk_wdiff,6,3),' [m/s]'])

title([' W bottom track [m/s] '])
xlabel(' W [m/s]')
ylabel('depth [m]')
axis(ax)

subplot(224); hold on;
plot(rwanbe,izm,'.k');				% all velocities
plot(rwan,izm,'.g');				% edited velocities
plot(rwan(iz),izm(iz),'.r');			% used-for-BT velocities
axis tight
ax = [0.5,1.2,-60,60];
title([' W bottom track factor '])
xlabel(' abs(W)/abs(W_{ref}) ')
ylabel('depth [m]')
axis(ax)
plot(ax(1:2),[0 0],'-')
plot([1 1],ax(3:4),'k-','linewidth',1.5)
text(tx2,10,' Bottom')
plot(1-polyval(p.wbslope,[40 0]),[-40 0],'-b')
plot(1-wbias,-zbias,'pb')
grid
text(0.54,-50,[' bias slope ',num3str(p.wbslope(1),6,3),' [1/m]'])
text(0.54,-58,[' offset  ',num3str(p.wbslope(2),6,3) ])
text(0.54,50,[' beam angle  ',int2str(d.down.Beam_angle)])
text(0.54,42,[' bin length  ',int2str(d.down.Cell_length/100),' [m]'])
text(0.54,34,[' Pings/Ens  ',int2str(d.down.Pings_per_Ensemble)])
text(0.54,26,[' Frequency  ',int2str(d.down.Frequency),' [kHz]'])

text(0.54,10,[' Zbottom  ',int2str(p.zbottom),' [m]'])

 
%streamer([p.name,' Figure 13']);
suplabel([p.name,' Figure 13'],'t');

orient tall
hgsave('tmp/13');
