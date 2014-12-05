function [d,p]=misc_fixcompass(d,p)
% 
% fix compass issues
%
% p.fix_compass=1  apply offset only
% p.fix_compass=2  fix up looking compass to down looking compass + offset
% p.fix_compass=3  fix down looking compass to up looking compass + offset
%
% p.hdg_offset will apply a fixed rotation to the compass. In most cases you
%        will want to set the compass that you trust to 0 and rotate the other one
%        by the difference in direction between the two beam 3 headings as mounted 
%        on the frame
%        : [downoffset  upoffset]
%
% once you done this rotup2down gets disabled since I assume you know what the offset was
%
% version 0.2  last change 29.08.2011

% M. Visbeck LDEO 2004

% small bugs as this has for a long time not been used   GK, 29.08.2011  0.1-->0.2

disp(' ')
disp('FIXCOMPASS:  adjust compass')

% use tilt sensors to figure out allingment between instruments
if length(d.izu)>0
  diary off
  hoff2 = fminsearch('checktilt',0,[],[d.rol(2,:);d.pit(2,:);d.rol(1,:);d.pit(1,:)]);
  diary on
  disp(['  Mean instrument alignment based on tilt is ',num2str(hoff2),' deg'])
end

[lb,lt] = size(d.ru);
if p.fix_compass==2
  p = setdefv(p,'hdg_offset',[0 hoff2]);
elseif p.fix_compass==3
  p = setdefv(p,'hdg_offset',[-hoff2 0]);
else
  p = setdefv(p,'hdg_offset',[0 0]);
end
p = setdefv(p,'hdg_offset_used',[0 0]);

rot = p.hdg_offset;

rotv(1,:) = ones(1,lt)*rot(1);
disp(['  Offset down looking compass by : ',num2str(rot(1))]) 
p.hdg_offset_used(1) = p.hdg_offset_used(1)+rot(1);

if p.fix_compass==3
  disp('  Using up looking compass on down looking ADCP')
  rotv(1,:) = d.hdg(1,:)-d.hdg(2,:)+rot(1);
  p.rotup2down = 0;
end

if length(d.izu)>0
  rotv(2,:) = ones(1,lt)*rot(2);
  disp(['  Offset up looking compass by : ',num2str(rot(2))]) 
  p.hdg_offset_used(2) = p.hdg_offset_used(2)+rot(2);
  if p.fix_compass==2
    disp('  Using down looking compass on up looking ADCP')
    rotv(2,:) = d.hdg(2,:)-d.hdg(1,:)+rot(2);
    p.rotup2down = 0;
  end
end

rotm = meshgrid(rotv(1,:),d.izd);
[d.ru(d.izd,:),d.rv(d.izd,:)] = uvrot(d.ru(d.izd,:),d.rv(d.izd,:),rotm);

if length(rotv(1,:))==size(d.bvel,1)
  [d.bvel(:,1),d.bvel(:,2)] = uvrot(d.bvel(:,1),d.bvel(:,2),rotv(1,:)');
else
  [d.bvel(1,:),d.bvel(2,:)] = uvrot(d.bvel(1,:),d.bvel(2,:),rotv(1,:));
end

if length(d.izu)>0
  rotm = meshgrid(rotv(2,:),d.izu);
  [d.ru(d.izu,:),d.rv(d.izu,:)] = uvrot(d.ru(d.izu,:),d.rv(d.izu,:),rotm);
end
