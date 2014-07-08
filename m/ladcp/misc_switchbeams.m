function  [d,p,values] = misc_switchbeams(d,p,values)
% function  [d,p,values] = misc_switchbeams(d,p,values)
% 
% simple routine to "correct" if beams are switched
% it only performs this in 2D which is only roughly correct
% 
% p.beam_switch=[ 1 1 1]  switch u and v beam pairs down looker
% p.beam_switch=[ 1 0 2]  switch u not v beam pairs up looker
% 
%
disp('SWITCHBEAMS: ******** switch beams ***********')
if isfield(p,'beam_switch')==0
 disp(' need to give switch matix p.beam_switch ')
 return
end

% first rotate to instrument coordinates
%
% heading of beam 3 is
if ~isfinite(values.magdev)
  values.magdev = 0;
end
hrot = values.magdev + d.hdg(p.beam_switch(3),:);

% check if called before
if isfield(d,'ru_orig')==0
  d.ru_orig = d.ru;
  d.rv_orig = d.rv;
  d.bvel_orig = d.bvel;
else
  disp(' found ru_orig, will use those to rotate from ')
end

if p.beam_switch(3)==1
  iz = d.izd;
elseif  p.beam_switch(3)==2
  iz = d.izu;
else
  disp(' which instrument shall I switch? ')
  return
end

% backrotate to instrument coordinates
[rui,rvi] = uvrot(d.ru_orig(iz,:),d.rv_orig(iz,:),repmat(-hrot,[length(iz),1]));
if p.beam_switch(3)==1
  [bui,bvi] = uvrot(d.bvel_orig(:,1),d.bvel_orig(:,2),-hrot');
end


% switch beams
if p.beam_switch(1)
  rui = -rui;
  bui = -bui;
  warn = (' switched u-beam pair in instrument coordinates');
  disp(warn)
  p.warnp(size(p.warnp,1)+1,1:length(warn)) = warn;
end

if p.beam_switch(2)
  rvi = -rvi;
  bvi = -bvi;
  warn = (' switched v-beam pair in instrument coordinates');
  disp(warn)
  p.warnp(size(p.warnp,1)+1,1:length(warn)) = warn;
end

% rotate to earth coordinates
[ru,rv] = uvrot(rui,rvi,repmat(hrot,[length(iz),1]));
if p.beam_switch(3)==1
  [bu,bv] = uvrot(bui,bvi,hrot');
end

d.ru(iz,:) = ru;
d.rv(iz,:) = rv;
if p.beam_switch(3)==1
  d.bvel(:,1) = bu;
  d.bvel(:,2) = bv;
end
