function [d,p]=uvwrot(d,p,iup)
% function [d,p]=uvwrot(d,p,iup)
% rotate velocities by p.drot p.dpit p.drol
% iup=1, up-looking data olny
% iup=0, down-looking data only

p.uwvrot = 1;
p = setdefv(p,'drot',0);
p = setdefv(p,'dpit',0);
p = setdefv(p,'drol',0);


% which bins to transform
if nargin<3,
  iz = [1:size(d.ru,1)];
  it = 1;
else
  if iup
    % up looking only
    iz=d.izu;
    it=2;
    disp(' up-looking data pitch-roll offset')
    disp([' pitch: ',num2str(p.dpit),' roll: ',num2str(p.dpit)])
  else
    % down looking only
    iz=d.izd;
    it=1;
    disp(' down-looking data pitch-roll offset')
   disp([' pitch: ',num2str(p.dpit),' roll: ',num2str(p.dpit)])
  end
end

d=setdefv(d,'drot',[0 0]);
d=setdefv(d,'dpit',[0 0]);
d=setdefv(d,'drol',[0 0]);

drot=p.drot-d.drot(it);
dpit=p.dpit-d.dpit(it);
drol=p.drol-d.drol(it);

d.drot(it)=p.drot;
d.dpit(it)=p.dpit;
d.drol(it)=p.drol;



% precompute some constants
d2r=pi/180;	% conversion from degrees to radians

[nb,ne]=size(d.ru);

%big loop
for ii=1:ne

 % get roll, pitch and heading
 roll=d.rol(it,ii);
 pitch=d.pit(it,ii);
 head=d.hdg(it,ii);

 roll2=roll+drol;
 pitch2=pitch+dpit;
 head2=head+drot;

 earth=[d.ru(:,ii),d.rv(:,ii),d.rw(:,ii)];

 % Step 1 - determine rotation angles from sensor readings
 % fixed sensor case
 % make sure everything is expressed in radians for MATLAB
 RR=roll.*d2r;
 KA=sqrt(1.0 - (sin(pitch.*d2r).*sin(roll.*d2r)).^2);
 PP=asin(sin(pitch.*d2r).*cos(roll.*d2r)./KA);
 HH=head.*d2r;

 RR2=roll2.*d2r;
 KA=sqrt(1.0 - (sin(pitch2.*d2r).*sin(roll2.*d2r)).^2);
 PP2=asin(sin(pitch2.*d2r).*cos(roll2.*d2r)./KA);
 HH2=head2.*d2r;

 % Step 2 - calculate trig functions and scaling factors
 CP=cos(PP); CR=cos(RR); 
 SP=sin(PP); SR=sin(RR); 

 CP2=cos(PP2); CR2=cos(RR2); 
 SP2=sin(PP2); SR2=sin(RR2); 

 CH=cos(HH); SH=sin(HH);

 CH2=cos(HH2); SH2=sin(HH2);

 % rotation matrix
 %	VX =  VXE.*(CH*CR + SH*SR*SP) + VYE.*SH.*CP + VZE.*(CH*SR - SH*CR*SP);
 %	VY = -VXE.*(SH*CR - CH*SR*SP) + VYE.*CH.*CP - VZE.*(SH*SR + CH*SP*CR);
 %  VZ = -VXE.*(SR*CP)            + VYE.*SP     + VZE.*(CP*CR);

	R(1,1) = (CH*CR + SH*SR*SP);
    R(1,2) = SH.*CP ;
    R(1,3) = (CH*SR - SH*CR*SP);
	R(2,1) = -(SH*CR - CH*SR*SP);
    R(2,2) = CH.*CP;
    R(2,3) = -(SH*SR + CH*SP*CR);
    R(3,1) = -(SR*CP);
    R(3,2) = SP;
    R(3,3) = (CP*CR);

    RI=inv(R);

	R2(1,1) = (CH2*CR2 + SH2*SR2*SP2);
    R2(1,2) = SH2.*CP2 ;
    R2(1,3) = (CH2*SR2 - SH2*CR2*SP2);
	R2(2,1) = -(SH2*CR2 - CH2*SR2*SP2);
    R2(2,2) = CH2.*CP2;
    R2(2,3) = -(SH2*SR2 + CH2*SP2*CR2);
    R2(3,1) = -(SR2*CP2);
    R2(3,2) = SP2;
    R2(3,3) = (CP2*CR2);


 for IB=iz,

    % Step 4: rotate to instument coordinates
    vb = RI*earth(IB,:)';

    % Step 4: rotate to earth coordinates
    ve = R2*vb;

    earth(IB,:)=ve';

 end	

 d.ru(iz,ii)=earth(iz,1);
 d.rv(iz,ii)=earth(iz,2);
 d.rw(iz,ii)=earth(iz,3);

end % Big Loop

