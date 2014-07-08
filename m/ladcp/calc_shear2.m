function [ds,dr,ps,p,messages] = getshear(d,p,ps,dr,messages)
%
% - compute shear profiles 
% - use only central difference
% - use 2*std editing
% 
%  Martin Visbeck, LDEO, 3/7/97


%
% general function info
%
disp(' ')
disp(['CALC_SHEAR2:  calculate a velocity profile based on shears only'])


%
% resolution of final shear profile in meter
%
idz = ps.dz;
disp(['    Averaging shear profile over (ps.dz) ',num2str(idz),' [m]'])


%
% average over how much percent of the data
%
stdf = ps.shear_stdf;
disp(['    Maximum std (stdf) ',num2str(stdf),' of data '])
if nargin<4, 
  dr.dummy = 0; 
end

%
% check if only one istrument is to be used
%
if ps.up_dn_looker==2
  % down looker only
  d.weight(d.izu,:)=nan;
elseif ps.up_dn_looker==3
  % up looker only
  d.weight(d.izd,:)=nan;
end

%
% weightmin
%
disp(['    Minimum weight  ',num2str(ps.shear_weightmin),' of data '])

is1 = sum(sum(isfinite(d.weight)));
w = double(d.weight>ps.shear_weightmin);
w = replace(w,w==0,nan);
is2 =  sum(sum(isfinite(w)));
disp(['    Will use ',int2str(is2/is1*100),' % of data '])


%
% compute central shear
%
% disp(['  compute central diff shear '])
iiok = [1:length(d.z)];
us = [NaN*d.ru(1,iiok);diff2(d.ru(:,iiok))./diff2(d.izm);NaN*d.ru(1,iiok)].*w;
vs = [NaN*d.rv(1,iiok);diff2(d.rv(:,iiok))./diff2(d.izm);NaN*d.rv(1,iiok)].*w;
ws = [NaN*d.rw(1,iiok);diff2(d.rw(:,iiok))./diff2(d.izm);NaN*d.rw(1,iiok)].*w;


%
% loop over final depth profile
%
if isfield(dr,'z')==1
  z = dr.z;
else
  izmax = -min(d.z);
  izmax = -nmax(-[izmax,p.zbottom]);
  z = [idz*0.5:idz:izmax];
end


ds.usm = zeros(length(z),1)*NaN;
ds.vsm = ds.usm;
ds.wsm = ds.usm;
ds.use = ds.usm;
ds.vse = ds.usm;
ds.wse = ds.usm;
ds.nn = ds.usm;
ds.z = ds.usm;


il = length(z);


iv = 0;
%disp([' average shear estimates '])
for in=1:il
  i = z(in);
  iv = iv+1;
  % up and down
  i2 = find(abs(-d.izm-i-idz/2) <= idz);
  i1 = i2(find(isfinite(us(i2)+vs(i2))));
  ds.nn(iv) = length(i1);
  if ds.nn(iv) > 2
    usmm = median(us(i1)');
    ussd1 = std(us(i1));
    ii = i1(find(abs(us(i1)-usmm)<stdf*ussd1)); 
    if length(ii)>1
      ds.usm(iv) = mean(us(ii));
      ds.use(iv) = std(us(ii));
    end

    vsmm = median(vs(i1)');
    vssd1 = std(vs(i1));
    ii = i1(find(abs(vs(i1)-vsmm)<stdf*vssd1)); 
    if length(ii)>1
      ds.vsm(iv) = mean(vs(ii));
      ds.vse(iv) = std(vs(ii));
    end

    wsmm = median(ws(i1)');
    wssd1 = std(ws(i1));
    ii = i1(find(abs(ws(i1)-wsmm)<stdf*wssd1)); 
    if length(ii)>1
      ds.wsm(iv) = mean(ws(ii));
      ds.wse(iv) = std(ws(ii));
    end

  end

  % save depth vector
  ds.z(iv) = i;

  % end of shear loop
end


%
% integrate shear profile (from bottom up)
%

%ii = find(isnan(ds.usm));
%ds.usm(ii) = 0;
%ii = find(isnan(ds.vsm));
%ds.vsm(ii) = 0;
%ii = find(isnan(ds.wsm));
%ds.wsm(ii) = 0;
ds.usm = replace( ds.usm, isnan(ds.usm), 0 );
ds.vsm = replace( ds.vsm, isnan(ds.vsm), 0 );
ds.wsm = replace( ds.wsm, isnan(ds.wsm), 0 );


ds.ur = flipud(cumsum(flipud(ds.usm)))*idz;
ds.vr = flipud(cumsum(flipud(ds.vsm)))*idz;
ds.wr = flipud(cumsum(flipud(ds.wsm)))*idz;
ds.ur = ds.ur-mean(ds.ur);
ds.vr = ds.vr-mean(ds.vr);
ds.wr = ds.wr-mean(ds.wr);

if isfield(d,'down')
  dz = 2*abs(mean(diff(d.zd)));
  fac = 1/tan(d.down.Beam_angle*pi/180)*sqrt(2)*dz;
  ds.ensemble_vel_err = ds.wse*fac;
  dr.ensemble_vel_err = ds.wse*fac;
end

if nargin>3
  dr.u_shear_method = ds.ur;
  dr.v_shear_method = ds.vr;
  dr.w_shear_method = ds.wr;
  uds = nstd(dr.u-mean(dr.u)-ds.ur);
  vds = nstd(dr.v-mean(dr.v)-ds.vr);
  uvds = sqrt(uds.^2+vds.^2);
  if uvds>nmean(dr.uerr)
    warn = ('>   Increased error because of shear - inverse difference');
    disp(warn)
    if uvds>1.5*nmean(dr.uerr) & uvds>0.1
      messages.warn = strvcat(messages.warn,warn);
    end
    dr.uerr = dr.uerr/nmean(dr.uerr)*uvds/1.5;
  end
end

%--------------------------------------------------

function x = diff2(x,k,dn)
%DIFF2   Difference function.  If X is a vector [x(1) x(2) ... x(n)],
%       then DIFF(X) returns a vector of central differences between
%       elements [x(3)-x(1)  x(4)-x(2) ... x(n)-x(n-2)].  If X is a
%	matrix, the differences are calculated down each column:
%       DIFF(X) = X(3:n,:) - X(1:n-2,:).
%	DIFF(X,n) is the n'th difference function.

%	J.N. Little 8-30-85
%	Copyright (c) 1985, 1986 by the MathWorks, Inc.

if nargin < 2,	k = 1; end
if nargin < 3,	dn = 2; end
for i=1:k
	[m,n] = size(x);
	if m == 1
                x = x(1+dn:n) - x(1:n-dn);
	else
                x = x(1+dn:m,:) - x(1:m-dn,:);
	end
end

