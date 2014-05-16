function [data,params,values,messages]=calc_depth(data,params,values,messages)
% function [data,params,values,messages]=calc_depth(data,params,values,messages)
%
% Determine the depth of the LADCP system during the cast
% Different methods are employed according to the available
% data. Preferred method is to use parallel 
% high quality CTD pressure data.
% If that is not available, lower quality LADCP pressure data
% is used. If that isn't available either, vertical velocities
% are integrated by an inverse method.
%
% version 0.1	last change 27.7.2005

% Martin Visbeck December 2002, LDEO
% G.Krahmann, IFM-GEOMAR, Jul 2005
% major changes in Jul 2008   0.7 -> 0.8 
%

%%
% general function info
%
disp(' ')
disp('CALC_DEPTH:  determine the depth of the LADCP during cast')
values = setdefv(values,'ladcpdepth',-1);


%%
% calculate z the simple way by integrating w and correcting for
% preset start, end, and maximum depth
%
% these depths will not be used for any calculation but 
% are a good indicator for the quality of the LADCP record
%
disp(['    Calculating simple LADCP depth from w'])
if length(params.all_trusted_i)>1
  dw = meanmediannan( data.rw(params.all_trusted_i,:), 1 );
elseif length(params.all_trusted_i)==1
  dw = data.rw(params.all_trusted_i,:);
else
  dw = nmedian(data.rw,1);
end

%
%  set non finite w to zero
%
dw = replace( dw, ~isfinite(dw), 0 );

%
% get time difference for w-values
%
dt = diff(data.time_jul)*24*3600;
dt = mean([dt([1,1:end]);dt([1:end,end])]);

%
%  integrate result
%
zz = cumsum(dw.*dt);

%
%  make sure that start and end depth are as wanted
%
zz = zz-linspace(-nmax([0 params.zpar(1)]),...
	-nmax([0 params.zpar(3)])+zz(end),length(zz));

%
%  set maxdepth to be as requested
%
if isfinite(params.zpar(2))
  zz = zz/max(zz)*params.zpar(2);
end

%
% store these depths
%
data.ladcp_z_simple = -zz; 

%
% us as first guess if no CTD time series are present
%

if isempty(data.ctdtime_data) & values.ladcpdepth<0
 data.z=-zz;
end

%
% find index and depth of deepest point of profile
%
[values.maxdepth,ibot] = max(-data.z);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do the serious depth calculation
%

%
% for very shallow stations turn off surface detection
%
% don't know why that   GK
%
%if values.maxdepth<100 
%  params.surfdist = 0; 
%  disp('>   Found shallow station, am disabling surface detection ')
%end


figure(2)
clf

%
% first detect the surface of the down cast
%
% 

ind = [1:ibot];
iok = ind(find(data.z(ind)>-200));
if length(iok)>2
  iok2 = [1:iok(end)];
  subplot(221)
  plot(iok2,data.z(iok2),'.r',iok2,iok2*0,'-k'),
  hold on
  title('Best Depth (.red) Surface Signal (.blue)')
  xlabel('time in ensembles')
  ylabel('depth in meter')
  axis tight
  ax = axis;
  ax(4) = max(ax(4),20);
  axis(ax);
end

% then up cast
ind = [ibot:length(data.z)];
iok = ind(find(data.z(ind)>-200));
if length(iok)>2
  iok2 = iok(1):length(data.z);
  subplot(222)
  plot(iok2,data.z(iok2),'.r',iok2,iok2*0,'-k'),
  hold on
  title('Best Depth (.red) Surface Signal (.blue)')
  xlabel('time in ensembles')
  ylabel('depth in meter')
  axis tight
  ax = axis;
  ax(4) = max(ax(4),20);
  axis(ax);
end


%
% extract surface distance to find start depth
%
if isfield(data,'hsurf') & params.surfdist
  if sum(isfinite(data.hsurf))>10
    [zmax,ibot] = max(-data.z);
    % start depth
    ii = find(isfinite(data.hsurf(1:ibot)));
    % limit to surface detections within 30m of answer
    iok = ii(find(data.z(ii)>-200 & data.z(ii)<-40 &...
         abs(data.hsurf(ii)+data.z(ii))<30));
    if ~isempty(iok)
      ddoz = data.hsurf(iok)';
      ddoi = iok';
      % plot of surface detection
      subplot(221)
      plot(iok,data.hsurf(iok)+data.z(iok),'b.'),     
      axis tight
      ax = axis;
      ax(4) = max(ax(4),20);
      axis(ax);
    else
      disp('>   Found no surface reflections down-cast')
    end

    % surface distance to find end depth
    ii = find(isfinite(data.hsurf(ibot:end)))+ibot-1;
    % limit to surface detections within 30m of answer
    iok = ii(find(data.z(ii)>-200 & data.z(ii)<-40 &...
                abs(data.hsurf(ii)+data.z(ii))<30));
    if ~isempty(iok)
      dupz = data.hsurf(iok)';
      dupi = iok';
      % plot of surface detection
      subplot(222)
      plot(iok,data.hsurf(iok)+data.z(iok),'.'),
      axis tight
      ax = axis;
      ax(4) = max(ax(4),20);
      axis(ax);
    else
      disp('>   Found no surface reflections up-cast')
    end
  end 
end 

disp(['    Maximum depth from int W is : ',int2str(nmax(-data.z))])
disp(['    Preset maximum depth is     : ',int2str(params.zpar(2))])


%
% determine bottom depth
%
subplot(212)
iok = find((max(-data.z)+data.z)<200);
iok2 = [iok(1):iok(end)];
plot(iok2,data.z(iok2),'.r'),
%if params.zbottom==0
if 1
  if sum(isfinite(data.hbot))>10
    % look for bottom only close to deepest CTD depth
    iok = find((max(-data.z)+data.z)<200 & data.hbot>10 );
    if length(iok)>5
      % fit polynom to bottom depth time series
      zbottom = median(data.hbot(iok)-data.z(iok));
      zbottomerr = zbottom-(data.hbot(iok)-data.z(iok)) ;
      [dum,is] = sort(abs(zbottomerr));
      is = is(1:fix(length(is)/2));
% the iok-ibot stuff is new, GK
      c1 = polyfit(iok(is)-ibot,data.hbot(iok(is))-data.z(iok(is)),1);
      zbottomerr = polyval(c1,iok-ibot)-(data.hbot(iok)-data.z(iok)) ;
      [dum,is] = sort(abs(zbottomerr));
      is = is([1:fix(length(is)/2)]);
      iok = iok(find(abs(zbottomerr)<2*std(zbottomerr(is)) |...
		 abs(zbottomerr)<30 ));
      c = polyfit(iok-ibot,data.hbot(iok)-data.z(iok),2);
      zbottomerr = polyval(c,iok-ibot)-(data.hbot(iok)-data.z(iok)) ;
      zbottom = polyval(c,ibot-ibot);
      params.zbottom = zbottom;
      % save bottom distances for inversion
      if exist('values.depth_iter','var')
          dbotdz = (data.hbot(iok)+polyval(c,iok-ibot)+data.z(ibot))';
          values.depth_iter = values.depth_iter+1;
      else
          dbotdz = (data.hbot(iok)+polyval(c1,iok-ibot)+data.z(ibot))';
          values.depth_iter=1;
      end
      dboti = iok';
      params.zbottomerror = nmedian( abs(zbottomerr) ); 

      % temporary plot of bottom detection
      iok2 = [iok(1):iok(end)];
      plot(iok,-data.hbot(iok)+data.z(iok),'b.',iok2,data.ladcp_z_simple(iok2),'-g'),
      hold on
      
      plot(iok2,iok2*0-params.zbottom,'--k') 
      plot(iok2,-polyval(c,iok2-ibot),'-b')
      title('Bottom calculated (--k) LADCP depth (.r) Bottom Signal (.b) Integrated W depth (-g)')
      xlabel('time in ensembles')
      ylabel('depth in meter')

      % remove suspicious bottom track data
      axis tight
      ax = axis;
      text(ax(1)+abs(diff(ax(1:2)))*.15,ax(3)+abs(diff(ax(3:4)))*.8,...
        ['bottom at: ',int2str(params.zbottom),' [m]   ADCP was ',...
        int2str(params.zbottom-max(-data.z)),' m above bottom'])
      ibad = [1:length(data.hbot)];
      % good data are
      ibad(iok) = [];
      data.hbot(ibad) = NaN;
      data.bvel(ibad,:) = NaN;
    else
      disp('    Not enough data to determine water depth.')
      zbottom = nan;
      params.zbottom = nan;
      params.zbottomerror = nan;
    end
  else
    disp('    Not enough data to determine water depth.')
    zbottom = NaN;
    params.zbottom = nan;
    params.zbottomerror = nan;
  end
else
  disp('    Not determining water depth. (no data or turned off)')
  zbottom = NaN;
  params.zbottom = NaN;
  params.zbottomerror = NaN;
end


%
% check if bottom is shallower that maxctd-depth an
%
if ((zbottom-values.maxdepth<-(values.maxdepth*0.01+10) &...
	 isfinite(zbottom)) |...
       	params.zbottomerror > 20 )
  disp('    Found no bottom')
  disp(['      given maximum profile depth : ',int2str(values.maxdepth)])
  disp(['      extracted bottom depth      : ',int2str(zbottom)])
  disp(['      bottom depth error          : ',int2str(params.zbottomerror)])
  params.zbottom = NaN;
end
params.zbottom = zbottom;
disp(['    Found bottom at ',int2str(params.zbottom),' +/- ',...
                                int2str(params.zbottomerror),' m'])
if (params.zbottom<values.maxdepth)
  disp('    Extracted a bottom within 20m above given maximum profile depth')
end


%
% first look whether CTD time data exists
%
if ~isempty(data.ctdtime_data)

  disp('    Using CTD pressure record for depth')

  data.z_ladcp = -zz;
  dz = data.z_ladcp-data.z;
  params.ladcpr_CTD_depth_std = [nmean(dz), nstd(dz)];
  disp(['    LADCP minus CTD depth mean: ',...
        num2str(params.ladcpr_CTD_depth_std(1)),...
          '  std: ',num2str(params.ladcpr_CTD_depth_std(2))]);
  values.ladcpdepth = 0;
  
%
% no CTD time data, continue with other ways
%
else

  %
  % some instruments contain a pressure sensor
  % first check for that one
  % 
  if rms(data.pres)>0

    disp(['    Using LADCP pressure data'])
    zz = sw_dpth(data.pres,values.lat*ones(size(data.pres)));
    data.z = -zz;
    data.z_ladcp = -zz;
    data.z_ladcp_p = -zz;
    data.izmflag = data.rw*0;
    values.ladcpdepth = 3;
    plot(iok2,data.z(iok2),'-c')
  else
    params.weight_ladcp_pressure=0;
  end

  %
  % last resort, use the LADCP itself to get the depth
  %    
  if rms(data.pres)==0 | params.weight_ladcp_pressure>0
    disp(['    Calculating inverted LADCP depth from w'])

    % initialize a matrix where below bottom and above surface data
    % get flagged as bad
    data.izmflag = data.rw*0;

    % set base matrix
    dw = data.rw(params.all_trusted_i,:)+data.izmflag(params.all_trusted_i,:);
    [d1,A1,ibot] = dinset(dw,dt);

    % set boundary conditions for inversion
    if isfinite(params.zpar(1))
      d1 = [d1;params.zpar(1)];
    else 
      d1 = [d1;params.cut];
    end
    A1(length(d1),1) = 1;

    if isfinite(params.zpar(2))
      d1 = [d1;params.zpar(2)*10];
      A1(length(d1),ibot) = 10;
    end

    if isfinite(params.zpar(3))
      d1 = [d1;params.zpar(3)];
    else
      d1 = [d1;params.cut];
    end
    A1(length(d1),end) = 1;

    % add surface/bottom reflections if present
 
    if exist('ddoz','var')
      sfac = 0.01;
      [ld,lz] = size(A1);
      d1 = [d1;ddoz*0.1];
      A2 = sparse(1:length(ddoi),ddoi,sfac);
      A2(1,lz) = 0;
      A1 = [A1;A2];
      disp(['    Used ',int2str(length(ddoi)),...
	' surface detections down cast'])
    end

    if exist('dupz','var')
      sfac = 0.01;
      [ld,lz] = size(A1);
      d1 = [d1;dupz*0.1];
      A2 = sparse(1:length(dupi),dupi,sfac);
      A2(1,lz) = 0;
      A1 = [A1;A2];
      disp(['    Used ',int2str(length(dupi)),...
	' surface detections up cast'])
    end

    if exist('dbotdz','var') 
      bfac = 0.1;
      [ld,lz] = size(A1);
      d1 = [d1;dbotdz*bfac];
      ix = [1:length(dboti)]';
      A2 = sparse([ix,ix],[dboti;ix*0+ibot],[ix*0-bfac;ix*0+bfac]);
      A2(1,lz) = 0;
      A1 = [A1;A2];
      disp(['    Used ',int2str(length(dboti)),...
	' bottom distances'])
    end
    
    
    % add surface/bottom reflections if present
    
    if params.weight_ladcp_pressure>0
      [ld,lz] = size(A1);
      ix = find(isfinite(data.z_ladcp_p));
      d1 = [d1;-data.z_ladcp_p(ix)'*params.weight_ladcp_pressure];
      A2 = sparse(1:length(ix),ix,params.weight_ladcp_pressure);
      A2(1,lz) = 0;
      A1 = [A1;A2];
      disp(['    Used ',int2str(length(ix)),...
	' ADCP pressure data'])
    end

    % require z to be smooth (needed otherwise ill constrained)
    [A1,d1] = dismoo(A1,d1,0.05);

    % solve for best depth time series
    zz = lesqchol(d1,A1)';
    iok = find((max(zz)-zz)<200);
    iok2 = [iok(1):iok(end)];
    plot(iok2,-zz(iok2),'.g'),
   
    
    data.z = -zz;
    data.z_ladcp = -zz;
    
    

% ??? wrong sign in the above two lines ?
%    data.z = zz;
%    data.z_ladcp = zz;

    values.ladcpdepth = 2;
  end
 
end
plot(iok2,data.z(iok2),'.r')

%-----------------------------------------------------------------
function [d,A,izbot]=dinset(dw,dt)
% function [d,A]=dinset(dw,dt)
% set up sparse Matrix for depth inversion
%
[nb,nt] = size(dw);

% find bottom roughly and devide cast in down and up trace
if nb>1
  wm = nmedian(dw);
else
  wm = dw;
end
ii = find(~isfinite(wm));
wm(ii) = 0;
zz = cumsum(wm.*dt);

[zbot,izbot] = nmax(zz);
%disp(['  using as deepest point of profile: ',int2str(zbot),'m at index ',...
%	int2str(izbot)])
%
% I think this one is misleading   GK


ido = [1:izbot];
iup = [izbot+1:length(dt)];

% 
dtm = repmat(dt,[nb 1]);
izm = repmat([1:nt],[nb 1]);

d = reshape(dw.*dtm,nb*nt,1);
izv = reshape(izm,nb*nt,1);

ibad = find((izv-1)<1 | (izv+1)>nt);
d(ibad) = [];
izv(ibad) = [];

iweak = find(~isfinite(d) | (izv-1)<1 | (izv+1)>nt);

it = [1:length(d)]';
i1 = it*0+0.5;
d(iweak) = 0;
i1(iweak) = 0.01;

A = sparse([it;it],[izv+1;izv-1],[i1;-i1]);
A(1,nt) = 0;


%-------------------------------------------------------------------
function [A,d]=dismoo(A,d,fs0,cur);
% function [A,d]=dismoo(A,d,fs0,cur);
%
% smooth results by minimizing curvature
% also smooth if elements are not constrained
%
if nargin<3 
  fs0 = 1; 
end
if nargin<4 
  cur = [-1 2 -1]; 
end

[ld,ls] = size(A);
fs = full(sum(abs(A)));
fsm = max( median(fs), 0.1 );
% find ill constrained data
ibad = find(fs<fsm*0.1);

% increase weight for poorly constrained data
fs = max(fs,fsm*0.1);
fs = fsm^2./fs * fs0(1);


if length(ibad)>0
  % set ill constrainded data to a minimum weight
  fs(ibad) = max(fs(ibad),0.5);
  if fs0==0
    disp(['    Found ',int2str(length(ibad)),...
	' ill constrained elements will smooth '])
  else
    disp(['    Found ',int2str(length(ibad)),...
	' ill constrained elements'])
  end
end

if sum(fs>0)>0

  cur = cur-mean(cur);

  lc = length(cur);
  lc2 = fix(lc/2);
  fs2 = fs([lc2+1:end-lc2]);
  inc = [1:length(cur)]-lc2;

  ii = find(fs2>0);
  % find how many smooth constraints to apply

  if length(ii)>0 
    [i1,i2] = meshgrid(inc,ii+lc2-1);
    [curm,fsm] = meshgrid(cur,fs2(ii));

    As = sparse(i2,i1+i2,curm.*fsm);
    [lt,lm] = size(A);
    if size(As,2)<lm 
      As(1,lm) = 0;
    end
    A = [A;As];

  end

  % smooth start and end of vector
  for j=1:lc2
    j0 = j-1;
    [lt,lm] = size(A);
    if fs(1+j0)>0
      A(lt+1,[1:2]+j0) = [2 -2]*fs(1+j0);
    end
    if fs(end-j0)>0
      A(lt+2,end-[1,0]-j0) = [-2 2]*fs(end-j0);
    end
  end

  [lt,lm] = size(A);
  d(lt) = 0;
else
  disp('    No smoothness constraint applied ')
end

