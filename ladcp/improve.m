function [data,params,values,messages] = improve(data,params,values,messages)
% function [data,params,values,messages] = improve(data,params,values,messages)
%
% improve the data
% fix all kinds of errors and problems, apply rotations, etc.
%
% input  :	data		- LADCP data structure
%		params		- LADCP parameter structure
%		values		- LADCP value structure
%		messages	- array of warnings
%
% output :	data		- modified LADCP data structure
%		params		- LADCP parameter structure
%		values		- modified LADCP value structure
%		messages	- array of warnings
%
% version 0.5	last change 20.09.2007

% G.Krahmann, LDEO Nov 2004

% added time dependence in magdev		GK, Mar 2007  	0.1-->0.2
% error check for empty raw data		GK, Aug 2007  	0.2-->0.3
% threw out too much wlim data			GK, Sep 2007	0.3-->0.4
% stop processing, if too many tilt>tiltmax	GK, Sep 2007	0.4-->0.5


%
% general function start info
%
disp(' ')
disp('IMPROVE:  apply various improvements to the data')


%
% check for outliers in the whole data set
%
[data,params] = outlier(data,params,values);


%
% apply magnetic deviation
%
if values.lat~=0 | values.lon~=0
  [a,b] = gregoria(values.start_time);
  decyear = a+(values.start_time-julian([a,1,1,0,0,0]))/365.25;
  values.magdev = magdev(values.lat,values.lon,0,decyear);
  [data.ru,data.rv] = uvrot(data.ru,data.rv,values.magdev);
  [data.bvel(:,1),data.bvel(:,2)] = ...
	uvrot(data.bvel(:,1),data.bvel(:,2),values.magdev);
  disp(['    Applying magnetic deviation of: ',num2str(values.magdev)])
else
  disp('>   NOT correcting for magnetic deviation')
  disp('>     found lat=0 & lon=0 , assuming no position info available')
  disp('>     else, set position in cast_params.m')
  values.magdev = 0;
end


%
% apply vertical velocity limit
%
ind = findany(params.bins_d,params.trusted_i);
data.wd = meanmediannan(data.rw(ind,:),1);
w(data.izd,:) = meshgrid(data.wd,data.izd,2);
if values.up==1
  ind = findany(params.bins_u,params.trusted_i);
  data.wu = meanmediannan(data.rw(ind,:),1);
  w(data.izu,:) = meshgrid(data.wu,data.izu,2);
end
j = find(abs(data.rw-w) > params.wlim);
disp(['    Removing ',int2str(length(j)),...
	' values because of vertical speed deviates > ',...
	num2str(params.wlim),' m/s'])
data.ru(j) = NaN;
data.rv(j) = NaN;
data.rw(j) = NaN;
j = find(abs(data.bvel(:,3)'-w(data.izd(1),:)) > params.wlim);
data.bvel(j,:) = NaN;


% 
% apply horizontal velocity limit
%
vel = sqrt(data.ru.^2+data.rv.^2);
j = find(vel > params.vlim);
data.ru(j) = NaN;
data.rv(j) = NaN;
data.rw(j) = NaN;
if length(j) > 0
  disp(['    Removing ',int2str(length(j)),...
	' values because of horizontal speed > ',num2str(params.vlim),' m/s'])
end
vel = sqrt(data.bvel(:,1).^2+data.bvel(:,2).^2);
j = find(vel > params.vlim);
data.bvel(j,:) = NaN;


% 
% remove empty profiles at beginning and end
%
ind = zeros(1,size(data.rw,2));
jj = find(~isnan(nmean(data.rw)));
if isempty(jj)
  error('no not-NaN data in rw')
end
ind(jj(1):jj(end)) = 1;
values.firstlastindex = [jj(1),jj(end)];
data = cutstruct(data,ind);


%
% fix problems with switched beams on instrument
%
if isprop(params,'beam_switch')
  [data,params,values] = misc_switchbeams(data,params,values);
end


%
% fix problems with a compass
%
if params.fix_compass>0
  [data,params] = misc_fix_compass(data,params);
end


%
% apply pitch/tilt corrections
%
if length(params.tiltcor)>1
  pd.dpit = params.tiltcor(1);
  pd.drol = params.tiltcor(2);
  [data,params] = uvwrot(data,pd,1);
end

if length(params.tiltcor)>2
  pu.dpit = params.tiltcor(3);
  pu.drol = params.tiltcor(4);
  [data,params] = uvwrot(data,pu,0);
end


%
% remove profiles with too large tilts
%
ind = find( data.tilt > params.tiltmax(1) );
disp(['    Removing ',int2str(length(ind)),...
      ' profiles due to tilt larger than ',int2str(params.tiltmax(1))])
if length(ind)>length(data.tilt)*0.1;
  warn=(['>   ',int2str(length(ind)*100/length(data.tilt)),...
         '%  tilt > ',int2str(params.tiltmax(1))]);
  disp(warn)
  messages.warn = strvcat(messages.warn,warn)
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
  disp('  Continued processing will likely fail as there are too')
  disp('  many ensembles with a too large tilt.')
  disp('  You can force the continued processing by increasing the')
  disp('  values in params.tiltmax')
  disp(['  At the moment ensembles with angles larger than ',...
	int2str(params.tiltmax(1))])
  disp(['  subsequent angles differing more than ',...
	int2str(params.tiltmax(2))])
  disp('  are being removed. Increase these settings either in')
  disp('  cruise_params.m or cast_params.m')
  disp('  This will enable continued processing but possibly result')
  disp('  in not trustworthy profiles !')
  disp(' ')
  error('Stopping execution of LADCP processing')
end
data.weight(:,ind)=NaN;


%
% remove profiles with too large consecutive tilt differences
%
ind = find( data.tiltd > params.tiltmax(2) );
disp(['    Removing ',int2str(length(ind)),...
      ' profiles due to tilt difference larger than ',...
	int2str(params.tiltmax(2))])
data.weight(:,ind)=NaN;


