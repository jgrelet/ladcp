function [hdr,data,hdr1,data1] = read_sbe_cnv(fname,start_scan);
% function [hdr,data,hdr1,data1] = read_sbe_cnv(fname,[start_scan]);
%
% load a Seabird CTD file (CNV)
%
% extracts info from header and from data array
% the variables  p,t,c,o,chlorophyl,timeJ,scan,lat,lon,pump, and flag
% are being extracted
% from the header position and time both from NMEA and from 
% typed info are extracted
%
% input  :  fname            - file name
%           start_scan  [1]  - scan at which the downcast begins
%
% output :  hdr              - structure containing the info of the file
%           data             - structure containing the data of the file
%           hdr1             - character array of the file header (like hdrload)
%           data1            - data array (like hdrload)
%   
% uses :	hdrload.m, nans.m
%
% version 0.11   last change 17.01.2009

% G.Krahmann, IFM-GEOMAR, Dec 2006

% added HaardtC as variable             GK, 30.1.2007   0.1-->0.2
% reorganisation of data                GK, 30.5.2007   0.2-->0.3
% added start_scan                      GK, 05.06.2007  0.3-->0.4
% flagged values to NaN                 GK, 27.02.2008  0.4-->0.5
% changed output fields                 GK, 02.03.2008  0.5-->0.6
% read sensor S/N                       GK, 30.09.2008  0.6-->0.7
% single sensor v5 diff                 GK, 03.11.2008  0.7-->0.8
% new time variable data.datenum that is consistent with
% NMEA date and monotonous in time      GK, 13.01.2009  0.8-->0.9
% hdr.badflag now a number              GK, 14.01.2009  0.9-->0.10
% 'turn off' pump 20 scans earlier      GK, 17.01.2009  0.10-->0.11


%
% give help
%
if nargin==0
  help read_sbe_cnv
  return
end


%
% load data and header
%
[hdr1,data1] = hdrload(fname);
for n=1:size(hdr1,1)
  ind = findstr(hdr1(n,:),'\');
  if ~isempty(ind)
    hdr1(n,ind) = ' ';
  end
end


%
% cut to data starting at start_scan
%
if nargin>1
    data1 = data1(start_scan:end,:);
end


%
% parse Seabird header
%
hdr.nmea_lat_str = extract_string('* NMEA Latitude =',hdr1);
hdr.nmea_lon_str = extract_string('* NMEA Longitude =',hdr1);
hdr.nmea_utc_str = extract_string('* NMEA UTC (Time) =',hdr1);
hdr.lat = extract_string('** Latitude :',hdr1);
hdr.lon = extract_string('** Longitude :',hdr1);
hdr.station = extract_string('** Station :',hdr1);
hdr.waterdepth = extract_string('** Water Depth :',hdr1);
hdr.badflag = str2num(extract_string('# bad_flag =',hdr1));
hdr.system_upload = extract_string('* System UpLoad Time =',hdr1);
if isempty(hdr.nmea_lat_str)
  hdr.nmea_lat = [];
  hdr.nmea_lon = [];
  hdr.nmea_utc = [];
  hdr.nmea_utc_datenum = [];
else
  hdr = parse_nmea(hdr);
end


%
% parse System Upload Time
%
strs = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
dat = hdr.system_upload;
dummy(2) = strmatch(dat(1:3),strs);
dummy(1) = sscanf(dat(8:11),'%d')';
dummy(3) = sscanf(dat(5:6),'%d')';
dummy(4) = sscanf(dat(13:14),'%d')';
dummy(5) = sscanf(dat(16:17),'%d')';
dummy(6) = sscanf(dat(19:20),'%d')';
hdr.system_upload_datenum = datenum(dummy);
hdr.system_upload_day = datenum(dummy.*[0,1,1,0,0,0]);


%
% replace flagged values by NaN
%
data1 = nans(data1,nan,hdr.badflag,'==');


%
% loop over first 30 entries in data description, as these contain the names and columns of the variables
%
mll_pri = 0;
mll_sec = 0;
for n=0:30
    %
    % extract variable names
    %
    hdr = setfield(hdr,['name',int2str(n)],...
          extract_string(['# name ',int2str(n),' ='],hdr1));

    %
    % extract the data for
    fld = getfield(hdr,['name',int2str(n)]);

    %
    % check whether the new variable is a standard one
    % if so, extract the data
    %
    if ~isempty(findstr(fld,'prDM'))
        data.p = data1(:,n+1);
    elseif ~isempty(findstr(fld,'sal00'))
        data.s_pri = data1(:,n+1);
    elseif ~isempty(findstr(fld,'sal11'))
        data.s_sec = data1(:,n+1);
    elseif ~isempty(findstr(fld,'t090C'))
        data.t_pri = data1(:,n+1);
    elseif ~isempty(findstr(fld,'t190C'))
        data.t_sec = data1(:,n+1);
    elseif ~isempty(findstr(fld,'c0S/m'))
        data.c_pri = data1(:,n+1);
    elseif ~isempty(findstr(fld,'c1S/m'))
        data.c_sec = data1(:,n+1);
    elseif ~isempty(findstr(fld,'sbeox0ML/L'))
        data.o_pri = data1(:,n+1);
        mll_pri = 1;
    elseif ~isempty(findstr(fld,'sbeox1ML/L'))
        data.o_sec = data1(:,n+1);
        mll_sec = 1;
    elseif ~isempty(findstr(fld,'timeJ'))
        data.timej = data1(:,n+1);
    elseif ~isempty(findstr(fld,'pumps'))
        data.pumps = data1(:,n+1);
        % 
        % here we turn off the pump 20 scans earlier than in
        % in the raw data, as the 'out of water' recognition reacts
        % late and thus often has wrong conductivities when the
        % CTD is already out of the water
        % 
        ind = find(data.pumps(1:end-1)==1 & data.pumps(2:end)==0);
        if ~isempty(ind) 
          if ind(end)>20
            data.pumps(ind(end)+[-20:0]) = 0;
          end
        end  
    elseif ~isempty(findstr(fld,'scan'))
        data.scan = data1(:,n+1);
    elseif ~isempty(findstr(fld,'latitude'))
        data.latitude = data1(:,n+1);
    elseif ~isempty(findstr(fld,'longitude'))
        data.longitude = data1(:,n+1);
    elseif ~isempty(findstr(fld,'flag'))
        data.flag = data1(:,n+1);
    elseif ~isempty(findstr(fld,'haardtC'))
        data.haardtc = data1(:,n+1);
    end
end


%
% create dummy data, if only single sensors have been used
%
if ~isfield(data,'t_sec')
  data.t_sec = data.t_pri;
  data.t_sec_is_dummy = 1;
end
if ~isfield(data,'c_sec')
  data.c_sec = data.c_pri;
  data.c_sec_is_dummy = 1;
  data.s_sec = nan*data.c_sec;
end
if ~isfield(data,'o_sec')
  data.o_sec = data.o_pri;
  data.o_sec_is_dummy = 1;
end


%
% convert to oxygen to mumol/kg from ml/l
%
if mll_pri==1
  sth = sw_pden(data.s_pri,data.t_pri,data.p,0);
  data.o_pri = data.o_pri*44600./sth;
end
if mll_pri==1
  sth = sw_pden(data.s_sec,data.t_sec,data.p,0);
  data.o_sec = data.o_sec*44600./sth;
end


%
% set all values to NaN when the pump is off
%
if isfield(data,'pumps')
  bad = find(data.pumps==0);
  data.t_pri(bad) = nan;
  data.t_sec(bad) = nan;
  data.c_pri(bad) = nan;
  data.c_sec(bad) = nan;
  data.o_pri(bad) = nan;
  data.o_sec(bad) = nan;
end


%
% extract information of sensor IDs
%
hdr.p_sn = 0;
hdr.t_pri_sn = 0;
hdr.t_sec_sn = 0;
hdr.c_pri_sn = 0;
hdr.c_sec_sn = 0;
hdr.o_pri_sn = 0;
hdr.o_sec_sn = 0;
hdr.chl_sn = 0;
hdr.tr_sn = 0;
for n=1:size(hdr1,1)
  if strcmp(hdr1(n,3:8),'sensor')
    if ~isempty(findstr(hdr1(n,:),'pressure'))
      ind = findstr(hdr1(n,:),',');
      hdr.p_sn = str2num(hdr1(n,ind(1)+1:ind(2)-1));
    elseif ~isempty(findstr(hdr1(n,:),'temperature, primary'))
      ind = findstr(hdr1(n,:),',');
      hdr.t_pri_sn = str2num(hdr1(n,ind(2)+1:ind(3)-1));
    elseif ~isempty(findstr(hdr1(n,:),'temperature, secondary'))
      ind = findstr(hdr1(n,:),',');
      hdr.t_sec_sn = str2num(hdr1(n,ind(2)+1:ind(3)-1));
    elseif ~isempty(findstr(hdr1(n,:),'temperature,'))
      ind = findstr(hdr1(n,:),',');
      hdr.t_pri_sn = str2num(hdr1(n,ind(1)+1:ind(2)-1));
    elseif ~isempty(findstr(hdr1(n,:),'conductivity, primary'))
      ind = findstr(hdr1(n,:),',');
      hdr.c_pri_sn = str2num(hdr1(n,ind(2)+1:ind(3)-1));
    elseif ~isempty(findstr(hdr1(n,:),'conductivity, secondary'))
      ind = findstr(hdr1(n,:),',');
      hdr.c_sec_sn = str2num(hdr1(n,ind(2)+1:ind(3)-1));
    elseif ~isempty(findstr(hdr1(n,:),'conductivity,'))
      ind = findstr(hdr1(n,:),',');
      hdr.c_pri_sn = str2num(hdr1(n,ind(1)+1:ind(2)-1));
    elseif ~isempty(findstr(hdr1(n,:),'Oxygen, SBE, primary'))
      ind = findstr(hdr1(n,:),',');
      hdr.o_pri_sn = str2num(hdr1(n,ind(3)+1:ind(4)-1));
    elseif ~isempty(findstr(hdr1(n,:),'Oxygen, SBE, secondary'))
      ind = findstr(hdr1(n,:),',');
      hdr.o_sec_sn = str2num(hdr1(n,ind(3)+1:ind(4)-1));
    elseif ~isempty(findstr(hdr1(n,:),'Oxygen, SBE,'))
      ind = findstr(hdr1(n,:),',');
      hdr.o_pri_sn = str2num(hdr1(n,ind(2)+1:ind(3)-1));
    end
  end
end


%
% compare Seabird Julian day (system upload time) with NMEA day
% and create a new time vector on the basis of Matlab's datenum
% also interpolate Seabird's time so that it is continuous and
% not always 2 same-times after another (too few digits in timej)
%
if ~isempty(hdr.nmea_utc)
  d0 = data.timej(1);
  d1 = hdr.system_upload_datenum-datenum(datevec(hdr.system_upload_datenum).*...
    [1,0,0,0,0,0]);
  d2 = datenum(hdr.nmea_utc)-datenum([hdr.nmea_utc(1),0,0,0,0,0]);
  data.datenum = data.timej-hdr.system_upload_day + datenum(hdr.nmea_utc(1:3));
  if abs(d1-d2)>1/1440
    fprintf(1,'\n',[]);
    disp(datestr(d1,31))
    disp(datestr(d2,31))
    warning('NMEA and CTD computer time differ by more than 1 minute')
    data.timej = data.timej-hdr.system_upload_day + datenum([0,hdr.nmea_utc(2:3)]);
  end
  jump = find(diff(data.datenum(:)')>1e-5);
  lim = [[1,jump];[jump+1,length(data.datenum)]];
  for n=1:size(lim,2)
    data.datenum(lim(1,n):lim(2,n)) = linspace(data.datenum(lim(1,n)),...
      data.datenum(lim(2,n)),lim(2,n)-lim(1,n)+1);
  end
else
  warning('no NMEA data in header, proper time handling not implemented!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = extract_string(str,strs)
ind = strmatch(str,strs);
ind2 = findstr(strs(ind,:),'=');
if ~isempty(ind2)
    data = strs(ind,ind2+1:end);
    data = deblank(data);
    data = fliplr(deblank(fliplr(data)));
else
    ind2 = findstr(strs(ind,:),':');
    if ~isempty(ind2)
        data = strs(ind,ind2+1:end);
        data = deblank(data);
        data = fliplr(deblank(fliplr(data)));
    else
        data = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr] = parse_nmea(hdr)

if isfield(hdr,'nmea_lat_str')
  dummy = sscanf(hdr.nmea_lat_str,'%d %f')';
  if findstr(hdr.nmea_lat_str,'N')
    hdr.nmea_lat = dummy(1)+dummy(2)/60;
  else
    hdr.nmea_lat = -dummy(1)-dummy(2)/60;
  end
end
if isfield(hdr,'nmea_lon_str')
  dummy = sscanf(hdr.nmea_lon_str,'%d %f')';
  if findstr(hdr.nmea_lon_str,'E')
    hdr.nmea_lon = dummy(1)+dummy(2)/60;
  else   
    hdr.nmea_lon = -dummy(1)-dummy(2)/60;
  end
end
if isfield(hdr,'nmea_utc_str')
  strs = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
  dat = hdr.nmea_utc_str;
  dummy(2) = strmatch(dat(1:3),strs);
  dummy(1) = sscanf(dat(8:11),'%d')';
  dummy(3) = sscanf(dat(5:6),'%d')';
  dummy(4) = sscanf(dat(14:15),'%d')';
  dummy(5) = sscanf(dat(17:18),'%d')';
  dummy(6) = sscanf(dat(20:21),'%d')';
  hdr.nmea_utc_datenum = datenum(dummy);
  hdr.nmea_utc = dummy;
end
