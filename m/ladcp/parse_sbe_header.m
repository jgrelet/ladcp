function [hdr] = parse_sbe_header(fname,hdr1);
% function [hdr] = parse_sbe_header(fname,hdr1);
%
% parse the header from a Seabird file, either cnv or btl
%
% input  :  fname            - file name
%           hdr1             - header as a text array
%
% output :  hdr              - structure containing the info of the file
%   
% uses :	nans.m
%
% version 0.9   last change 10.12.2012

% G.Krahmann, IFM-GEOMAR, Jun 2010

% taken out of read_sbe_cnv.m           GK, 09.06.2010  0.1
% bug in SN oxygen                      GK, 10.09.2010  0.1-->0.2
% handle missing SN in oxygen           GK, 21.09.2010  0.2-->0.3
% use datevec to interpret date string  GK, 22.09.2010  0.3-->0.4
% change missing SN in oxygen nan to 0  GK, 22.09.2010  0.4-->0.5
% fix problem when without NMEA         GK, 19.10.2010  0.5-->0.6
% more handling of no NMEA              GK, 11.04.2012  0.6-->0.7
% Islandia CTD sensors                  GK, 31.08.2012  0.7-->0.8
% bug in NMEA date                      GK, 10.12.2012  0.8-->0.9

%
% give help
%
if nargin==0
  help parse_sbe_header
  return
end


%
% parse Seabird header
%
hdr.nmea_lat_str = extract_string('**  NMEA Latitude =',hdr1);
if isempty(hdr.nmea_lat_str)
   hdr.nmea_lat_str = extract_string('* NMEA Latitude =',hdr1);
end
hdr.nmea_lon_str = extract_string('**  NMEA Longitude =',hdr1);
if isempty(hdr.nmea_lon_str)
   hdr.nmea_lon_str = extract_string('* NMEA Longitude =',hdr1);
end
hdr.nmea_utc_str = extract_string('**  NMEA UTC (Time) =',hdr1);
if isempty(hdr.nmea_utc_str)
   hdr.nmea_utc_str = extract_string('* NMEA UTC (Time) =',hdr1);
end
hdr.lat = extract_string('** Latitude :',hdr1);
hdr.lon = extract_string('** Longitude :',hdr1);
hdr.station = extract_string('** Station :',hdr1);
hdr.waterdepth = extract_string('** Water Depth :',hdr1);
hdr.bad_flag = extract_string('# bad_flag =',hdr1);
if ~isempty(findstr(hdr1(1,:),'19'))
  hdr.sbe_model = 19;
else
  hdr.sbe_model = 9;
end
if ~isempty(hdr.bad_flag)
  hdr.bad_flag = str2num(hdr.bad_flag);
end
hdr.system_upload = extract_string('* System UpLoad Time =',hdr1);
hdr.start_time    = extract_string('# start_time =',hdr1);
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
hdr.system_upload = hdr.system_upload(1:20);
dummy = datevec(hdr.system_upload,'mmm dd yyyy HH:MM:SS');
hdr.system_upload_datenum = datenum(dummy);
hdr.system_upload_sbe_doy = datenum(dummy.*[1,1,1,0,0,0])-datenum(dummy.*[1,0,0,0,0,0]);

%
%  parse System Start Time
%
hdr.start_time = hdr.start_time(1:20);
dummy = datevec(hdr.start_time,'mmm dd yyyy HH:MM:SS');
hdr.start_time_datenum = datenum(dummy);
hdr.start_time_doy = datenum(dummy.*[1,1,1,0,0,0])-datenum(dummy.*[1,0,0,0,0,0]);

%
% figure out the information to use for CTD's position and start
%
hdr.use_lat = 0;
hdr.use_lon = 0;
if ~isempty(hdr.nmea_lat)
  hdr.use_lat = hdr.nmea_lat;
  hdr.use_lon = hdr.nmea_lon;
  hdr.use_datenum = hdr.nmea_utc_datenum;
  if isnan(hdr.use_datenum)
    hdr.use_datenum = hdr.system_upload_datenum;
  end
elseif exist('mat/lat_lon.mat') | exist('mat\lat_lon.mat') | exist('lat_lon.txt')
  if exist('lat_lon.txt')
    latlon = load('lat_lon.txt');
  else
    load('mat/lat_lon.mat');
  end
  ind = findstr(fname,'_');
  stn = str2num(fname(ind(end)+[1:3]));
  ind = find(latlon(:,1)==stn);
  if isempty(ind)
    error('could not find position in lat_lon.txt')
  end
  hdr.use_lat = latlon(ind,2);
  hdr.use_lon = latlon(ind,3);
  hdr.use_datenum = hdr.system_upload_datenum;
else
  disp('create file lat_lon.txt')
  error('Found no position info. Using 0N 0E')
  hdr.use_lat = 0;
  hdr.use_lon = 0;
  hdr.use_datenum = hdr.system_upload_datenum;
end
hdr.use_datevec = datevec(hdr.use_datenum);


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
hdr.chl2_sn = 0;
hdr.par_sn = 0;
hdr.spar_sn = 0;
hdr.tr_sn = 0;

%
% check for old or new version of sensor information
%
new_style = 0;
for n=1:size(hdr1,1)
  if ~isempty(findstr(hdr1(n,:),'# <Sensors count='))
    new_style = 1;
  end
end

%
% the following lines are for output from SBE data processing <=7.19
%
if new_style==0
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
      if length(ind)>2
        hdr.t_sec_sn = str2num(hdr1(n,ind(2)+1:ind(3)-1));
      else
        hdr.t_sec_sn = 0;
      end
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
      if length(ind)>2
        hdr.o_pri_sn = str2num(hdr1(n,ind(3)+1:ind(4)-1));
      else
        hdr.o_pri_sn = 0;
      end
    elseif ~isempty(findstr(hdr1(n,:),'Oxygen, SBE, secondary'))
      ind = findstr(hdr1(n,:),',');
      if length(ind)>2
        hdr.o_sec_sn = str2num(hdr1(n,ind(3)+1:ind(4)-1));
      else
        hdr.o_sec_sn = 0;
      end
    elseif ~isempty(findstr(hdr1(n,:),'Oxygen, SBE,'))
      ind = findstr(hdr1(n,:),',');
      hdr.o_pri_sn = str2num(hdr1(n,ind(2)+1:ind(3)-1));
      hdr.o_sec_sn = 0;
    end
  end
 end
end


% the following lines are for output from SBE data processing >=7.18
if new_style==1
 for n=1:size(hdr1,1)
  if strcmp(hdr1(n,6:19),'sensor Channel')
    if ~isempty(findstr(hdr1(n+1,:),'Pressure, Digiquartz')) 
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.p_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Pressure, Strain Gauge'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.p_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Temperature --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.t_pri_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Temperature, 2 --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.t_sec_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Conductivity --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.c_pri_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Conductivity, 2 --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.c_sec_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Oxygen, SBE 43 --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.o_pri_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Oxygen, SBE 43, 2 --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.o_sec_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
      if isempty(hdr.o_sec_sn)
        hdr.o_sec_sn = 0;
      end
    elseif ~isempty(findstr(hdr1(n+1,:),'PAR/Irradiance, Biospherical/Licor --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.par_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'Fluorometer, WET Labs ECO-AFL/FL --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.chl2_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    elseif ~isempty(findstr(hdr1(n+1,:),'SPAR/Surface Irradiance --'))
      ind1 = findstr(hdr1(n+3,:),'>');
      ind2 = findstr(hdr1(n+3,:),'<');
      hdr.spar_sn = str2num(hdr1(n+3,ind1(1)+1:ind2(2)-1));
    end
  end
 end
end


%
% apparently this can be non-existent in btl files
%
if ~isempty(hdr.bad_flag)
  hdr.bad_flag = -9.990e-29;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = extract_string(str,strs)
data = [];
ind = strmatch(str,strs);
if ~isempty(ind)
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
else
  data = [];
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
if isfield(hdr,'nmea_utc_str') & length(hdr.nmea_utc_str)>=20
  dat = hdr.nmea_utc_str;
  dummy = datevec(dat,'mmm dd yyyy HH:MM:SS');
  hdr.nmea_utc_datenum = datenum(dummy);
  hdr.nmea_utc = dummy;
else
  hdr.nmea_utc_datenum = [];
  hdr.nmea_utc = [];
end
