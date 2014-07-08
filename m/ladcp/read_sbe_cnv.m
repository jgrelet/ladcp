function [hdr,data,hdr1,data1] = read_sbe_cnv(fname,start_scan,force_pump);
% function [hdr,data,hdr1,data1] = read_sbe_cnv(fname,[start_scan],[force_pump]);
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
%           force_pump  [0]  - use data even, if the pump status is 0
%
% output :  hdr              - structure containing the info of the file
%           data             - structure containing the data of the file
%           hdr1             - character array of the file header (like hdrload)
%           data1            - data array (like hdrload)
%   
% uses :	hdrload.m, nans.m, parse_sbe_header.m
%
% version 0.24   last change 25.09.2011

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
% new oxygen voltage vars               GK, 10.02.2009  0.11-->0.12
% mumol/kg now with in situ density     GK, 04.03.2009  0.12-->0.13
% read only for p-calibration           GK, 26.10.2009  0.13-->0.14
% add gross outlier removal             GK, 29.01.2010  0.14-->0.15
% SBE dataprocessing >=7.20 has different header for sensor information,
% handling of no NMEA data              GK, 09.06.2010  0.15-->0.16
% missing pumps allowed                 GK, 14.06.2010  0.16-->0.17
% catch 98.9762 which appears to be bad GK, 10.09.2010  0.17-->0.18
% SBE doy leap year trouble fixed       GK, 30.09.2010  0.18-->0.19
% variable voltage channels (Rinko O2)  GK, 17.10.2010  0.19-->0.20
% handle data without NMEA lat lon      GK, 16.12.2010  0.20-->0.21
% handle all-zero oxygen (no sensor?)   GK, 28.04.2011  0.21-->0.22
% changed missing s_pri handling        GK, 25.08.2011  0.22-->0.23a
% introduced force_pump                 GK, 25.09.2011  0.23a-->0.24

%
% give help
%
if nargin==0
  help read_sbe_cnv
  return
end


%
% parse input arguments
%
if nargin>1
  if isstr(start_scan)
    if strcmp(start_scan,'on')
      hyst = 'on';
    end
    start_scan = 1;
  end
else
  start_scan = 1;
  hyst = '';
end
if nargin<3
  force_pump = 0;
end


%
% load data and parse header
%
[hdr1,data1] = hdrload(fname);
for n=1:size(hdr1,1)
  ind = findstr(hdr1(n,:),'\');
  if ~isempty(ind)
    hdr1(n,ind) = ' ';
  end
end
[hdr] = parse_sbe_header(fname,hdr1);


%
% cut to data starting at start_scan
%
data1 = data1(start_scan:end,:);


%
% replace flagged values by NaN
%
data1 = nans(data1,nan,hdr.bad_flag,'==');


%
% replace bad values by NaN ?
% these number appears to be bad data, but is not flagged
%
data1 = nans(data1,nan,98.9762,'==');


%
% loop over first 30 entries in data description, as these contain the names and columns of the variables
%
mll_pri = 0;
mll_sec = 0;
for n=0:30


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
    elseif ~isempty(findstr(fld,'sbeox0V'))
        data.o_v_pri = data1(:,n+1);
    elseif ~isempty(findstr(fld,'sbeox1V'))
        data.o_v_sec = data1(:,n+1);
    elseif ~isempty(findstr(fld,'sbeox0dO'))
        data.o_dvdt_pri = data1(:,n+1);
    elseif ~isempty(findstr(fld,'sbeox1dO'))
        data.o_dvdt_sec = data1(:,n+1);
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
    elseif ~isempty(findstr(fld,'v1: Voltage 1'))
        data.voltage1 = data1(:,n+1);
    elseif ~isempty(findstr(fld,'v2: Voltage 2'))
        data.voltage2 = data1(:,n+1);
    elseif ~isempty(findstr(fld,'v3: Voltage 3'))
        data.voltage3 = data1(:,n+1);
    elseif ~isempty(findstr(fld,'v4: Voltage 4'))
        data.voltage4 = data1(:,n+1);
    elseif ~isempty(findstr(fld,'v5: Voltage 5'))
        data.voltage5 = data1(:,n+1);
    elseif ~isempty(findstr(fld,'v6: Voltage 6'))
        data.voltage6 = data1(:,n+1);
    elseif ~isempty(findstr(fld,'v7: Voltage 7'))
        data.voltage7 = data1(:,n+1);
    end
end


%
% create dummy data, if only single sensors have been used
%
if ~isfield(data,'t_pri')
  data.t_pri = nan*data.p;
  data.t_pri_is_dummy = 1;
end
if ~isfield(data,'t_sec')
  data.t_sec = data.t_pri;
  data.t_sec_is_dummy = 1;
end
if ~isfield(data,'c_pri')
  data.c_pri = nan*data.p;
  data.c_pri_is_dummy = 1;
  data.s_pri = data.c_pri;
end
if ~isfield(data,'c_sec')
  data.c_sec = data.c_pri;
  data.c_sec_is_dummy = 1;
  data.s_sec = data.s_pri;
end
if ~isfield(data,'s_pri')
  data.s_pri_is_dummy = 1;
  data.s_pri = nan*data.p;
end
if ~isfield(data,'s_sec')
  data.s_sec_is_dummy = 1;
  data.s_sec = nan*data.p;
end
if ~isfield(data,'o_pri')
  data.o_pri = nan*data.p;
  data.o_pri_is_dummy = 1;
end
if ~isfield(data,'o_sec')
  data.o_sec = data.o_pri;
  data.o_sec_is_dummy = 1;
end
if ~isfield(data,'haardtc')
  data.haardtc = nan*data.p;
  data.haardtc_is_dummy = 1;
end
if ~isfield(data,'pumps')
  data.pumps = 0*data.p+1;
  data.pumps_is_dummy = 1;
end
if ~isfield(data,'latitude')
  data.latitude = nan*data.p;
  data.longitude = nan*data.p;
end

%
% check for all-zero oxygen data. This seems to be missing data (M77-2)
%
if isfield(data,'o_pri')
  if all(data.o_pri==0)
    data.o_pri = nan*data.o_pri;
    disp('found all-zero primary oxygen, replacing with NaN')
  end
end
if isfield(data,'o_sec')
  if all(data.o_sec==0)
    data.o_sec = nan*data.o_sec;
    disp('found all-zero secondary oxygen, replacing with NaN')
  end
end


%
% correct oxygen hysteresis
% convert to oxygen to mumol/kg from ml/l
%
if mll_pri==1
  sig = sw_dens(data.s_pri,data.t_pri,data.p);
  data.o_pri = data.o_pri*44600./sig;
else   % correct Seabird mumol/kg to Kiel mumol/kg
  sth = sw_pden(data.s_pri,data.t_pri,data.p,0);
  sig = sw_dens(data.s_pri,data.t_pri,data.p);
  data.o_pri = data.o_pri*44600.*sth./sig;
end
if mll_pri==1
  sig = sw_dens(data.s_sec,data.t_sec,data.p);
  data.o_sec = data.o_sec*44600./sig;
else   % correct Seabird mumol/kg to Kiel mumol/kg
  sth = sw_pden(data.s_sec,data.t_sec,data.p,0);
  sig = sw_dens(data.s_sec,data.t_sec,data.p);
  data.o_sec = data.o_sec*44600.*sth./sig;
end


%
% set all values to NaN when the pump is off
%
if isfield(data,'pumps') & force_pump==0
  warning('Setting data to NaN when pump is off, is disabled !')
  bad = find(data.pumps==0);
  data.t_pri(bad) = nan;
  data.t_sec(bad) = nan;
  data.c_pri(bad) = nan;
  data.c_sec(bad) = nan;
  data.o_pri(bad) = nan;
  data.o_sec(bad) = nan;
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
  data.datenum = data.timej-hdr.system_upload_sbe_doy + datenum(hdr.nmea_utc(1:3));
  if abs(d1-d2)>1/1440
    fprintf(1,'\n',[]);
    disp(datestr(d1,31))
    disp(datestr(d2,31))
    warning('NMEA and CTD computer time differ by more than 1 minute')
    data.timej = data.timej-hdr.system_upload_sbe_doy + datenum([0,hdr.nmea_utc(2:3)]);
  end
  jump = find(diff(data.datenum(:)')>1e-5);
  lim = [[1,jump];[jump+1,length(data.datenum)]];
  for n=1:size(lim,2)
    data.datenum(lim(1,n):lim(2,n)) = linspace(data.datenum(lim(1,n)),...
      data.datenum(lim(2,n)),lim(2,n)-lim(1,n)+1);
  end
else
  warning('no NMEA data in header, proper times are not ensured!')
  upload = datevec(hdr.system_upload_datenum);
  data.datenum = data.timej + datenum(upload.*[1,0,0,0,0,0]);
  jump = find(diff(data.datenum(:)')>1e-5);
  lim = [[1,jump];[jump+1,length(data.datenum)]];
  for n=1:size(lim,2)
    data.datenum(lim(1,n):lim(2,n)) = linspace(data.datenum(lim(1,n)),...
      data.datenum(lim(2,n)),lim(2,n)-lim(1,n)+1);
  end
end


%
% outlier removal
%
[data.p,ind] = nans(data.p,nan,-10,'<');
if ~isempty(ind)
  disp([fname,' : ','pressure has ',int2str(length(ind)),' gross negative outliers'])
  data.p_neg_outliers = ind;
else 
  data.p_neg_outliers = [];
end
[data.p,ind] = nans(data.p,nan,7000,'>');
if ~isempty(ind)
  disp([fname,' : ','pressure has ',int2str(length(ind)),' gross positive outliers'])
  data.p_pos_outliers = ind;
else 
  data.p_pos_outliers = [];
end
[data.t_pri,ind] = nans(data.t_pri,nan,-3,'<');
if ~isempty(ind)
  disp([fname,' : ','primary t has ',int2str(length(ind)),' gross negative outliers'])
  data.t_pri_neg_outliers = ind;
else 
  data.t_pri_neg_outliers = [];
end
[data.t_pri,ind] = nans(data.t_pri,nan,35,'>');
if ~isempty(ind)
  disp([fname,' : ','primary t has ',int2str(length(ind)),' gross positive outliers'])
  data.t_pri_pos_outliers = ind;
else 
  data.t_pri_pos_outliers = [];
end
[data.t_sec,ind] = nans(data.t_sec,nan,-3,'<');
if ~isempty(ind)
  disp([fname,' : ','secondary t has ',int2str(length(ind)),' gross negative outliers'])
  data.t_sec_neg_outliers = ind;
else 
  data.t_sec_neg_outliers = [];
end
[data.t_sec,ind] = nans(data.t_sec,nan,35,'>');
if ~isempty(ind)
  disp([fname,' : ','secondary t has ',int2str(length(ind)),' gross positive outliers'])
  data.t_sec_pos_outliers = ind;
else 
  data.t_sec_pos_outliers = [];
end
[data.c_pri,ind] = nans(data.c_pri,nan,-1,'<');
if ~isempty(ind)
  disp([fname,' : ','primary c has ',int2str(length(ind)),' gross negative outliers'])
  data.c_pri_neg_outliers = ind;
  data.s_pri(ind) = nan;
else 
  data.c_pri_neg_outliers = [];
end
[data.c_pri,ind] = nans(data.c_pri,nan,10,'>');
if ~isempty(ind)
  disp([fname,' : ','primary c has ',int2str(length(ind)),' gross positive outliers'])
  data.c_pri_pos_outliers = ind;
  data.s_pri(ind) = nan;
else 
  data.c_pri_pos_outliers = [];
end
[data.c_sec,ind] = nans(data.c_sec,nan,-1,'<');
if ~isempty(ind)
  disp([fname,' : ','secondary c has ',int2str(length(ind)),' gross negative outliers'])
  data.c_sec_neg_outliers = ind;
  data.s_sec(ind) = nan;
else 
  data.c_sec_neg_outliers = [];
end
[data.c_sec,ind] = nans(data.c_sec,nan,10,'>');
if ~isempty(ind)
  disp([fname,' : ','secondary c has ',int2str(length(ind)),' gross positive outliers'])
  data.c_sec_pos_outliers = ind;
  data.s_sec(ind) = nan;
else 
  data.c_sec_pos_outliers = [];
end
[data.o_pri,ind] = nans(data.o_pri,nan,-10,'<');
if ~isempty(ind)
  disp([fname,' : ','primary o has ',int2str(length(ind)),' gross negative outliers'])
  data.o_pri_neg_outliers = ind;
else 
  data.o_pri_neg_outliers = [];
end
[data.o_pri,ind] = nans(data.o_pri,nan,500,'>');
if ~isempty(ind)
  disp([fname,' : ','primary o has ',int2str(length(ind)),' gross positive outliers'])
  data.o_pri_pos_outliers = ind;
else 
  data.o_pri_pos_outliers = [];
end
[data.o_sec,ind] = nans(data.o_sec,nan,-10,'<');
if ~isempty(ind)
  disp([fname,' : ','secondary o has ',int2str(length(ind)),' gross negative outliers'])
  data.o_sec_neg_outliers = ind;
else 
  data.o_sec_neg_outliers = [];
end
[data.o_sec,ind] = nans(data.o_sec,nan,500,'>');
if ~isempty(ind)
  disp([fname,' : ','secondary o has ',int2str(length(ind)),' gross positive outliers'])
  data.o_sec_pos_outliers = ind;
else 
  data.o_sec_pos_outliers = [];
end
[data.haardtc,ind] = nans(data.haardtc,nan,-1,'<');
if ~isempty(ind)
  disp([fname,' : ','chlorophyll has ',int2str(length(ind)),' gross negative outliers'])
  data.haardtc_neg_outliers = ind;
else 
  data.haardtc_neg_outliers = [];
end
[data.haardtc,ind] = nans(data.haardtc,nan,50,'>');
if ~isempty(ind)
  disp([fname,' : ','chlorophyll has ',int2str(length(ind)),' gross positive outliers'])
  data.haardtc_pos_outliers = ind;
else 
  data.haardtc_pos_outliers = [];
end