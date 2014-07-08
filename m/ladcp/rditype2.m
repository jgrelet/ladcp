function [d] = rditype(fid);
% function [d] = rditype(fid);
%
%	Read the fixed leader data from a raw ADCP
%	data file.
%	Returns the contents of the fixed leader
%	as 31 elements of the vector 'data' or an
%	empty matrix if the fixed leader ID is not
%	identified (error condition)

% get fixed leader data

% rewind file 
status = fseek(fid,0,'bof');

% header identification ID
hid = 127;
[id,n] = fread(fid,2,'uint8');
status = fseek(fid,0,'bof');

% if not 127, it must be an old narrowband instrument
if id(1) ~= hid

  % Fix parameter for NB-ADCP Kiel
  data = zeros(1,32);
  d.nbytes = NaN;
  d.ntypes = NaN;
  d.Firm_Version = 0;
  d.Frequency = 150;
  d.Up = 0;
  d.Beam_angle = 20;
  d.Cell_length = 1600;
  d.Pings_per_Ensemble = 8;
  d.Pulse_length = 1600;
  d.Blank = NaN;
  d.Mode = NaN;
  d.Time_Pings = NaN;
  d.Coordinates = 3;
  d.NarrowBand = 1;
  return

end
[d.nbytes, d.ntypes, offsets] = rdihead(fid);

 % ---- now at start of fixed leader data

% originally 
% Written by Marinna Martini
% for the U.S. Geological Survey
% Atlantic Marine Geology, Woods Hole, MA
% 1/7/95

data = zeros(1,32);
fld = 1;
  
% make sure we're looking at the beginning of
% the fixed leader record by testing for it's ID
data(fld) = fread(fid,1,'ushort');
if(data(fld)~=0),
  disp('Fixed Leader ID not found');
  data = [];
  return;
end
fld = fld+1;

% version number of CPU firmware
data(fld) = fread(fid,1,'uchar');
fld = fld+1;

% revision number of CPU firmware
data(fld) = fread(fid,1,'uchar');
d.Firm_Version = data(fld-1)+data(fld)/100;
fld = fld+1;

% configuration, uninterpreted
data(fld) = fread(fid,1,'uchar');
b = dec2binv(data(fld));
freqs = [75 150 300 600 1200 2400];
junk = bin2decv(b(6:8));
d.Frequency = freqs(junk+1);
d.Up = str2num(b(1));
fld = fld+1;

data(fld) = fread(fid,1,'uchar');
b = dec2binv(data(fld));
angles = [15 20 30 0];
junk = bin2decv(b(7:8));
d.Beam_angle = angles(junk+1);

% real (0) or simulated (1) data flag
fld = fld+1;
data(fld) = fread(fid,1,'uchar');	

% undefined
fld = fld+1;
data(fld) = fread(fid,1,'uchar');	

% number of beams
fld=fld+1;
data(fld) = fread(fid,1,'uchar');	

% number of depth cells
fld=fld+1;
data(fld) = fread(fid,1,'uchar');
d.Depth_Cells = data(fld);

% pings per ensemble
fld = fld+1;
data(fld) = fread(fid,1,'ushort');
d.Pings_per_Ensemble = data(fld);

% depth cell length in cm
fld = fld+1;
data(fld) = fread(fid,1,'ushort');
d.Cell_length = data(fld);

% blanking distance (WF)
fld = fld+1;
data(fld) = fread(fid,1,'ushort');
d.Blank = data(fld);

% Profiling mode (WM)
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.Mode = data(fld);

% Minimum correlation threshold (WC)
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.Min_Correlation = data(fld);

% number of code repetitions
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.Code_rep = data(fld);

% Minimum percent good to output data (WG)
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.Min_Percgood = data(fld);

% Error velocity threshold (WE)
fld = fld+1;
d.Max_Errorvel =data(fld);
data(fld) = fread( fid,1,'ushort');

% time between ping groups (TP)
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.Time_Pings = data(fld)*60;

fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.Time_Pings = d.Time_Pings+data(fld);

fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.Time_Pings = d.Time_Pings+data(fld)/100;

% coordinate transformation (EX)
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
b = dec2binv(data(fld));
junk = bin2decv(b(4:5));
d.Coordinates = junk;
d.use_tilt = str2num(b(6));

% Heading Alignment (EA)
fld = fld+1;
data(fld) = fread(fid,1,'uint16');
d.headin_alignment = data(fld);

% Heading Bias (EB)
fld = fld+1;
data(fld) = fread(fid,1,'uint16');
d.headin_bias = data(fld);

% Sensor source (EZ)
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
b = dec2binv(data(fld));
d.sensor_source = b;

% Sensors available
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
b = dec2binv(data(fld));
d.sensor_avail = b;

% Bin 1 distance
fld = fld+1;
data(fld) = fread(fid,1,'ushort');

% xmit pulse length
fld = fld+1;
data(fld) = fread(fid,1,'ushort');
d.Pulse_length = data(fld);

% starting depth cell
fld = fld+1;
data(fld) = fread(fid,1,'uchar');

% ending depth cell
fld = fld+1;
data(fld) = fread(fid,1,'uchar');

% false target reject threshold
fld = fld+1;
data(fld) = fread(fid,1,'uchar');
d.target_max = data(fld);

% spare
fld = fld+1;
data(fld) = fread(fid,1,'uchar');

% transmit lag distance
fld = fld+1;
data(fld) = fread(fid,1,'ushort');
d.xmit_lag = data(fld);

% CPU board serial number
for n=1:8
  cpu(n) = fread(fid,1,'uchar');
end
d.cpu_board_serial_number = cpu;

% System bandwidth
fld = fld+1;
data(fld) = fread(fid,1,'ushort');
d.system_bandwidth = data(fld);

% System power
fld = fld+1;
data(fld) = fread(fid,1,'uint8');
d.system_power = data(fld);

% spare
fld = fld+1;
data(fld) = fread(fid,1,'uchar');

% Instrument Serial Number
for n=1:2
  serial_number(n) = fread(fid,1,'uchar');
end
d.instrument_serial_number = serial_number;


%--------------------------------------

function d = bin2decv(h)
%BIN2DEC  BIN2DEC('X') returns the decimal number corresponding to the
%        binary number in quotes.  
%		 For example, BIN2DEC([1 1 0 0]) returns 12.

n = length(h);
h = fliplr(h);
d = 0;
for i=1:n
  if isstr(h(i)),
    p(i) = str2num(h(i)).*(2^(i-1));
  else
    p(i) = h(i).*(2^(i-1));
  end
end
d = sum(p);

% ---------------------------------

function h = dec2binv(d)
%DEC2BIN DEC2BIN(d) returns the binary number corresponding to the decimal
%        number d.  For example, DEC2BIN(202) returns '11001010'.
%
%	
h = dec2bin(d);
bits = length(h)/8;
if bits ~= fix(bits),
  % not an even multiple of 8, so pad
  npad = 8-((bits-fix(bits))*8);
  h0(1:npad)='0';
  h = [h0,h];
end
