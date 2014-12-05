function varargout = nbread(fid)
% function varargout = nbread(fid)
%
%NBREAD
%
% input  :  fid            - file identifier
%
% output :  ?
%
% version 0.2  last change 13.07.2012

% renamed cosd and sind to cos_d and sin_d             GK, 13.07.2012  0.1-->0.2

% rewind to beginning of file and read header to get the number of bytes
% in each ensemble, the number of data types, and the address offsets
status = fseek(fid,0,'bof');
[nbytes,dtype,offset] = nbhead(fid);
ntypes = length(dtype);

% get the number of ensembles from file size; each ensemble has nbytes
% plus two bytes for the checksum
status = fseek(fid,0,'eof');
ftell(fid);
m = floor(ftell(fid)/(nbytes + 2));
status = fseek(fid,0,'bof');

% number of bins
nbin = (offset(4) - offset(3))/6;

% data parameters
varid = [1:7];
scale = [NaN,NaN,0.0025,1,1,1,1];
precision = {'','','bit12','uint8','uint8','uint8','bit4'};
bad = scale.*[NaN,NaN,NaN,NaN,NaN,NaN,NaN];

% initialize output variables
for k = 1:length(varid)
  if k == 1
    varargout{k} = NaN*ones(1,1,14);
  elseif k == 2
    varargout{k} = NaN*ones(m,1,13);
  else
    varargout{k} = NaN*ones(m,nbin,4);
  end
end

% read fixed leader
status = fseek(fid,offset(1),'bof');
varargout{1} = nbflead(fid);

for i = 1:m

  % read ensemble to verify the checksum
  status = fseek(fid,(i-1)*(nbytes+2),'bof');
  buffer = fread(fid,nbytes,'uint8');
  checksum = fread(fid,1,'uint16');

  % read ensemble if checksum is ok
  if checksum == rem(sum(buffer),65536)
    for k = dtype(2:end)
      % set file pointer to beginning of data
      status = fseek(fid,(i-1)*(nbytes+2)+offset(k),'bof');
      switch varid(k)
        case varid(2)
          % variable leader data
          varargout{k}(i,1,:) = nbvlead(fid);
        otherwise
          % velocity, spectral width, amplitude, percent-good, or status data
          a = fread(fid,4*nbin,precision{k});
          varargout{k}(i,:,:) = scale(k)*reshape(a,4,nbin)';
      end
    end
  end

end


% scale pitch, roll, and heading
varargout{2}(:,1,2:4) = varargout{2}(:,1,2:4)*360/65536;

% scale temperature
varargout{2}(:,1,5) = 45 - varargout{2}(:,1,5)*50/4096;


%-------------------------------------------------------------------------------

function [nbytes,dtype,offset] = nbhead(fid)
%NBHEAD Read header data from raw narrow-band ADCP data file.
%  [NBYTES,DTYPE,OFFSET] = nbhead(FID)

h = fread(fid,7,'uint16')';

% number of bytes
nbytes = h(1);

% address offsets and data types
varid = [1:7];
offset = [14 h(2:end)];
dtype = varid(offset ~= 0);
offset = [14 cumsum(offset(1:end-1))];


%-------------------------------------------------------------------------------

function fl = nbflead(fid)
%NBFLEAD Read fixed leader data from raw narrow-band ADCP data file.

fl = zeros(1,6);
f = cos_d(20)/cos_d(30);
fseek(fid,8,'cof');

% pings per ensemble
fl(2) = fread(fid,1,'uint16');

% bins per ping
fl(1) = fread(fid,1,'uint8');

% bin length
fl(3) = fread(fid,1,'uint8');
fl(3) = f*2^fl(3);

% transmit interval
fl(6) = fread(fid,1,'uint8');

% delay
fl(4) = f*fread(fid,1,'uint8');

% bin 1 distance
fl(5) = fl(4) + fl(3)/2;

% attenuation
fl(7) = 0.039;

% source level;
fl(8) = 100;

% serial number;
fl(9:14)=[3 4 5 6 7 8];


%-------------------------------------------------------------------------------

function vl = nbvlead(fid)
%NBVLEAD Read variable leader data from raw narrow-band ADCP data file.

vl = zeros(1,13);

% time of ensemble (mm/dd hh:mm:ss)
a = fread(fid,5,'uint8');
a = dec2hex(a);
c(1)=1900;
c(2:6) = str2num(a);
vl(1) = julian(c(1),c(2),c(3),c(4)+c(5)/60+c(6)/3600);

% pitch, roll, heading, and temperature
fseek(fid,16,'cof');
vl(2:5) = fread(fid,4,'uint16');
vl(2:3) = vl(2:3) - floor(vl(2:3)/(183*180))*360*182;

vl(6) = 35;
vl(7) = 1536;
vl(8:10)=[nan nan nan];


%-------------------------------------------------------------------------------

function bt = nbbtrack(fid)
%RDBTRACK Read the bottom track data from a raw ADCP data file.
%  BT = NBBTRACK(FID)

fseek(fid,14,'cof');

% range
%bt(1:4) = 0.01*fread(fid,4,'uint16');

% velocity
%bt(5:8) = 0.001*fread(fid,4,'int16');

% correlation magnitude and percent good
%bt(9:16) = fread(fid,8,'uint8');

% not implemented
bt(1:16) = NaN;

%-------------------------------------------------------------------------------


function d = y2k(d)
% fix date string
if d<80, 
  d=2000+d; 
end
if d<100, 
  d=1900+d; 
end
