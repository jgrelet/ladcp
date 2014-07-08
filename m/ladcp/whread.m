function varargout = rdread(fid)
%RDREAD Read RDI BB data.
%  RDREAD(FID)
%
% version 0.3  last change 10.06.2011

%  Christian Mertens, IfM Kiel

% read instrument serial number                   GK, 13.05.2011  xx-->0.2
% read NB/BB mode                                 GK, 10.06.2011  0.2-->0.3


% rewind to beginning of file and read header to get the number of bytes
% in each ensemble, the number of data types, and the address offsets
status = fseek(fid,0,'bof');
[nbytes,dtype,offset] = rdihead(fid);
% disp([' raw data has ',int2str(nbytes),'+2 bytes per ensemble'])
ntypes = length(dtype);

% get the number of ensembles from file size; each ensemble has nbytes
% plus two bytes for the checksum
status = fseek(fid,0,'eof');
m = floor(ftell(fid)/(nbytes + 2));
status = fseek(fid,0,'bof');

% number of bins is the offset difference (minus 2 bytes for the ID code)
% between velocity data and correlation magnitude devided by 4 beams
% times 2 bytes
n = (offset(4) - offset(3) - 2)/(2*4);

% data parameters
scale = [NaN,NaN,0.001,1,0.45,1,0.001];
precision = {'','','int16','uint8','uint8','uint8',''};
varid = [0,128,256,512,768,1024,1536];
badvals = scale.*[NaN,NaN,-32768,0,NaN,NaN,-32768];

% initialize output variables
for k = 1:length(varid)
  if k == 1
    varargout{k} = NaN*ones(1,1,14);
  elseif k == 2
    varargout{k} = NaN*ones(m,1,12);
  elseif (k >= 3 & k <= 6)
    varargout{k} = NaN*ones(m,n,4);
  elseif k == 7
    varargout{k} = NaN*ones(m,1,16);
  end
end

% read fixed leader data
status = fseek(fid,offset(1)+2,'bof');
varargout{1} = rdflead(fid);

icheck = 0;

for i = 1:m
  % read ensemble to verify the checksum
  status = fseek(fid,(i-1)*(nbytes+2),'bof');
  buffer = fread(fid,nbytes,'uint8');
  checksum = fread(fid,1,'uint16');

  % read ensemble if checksum is ok
  if checksum == rem(sum(buffer),2^16);
    for kk = 2:length(dtype)
      k = dtype(kk);
      % set file pointer to beginning of data
      status = fseek(fid,(i-1)*(nbytes+2)+offset(kk)+2,'bof');
      switch varid(k)
        case varid(2)
          % variable leader data
          varargout{k}(i,1,:) = rdvlead(fid);
        case varid(7)
          % bottom track data
          varargout{k}(i,1,:) = rdbtrack(fid);
        otherwise
          % velocity, correlation, echo intensity, or percent-good data
          a = fread(fid,4*n,precision{k});
          varargout{k}(i,:,:) = scale(k)*reshape(a,4,n)';
      end
    end
  else
   icheck=icheck+1;
  end
end

if icheck > m*0.01 
  disp(['>   Found ',int2str(icheck),' ensembles with bad checksum '])
end

%
% introduced the following after errors in a BB profile  GK, June 2006
%
bad = find( isnan( varargout{2}(:,1,1) ) );
if ~isempty( bad )
  good = find( ~isnan( varargout{2}(:,1,1) ) );
  disp(['>   Found ',int2str(length(bad)),' ensembles with NaN as time '])
  disp('>     removing those and continuing with the rest')
  varargout{2} = varargout{2}(good,:,:);
  varargout{3} = varargout{3}(good,:,:);
  varargout{4} = varargout{4}(good,:,:);
  varargout{5} = varargout{5}(good,:,:);
  varargout{6} = varargout{6}(good,:,:);
end

% check for bad values
i = find(varargout{3} == badvals(3));
varargout{3}(i) = NaN;
i = find(varargout{4} == badvals(4));
varargout{4}(i) = NaN;
% bottom track
i = find(varargout{7}(:,1,5) == badvals(7));
varargout{7}(i,1,1:8) = NaN;
i = find(varargout{7}(:,1,6) == badvals(7));
varargout{7}(i,1,1:8) = NaN;
i = find(varargout{7}(:,1,7) == badvals(7));
varargout{7}(i,1,1:8) = NaN;
i = find(varargout{7}(:,1,8) == badvals(7));
varargout{7}(i,1,1:8) = NaN;


%-------------------------------------------------------------------------------

function fl = rdflead(fid);
%RDFLEAD Read the fixed leader data from a raw ADCP data file.
%  FL = RDFLEAD(FID)

fseek(fid,7,'cof');

% number of depth cells
fl(1) = fread(fid,1,'uint8');

% pings per ensemble, depth cell length in cm, blank after transmit
fl(2:4) = fread(fid,3,'uint16');
fl(3) = 0.01*fl(3);
fl(4) = 0.01*fl(4);

fseek(fid,16,'cof');

% Bin 1 distance, xmit pulse length
fl(5:6) = 0.01*fread(fid,2,'ushort');

fseek(fid,6,'cof');

% Serial Number of CPU board
fl(7:14) = fread(fid,8,'uint8');

% Bandwith 1=NB  0=BB
fl(15) = fread(fid,1,'uint16');

fseek(fid,2,'cof');

% Serial Number of Instrument
fl(16:19) = fread(fid,4,'uchar');


%-------------------------------------------------------------------------------

function vl = rdvlead(fid)
%RDVLEAD Read the variable leader data from a raw ADCP data file.
%  VL = RDVLEAD(FID)

fseek(fid,2,'cof');

% time of ensemble
c = fread(fid,7,'uint8');

% fix y2k problem because of RDI only storing 2 digits for the year
if c(1)<80
  c(1) = c(1)+2000;
else
  c(1) = c(1)+1900;
end
vl(1) = julian(c(1),c(2),c(3),c(4)+c(5)/60+c(6)/3600+c(7)/360000);

fseek(fid,3,'cof');

% speed of sound (EC)
vl(7) = fread(fid,1,'uint16');
fseek(fid,2,'cof');

% heading (EH)
vl(4) = 0.01*fread(fid,1,'uint16');

% pitch (EP) and roll (ER)
vl(2:3) = 0.01*fread(fid,2,'int16');

% salinity (ES)
vl(6) = 0.001*fread(fid,1,'uint16');

% temperature (ET)
vl(5) = 0.01*fread(fid,1,'int16');
fseek(fid,6,'cof');

% ADC channels
% Transmit Current
vl(8) = fread(fid,1,'uint8');
% Transmit Volt
vl(9) = fread(fid,1,'uint8');
% Internal Temperature
vl(10) = fread(fid,1,'uint8');

fseek(fid,11,'cof');

% Internal Pressure
vl(11) = fread(fid,1,'uint32')/1000;
% Internal Pressure Variance
vl(12) = fread(fid,1,'uint32')/1000;



%-------------------------------------------------------------------------------

function bt = rdbtrack(fid)
%RDBTRACK Read the bottom track data from a raw ADCP data file.
%  BT = RDBTRACK(FID)

fseek(fid,14,'cof');

% range
bt(1:4) = 0.01*fread(fid,4,'uint16');

% velocity
bt(5:8) = 0.001*fread(fid,4,'int16');

% correlation magnitude and percent good
bt(9:16) = fread(fid,8,'uint8');
