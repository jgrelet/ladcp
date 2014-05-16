function [nbytes,dtype,offset] = rdihead(fid,rewind)
% function [nbytes,dtype,offset] = rdihead(fid,rewind)
%
% RDIHEAD Read the header data from a raw ADCP data file.

if nargin<2
  rewind = 0;
end

hid = 127;  % header identification byte
sid = 127;  % data source identification byte

% get file position pointer
fpos = ftell(fid);

% check header and data source identification bytes
[id,n] = fread(fid,2,'uint8');
if  (n < 2 | feof(fid))
  error('Unexpected end of file.')
end
if (id(1) ~= hid | id(2) ~= sid)
  error('Header identification byte not found.')
end

% read the number of bytes
nbytes = fread(fid,1,'uint16');

% skip spare byte
fseek(fid,1,'cof');

% read the number of data types
ndt = fread(fid,1,'uint8');

% read address offsets for data types
offset = fread(fid,ndt,'uint16');

% workaround between Martini and Mertens versions
pos_after_offset = ftell(fid);

% read variable identifiers
varid = [0 128 256 512 768 1024 1536];
for i = 1:ndt
  fseek(fid,fpos+offset(i),'bof');
  id = fread(fid,1,'uint16');
  dtype(i) = find(id == varid);
end

% workaround between Martini and Mertens versions
fseek(fid,pos_after_offset,'bof');

% rewind to the beginning of the ensemble
if rewind==1
  fseek(fid,fpos,'bof');
end
