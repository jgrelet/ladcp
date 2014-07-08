function [] = ladcp2cdf(fname,struct,st2,st3,st4,st5,st6,st7);
% function [] = ladcp2cdf(fname,struct,st2,st3,st4,st5,st6,st7);
%
% function to save LADCP data into a netcdf file
%
% input  :	fname		- filename
%		struct		- structure containing the data
%		struct2		- structure containing the attributes
%
% version 0.1	last change 08.03.2002

% Felix Tubiana Added Variable Attributes 12.12.03
% G.Krahmann, LDEO Mar 2002
% warning off MATLAB:DeprecatedLogicalAPI

% check arguments
if nargin<2
  error('need two input arguments')
end
if ~isstruct(struct)
  error('second argument must be a structure')
end

% check filename
fname = deblank(fname);

nc = netcdf(fname,'clobber');

if isempty(nc)
  error('can not open netcdf file')
end


% look for the three dimension bases
lbot = 0;
lz = 0;
ltim = 0;
lsadcp = 0;
if isfield(struct,'z');
  lz = length(getfield(struct,'z'));
end  
if isfield(struct,'tim');
  ltim = length(getfield(struct,'tim'));
end  
if isfield(struct,'zbot');
  lbot = length(getfield(struct,'zbot'));
end  
if isfield(struct,'z_sadcp');
  lsadcp = length(getfield(struct,'z_sadcp'));
end  

% define dimensions in netcdf file
if lz>0
  nc('z') = fix(lz);
end
if ltim>0
  nc('tim') = fix(ltim);
end
if lbot>0
  nc('zbot') = fix(lbot);
end
if lsadcp>0
  nc('z_sadcp') = fix(lsadcp);
end
nc('lat') = 1;
nc('lon') = 1;
nc('date') = 6;
nc('name') = length(getfield(struct,'name'));

% define standard variables
nc{'lat'} = ncfloat('lat');
nc{'lon'} = ncfloat('lat');
nc{'date'} = ncint('date');
nc{'name'} = ncchar('name');

% store standard variables
nc{'lat'}(:) = struct.lat;
nc{'lon'}(:) = struct.lon;
nc{'date'}(:) = struct.date;
nc{'name'}(:) = struct.name;

% parse fieldnames, define the proper variable and store it
fnames = fieldnames(struct);
for n=1:size(fnames,1)
  dummy = getfield(struct,fnames{n});
  if length(dummy)==lz
    eval(['nc{''',fnames{n},'''} = ncfloat(''z'');'])
    eval(['nc{''',fnames{n},'''}(:) = dummy;'])
  end
  if length(dummy)==ltim
    eval(['nc{''',fnames{n},'''} = ncfloat(''tim'');'])
    eval(['nc{''',fnames{n},'''}(:) = dummy;'])
  end
  if length(dummy)==lbot
    eval(['nc{''',fnames{n},'''} = ncfloat(''zbot'');'])
    eval(['nc{''',fnames{n},'''}(:) = dummy;'])
  end
  if length(dummy)==lsadcp
    eval(['nc{''',fnames{n},'''} = ncfloat(''z_sadcp'');'])
    eval(['nc{''',fnames{n},'''}(:) = dummy;'])
  end
end
  
% parse fieldnames and add them to the netcdf file  
for i=2:(nargin-1)
   eval(['a=st',int2str(i),';']);
   fnames = fieldnames(a);
   if isstruct(a)
      if ~isstruct(eval(['a.' fnames{1}])) % No SubStructure
	 dummy='New Structure';
	 eval(['nc.',['st',int2str(i)],'=dummy;'])
	 for n = 1:size(fnames,1)
	    dummy = getfield(a,fnames{n});
	    if size(dummy,1)==1
	       if isstr(dummy)
		  eval(['nc.',fnames{n},'=dummy;'])
	       else
		  eval(['nc.',fnames{n},'=dummy(:);'])
	       end    
	    end
	 end
      else % SubStructures -> Variable Attributes
	 for n = 1:size(fnames,1)
	    atts = eval(['fieldnames(a.' fnames{n} ');']);
	    for j = 1:size(atts,1)
	       dummy = eval(['a.' fnames{n} '.' atts{j} ';']);
	       if size(dummy,1) == 1
		  if isstr(dummy)
		     eval(['nc{''',fnames{n}, '''}.' atts{j} '=dummy;'])
		  else
		     eval(['nc{''',fnames{n}, '''}.' atts{j} '=dummy(:);'])
		  end    
	       end
	    end	       
	 end
      end
   else
      disp(' not structure')
   end
end

% close netcdf file
close(nc)
