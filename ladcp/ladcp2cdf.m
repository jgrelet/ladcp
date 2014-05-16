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

%% Check Netcdf library version
% -----------------------------
switch ( version('-release') )
  case { '11', '12', '13', '14', '2006a', '2006b', '2007a', '2007b', '2008a' }
    write_netcdf_toolbox;
  otherwise
    write_netcdf_native;
end

% use native toolbox +netcdf since R2008b
% ---------------------------------------
  function write_netcdf_toolbox
    
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
  end

% use native toolbox +netcdf since R2008b
% ---------------------------------------
  function write_netcdf_native
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
      z_dimID = netcdf.defDim(nc,'z',fix(lz));
    end
    if ltim>0
      tim_dimID = netcdf.defDim(nc,'tim',fix(ltim));
    end
    if lbot>0
      zbot_dimID = netcdf.defDim(nc,'zbot',fix(lbot));
    end
    if lsadcp>0
      z_sadcp_dimID = netcdf.defDim(nc,'zbot',fix(lbot));
    end
    lat_dimID = netcdf.defDim(nc,'lat',1);
    lon_dimID = netcdf.defDim(nc,'lon',1);
    date_dimID = netcdf.defDim(nc,'date',6);
    name_dimID = netcdf.defDim(nc,'name',length(getfield(struct,'name')));
    
    % define standard variables
    lat_varID = netcdf.defVar(nc,'lat','float',lat_dimID);
    lon_varID = netcdf.defVar(nc,'lon','float',lon_dimID);
    date_varID = netcdf.defVar(nc,'date','int',date_dimID);
    name_varID = netcdf.defVar(nc,'name','char',name_dimID);
    
    netcdf.endDef(nc);
    
    % store standard variables
    netcdf.putVar(nc,lat_varID,struct.lat);
    netcdf.putVar(nc,lon_varID,struct.lon);
    netcdf.putVar(nc,date_varID,struct.date);
    netcdf.putVar(nc,name_varID,struct.name);
    
    % parse fieldnames, define the proper variable and store it
    fnames = fieldnames(struct);
    for n=1:size(fnames,1)
      dummy = getfield(struct,fnames{n});
      if length(dummy)==lz
        netcdf.reDef(nc);
        varID = netcdf.defVar(nc,fnames{n},'float',z_dimID);
        netcdf.endDef(nc);
        netcdf.putVar(nc,varID,dummy);
        eval([fnames{n},'_varID = varID;']);
      end
      if length(dummy)==ltim
        netcdf.reDef(nc);
        varID = netcdf.defVar(nc,fnames{n},'float',tim_dimID);
        netcdf.endDef(nc);
        netcdf.putVar(nc,varID,dummy);
        eval([fnames{n},'_varID = varID;']);
      end
      if length(dummy)==lbot
        netcdf.reDef(nc);
        varID = netcdf.defVar(nc,fnames{n},'float',zbot_dimID);
        netcdf.endDef(nc);
        netcdf.putVar(nc,varID,dummy);
        eval([fnames{n},'_varID = varID;']);
      end
      if length(dummy)==lsadcp
        netcdf.reDef(nc);
        varID = netcdf.defVar(nc,fnames{n},'float',z_sadcp_dimID);
        netcdf.endDef(nc);
        netcdf.putVar(nc,varID,dummy);
        eval([fnames{n},'_varID = varID;']);
      end
    end
    
    % parse fieldnames and add them to the netcdf file
    netcdf.reDef(nc);
    
    for i=2:(nargin-1)
      eval(['a=st',int2str(i),';']);
      fnames = fieldnames(a);
      if isstruct(a)
        if ~isstruct(eval(['a.' fnames{1}])) % No SubStructure
          dummy='New Structure';
          netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),['st',int2str(i)],dummy);
          for n = 1:size(fnames,1)
            dummy = getfield(a,fnames{n});
            if size(dummy,1)==1
              if islogical(dummy)
                netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),fnames{n},int32(dummy));
              else
                netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),fnames{n},dummy);
              end
            end
          end
        else % SubStructures -> Variable Attributes
          [ndims,nvars] = netcdf.inq(nc);
          for n = 1:size(fnames,1)
            for varid=0:nvars-1
              if strcmp(fnames{n},netcdf.inqVar(nc,varid))
                atts = eval(['fieldnames(a.' fnames{n} ');']);
                for j = 1:size(atts,1)
                  dummy = eval(['a.' fnames{n} '.' atts{j} ';']);
                  if size(dummy,1) == 1
                    netcdf.putAtt(nc,varid,atts{j},dummy);
                  end
                end
              end
            end
          end
        end
      else
        disp(' not structure')
      end
    end
    netcdf.endDef(nc);
    
    % close netcdf file
    netcdf.close(nc);
  end

end
