function save6( varargin )
%SAVE6 - save command forcing Matlab version 6 .mat file
%
% When working w/ Matlab 6 and 7 special commands are required in
% Matlab
% 7 to make the .mat files readable in version 6. This function
% simplifies the process. Nothing fancy is needed for load version 6
% .mat files in version 7.
%
% % Example
% a = rand(10);
% b = rand(20);
%
% save6 testfile a b
%
% % now load testfile in version 6 Matlab
%
%
% Idea from Matlab news group post
%
% http://groups.google.com/group/comp.soft-sys.matlab/browse_frm/thread/5
% b399f102f9e65a9/22643309829877ac#22643309829877ac
%
% See also SAVE, LOAD.

%% Rob Evans 2005/12/21
%
% version 1.2
% birth
% swap path's / \ acc to system		GK, Sep 2007	1.0-->1.1
% changed from ver to version		GK, Sep 2007	1.1-->1.2

%
% change / \ according to unix/mac or win systems
%
tmp = varargin{1};
if ispc==1
  ind = findstr(tmp,'/');
  if ~isempty(ind)
    tmp(ind) = '\';
  end
else
  ind = findstr(tmp,'\');
  if ~isempty(ind)
    tmp(ind) = '/';
  end
end

%% version testing
%v = ver;
%v = v(1).Version;
% changed this
v = version;


% handel version numbers like 6.5.1
%idx=find(v=='.');
%v=v(1:idx(1)-1);
%v=str2num(v);
% changed this GK
v = str2num(v(1));
if v<4
  v2 = ver;
  save version_info v2 v
  warning('please send the file version_info.mat to gkrahmann@ifm-geomar.de')
end

% just incase this is called from version 6
if ( v == 6 )
  endcmd = ';';
elseif ( v >= 7 )
  endcmd = ' -v6;';
else
  error('This function requires Matlab 6 or greater');
end

%% make the input string
for c = 2 : nargin
  tmp = [tmp ' ' varargin{c}];
end

%% the real work is done here
disp(['save ' tmp endcmd]);
evalin('caller',['save ' tmp endcmd]);
