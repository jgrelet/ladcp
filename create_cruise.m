function create_cruise()
%
% LADCP processing software 10.0
%
% create a directory structure for a new cruise
%
% version 0.2   last change 4.12.2006

% G.Krahmann, Jul 2005

% windows incompatibility       Dec 2006 0.1->0.2


% check whether we are at the right place
if ~exist('m/ladcp') | ~exist('m/sw') | ~exist('m/initial_dir')
  disp('LADCP processing directories appear not to be in place.')
  disp('You need the directories  m/ladcp  m/sw  m/initial_dir')
  disp('to use the software.')
  return
end

% get a new cruise identifier
ok = 0;
while ok==0
  disp('This script will create a new directory with an')
  disp('identifier as name. All processing will take place')
  disp('within this directory and its subdirectories.')
  disp('It should not be longer than 4 characters')
  cruise_id = input('Please enter the new cruise identifier! ','s');
  if ~exist(cruise_id)
    ok = 1;
    mkdir(cruise_id)
  else
    disp(' ')
    disp('The given directory already exists.')
    disp(' ')
  end
end

% create all directories
mkdir([cruise_id,'/data'])
mkdir([cruise_id,'/data/ctdprof'])
mkdir([cruise_id,'/data/ctdtime'])
mkdir([cruise_id,'/data/nav'])
mkdir([cruise_id,'/data/sadcp'])
mkdir([cruise_id,'/data/raw_ctdtime'])
mkdir([cruise_id,'/data/raw_ctdprof'])
mkdir([cruise_id,'/data/raw_nav'])
mkdir([cruise_id,'/data/raw_sadcp'])
mkdir([cruise_id,'/data/raw_ladcp'])
mkdir([cruise_id,'/logs'])
mkdir([cruise_id,'/profiles'])
mkdir([cruise_id,'/plots'])
mkdir([cruise_id,'/tmp'])
mkdir([cruise_id,'/m'])

% copy some initial files
if ispc
  copyfile('m\initial_dir\cast_params.m',[cruise_id]);
  copyfile('m\initial_dir\cruise_params.m',[cruise_id]);
  copyfile('m\initial_dir\prepctdprof.m',[cruise_id,'\m']);
  copyfile('m\initial_dir\prepctdtime.m',[cruise_id,'\m']);
  copyfile('m\initial_dir\prepnav.m',[cruise_id,'\m']);
  copyfile('m\initial_dir\prepsadcp.m',[cruise_id,'\m']);
  copyfile('m\initial_dir\prepladcp.m',[cruise_id,'\m']);
  copyfile('m\initial_dir\batt_info.dat',[cruise_id,'\m']);
  copyfile('m\initial_dir\startup.m',cruise_id);
else
  copyfile('m/initial_dir/cast_params.m',[cruise_id]);
  copyfile('m/initial_dir/cruise_params.m',[cruise_id]);
  copyfile('m/initial_dir/prepctdprof.m',[cruise_id,'/m']);
  copyfile('m/initial_dir/prepctdtime.m',[cruise_id,'/m']);
  copyfile('m/initial_dir/prepnav.m',[cruise_id,'/m']);
  copyfile('m/initial_dir/prepsadcp.m',[cruise_id,'/m']);
  copyfile('m/initial_dir/prepladcp.m',[cruise_id,'/m']);
  copyfile('m/initial_dir/batt_info.dat',[cruise_id,'/m']);
  copyfile('m/initial_dir/startup.m',cruise_id);
end

disp(' ')
disp('Created the necessary directories')
disp('You now need to modify ALL cruise dependent m-files')
disp(['in ',cruise_id,'/m'])
disp(' ')
disp('You also need to quit matlab and start from the new directory')
disp(' ')

%
% workaround around a bug in Matlab 6.5
%
function [success, OSMessage]=copyfile(src,dest)
if ispc
  [Status, OSMessage] = dos(['copy ' src ' '  dest]);
elseif isunix
  [Status, OSMessage] = unix(['cp -r ' src ' ' dest]);
end
success = ~Status;
