function [] = prepsadcp(stn,values)
% function [] = prepsadcp(stn,values)
%
% prepare Ship-ADCP data
%
% this file is cruise specific for NBP 04-08
% 
% to create a file for your own cruise, copy this file and
% modify where necessary
%
% wf = 1	NB150
% wf = 2	OS38BB
% wf = 3	OS38NB
%
% all three files are there but only 'wf' is loaded into
% the processing routine
wf = 2;

% copy raw SADCP files
% this is data mounted via SMBMOUNT from //peale/adcp_home
% there are three different versions of SADCP data on NBP
% all of them are archived here into one file
% the one selected above with 'wf' has the standard name
% used by the LADCP processing.
% a small program could be written to copy one of the other
% data sets into the standard names
if exist('/nbp/peale/www')
  !cp /nbp/peale/www/data/* data/raw_sadcp

  % load this data and convert to standard format
  % updated regulary (by processing LADCP casts)
  % this will also be a backup of the whole
  % cruise SADCP data set (at least until the time of the last
  % processing)
  load data/raw_sadcp/311/nb150_cont_uv
  load data/raw_sadcp/311/nb150_cont_xy
  oldxyt = xyt;
  olduv = uv;
  load data/raw_sadcp/nb150_cont_uv
  load data/raw_sadcp/nb150_cont_xy
  xyt = [oldxyt,xyt];
  uv = [olduv,uv];
  lon_sadcp = xyt(1,:);
  lat_sadcp = xyt(2,:);
  tim_sadcp = xyt(3,:) + julian([2004 1 1 0 0 0]);
  u_sadcp = uv(:,1:2:end);
  v_sadcp = uv(:,2:2:end);
  z_sadcp = zc;
  save data/sadcp/sadcp_nb150

  load data/raw_sadcp/311/os38nb_cont_uv
  load data/raw_sadcp/311/os38nb_cont_xy
  oldxyt = xyt;
  olduv = uv;
  load data/raw_sadcp/os38nb_cont_uv
  load data/raw_sadcp/os38nb_cont_xy
  xyt = [oldxyt,xyt];
  uv = [olduv,uv];
  lon_sadcp = xyt(1,:);
  lat_sadcp = xyt(2,:);
  tim_sadcp = xyt(3,:) + julian([2004 1 1 0 0 0]);
  u_sadcp = uv(:,1:2:end);
  v_sadcp = uv(:,2:2:end);
  z_sadcp = zc;
  save data/sadcp/sadcp_nb38

  load data/raw_sadcp/311/os38bb_cont_uv
  load data/raw_sadcp/311/os38bb_cont_xy
  oldxyt = xyt;
  olduv = uv;
  load data/raw_sadcp/os38bb_cont_uv
  load data/raw_sadcp/os38bb_cont_xy
  xyt = [oldxyt,xyt];
  uv = [olduv,uv];
  lon_sadcp = xyt(1,:);
  lat_sadcp = xyt(2,:);
  tim_sadcp = xyt(3,:) + julian([2004 1 1 0 0 0]);
  u_sadcp = uv(:,1:2:end);
  v_sadcp = uv(:,2:2:end);
  z_sadcp = zc;
  save data/sadcp/sadcp_bb38
end

% load all data and restrict to cast
load data/sadcp/sadcp_nb150
good1 = find( tim_sadcp>values.start_time-0.1 & tim_sadcp<values.end_time+0.1);
nb150.tim_sadcp = tim_sadcp(good1);
nb150.lat_sadcp = lat_sadcp(good1);
nb150.lon_sadcp = lon_sadcp(good1);
nb150.u_sadcp = u_sadcp(:,good1);
nb150.v_sadcp = v_sadcp(:,good1);
nb150.z_sadcp = z_sadcp;
load data/sadcp/sadcp_nb38
good2 = find( tim_sadcp>values.start_time-0.1 & tim_sadcp<values.end_time+0.1);
nb38.tim_sadcp = tim_sadcp(good2);
nb38.lat_sadcp = lat_sadcp(good2);
nb38.lon_sadcp = lon_sadcp(good2);
nb38.u_sadcp = u_sadcp(:,good2);
nb38.v_sadcp = v_sadcp(:,good2);
nb38.z_sadcp = z_sadcp;
load data/sadcp/sadcp_bb38
good3 = find( tim_sadcp>values.start_time-0.1 & tim_sadcp<values.end_time+0.1);
bb38.tim_sadcp = tim_sadcp(good3);
bb38.lat_sadcp = lat_sadcp(good3);
bb38.lon_sadcp = lon_sadcp(good3);
bb38.u_sadcp = u_sadcp(:,good3);
bb38.v_sadcp = v_sadcp(:,good3);
bb38.z_sadcp = z_sadcp;

if length(good1)>3
z = [0:10:1200];
ll = min([length(good1),length(good2),length(good3)]);
nb150.u_sadcp = nb150.u_sadcp(:,1:ll);
nb150.v_sadcp = nb150.v_sadcp(:,1:ll);
nb150.tim_sadcp = nb150.tim_sadcp(:,1:ll);
nb150.lon_sadcp = nb150.lon_sadcp(:,1:ll);
nb150.lat_sadcp = nb150.lat_sadcp(:,1:ll);
nb38.u_sadcp = nb38.u_sadcp(:,1:ll);
nb38.v_sadcp = nb38.v_sadcp(:,1:ll);
nb38.tim_sadcp = nb38.tim_sadcp(:,1:ll);
nb38.lon_sadcp = nb38.lon_sadcp(:,1:ll);
nb38.lat_sadcp = nb38.lat_sadcp(:,1:ll);
bb38.u_sadcp = bb38.u_sadcp(:,1:ll);
bb38.v_sadcp = bb38.v_sadcp(:,1:ll);
bb38.tim_sadcp = bb38.tim_sadcp(:,1:ll);
bb38.lon_sadcp = bb38.lon_sadcp(:,1:ll);
bb38.lat_sadcp = bb38.lat_sadcp(:,1:ll);
u1 = interp1q(nb150.z_sadcp,nb150.u_sadcp,[0:10:1200]');
u2 = interp1q(bb38.z_sadcp,bb38.u_sadcp,[0:10:1200]');
u3 = interp1q(nb38.z_sadcp,nb38.u_sadcp,[0:10:1200]');
ss = size(u1);
u1 = u1(:);
u2 = u2(:);
u3 = u3(:);
u = meannan([u1,u2,u3],2);
u = reshape(u,ss);
v1 = interp1q(nb150.z_sadcp,nb150.v_sadcp,[0:10:1200]');
v2 = interp1q(bb38.z_sadcp,bb38.v_sadcp,[0:10:1200]');
v3 = interp1q(nb38.z_sadcp,nb38.v_sadcp,[0:10:1200]');
ss = size(v1);
v1 = v1(:);
v2 = v2(:);
v3 = v3(:);
v = meannan([v1,v2,v3],2);
v = reshape(v,ss);
%figure(3)
%subplot(2,1,1)
%pcolor(u)
%caxis([-0.1,0.1])
%shading flat
%subplot(2,1,2)
%pcolor(v)
%caxis([-0.1,0.1])
%shading flat
end

if wf==1
  tim_sadcp = nb150.tim_sadcp;
  lon_sadcp = nb150.lon_sadcp;
  lat_sadcp = nb150.lat_sadcp;
  u_sadcp = nb150.u_sadcp;
  v_sadcp = nb150.v_sadcp;
  z_sadcp = nb150.z_sadcp;
elseif wf==2
  tim_sadcp = nb38.tim_sadcp;
  lon_sadcp = nb38.lon_sadcp;
  lat_sadcp = nb38.lat_sadcp;
  u_sadcp = nb38.u_sadcp;
  v_sadcp = nb38.v_sadcp;
  z_sadcp = nb38.z_sadcp;
elseif wf==3
  tim_sadcp = bb38.tim_sadcp;
  lon_sadcp = bb38.lon_sadcp;
  lat_sadcp = bb38.lat_sadcp;
  u_sadcp = bb38.u_sadcp;
  v_sadcp = bb38.v_sadcp;
  z_sadcp = bb38.z_sadcp;
end

eval(['save data/sadcp/sadcp',int2str0(stn,3),' bb38 nb38 nb150 ',...
	'tim_sadcp lon_sadcp lat_sadcp u_sadcp v_sadcp z_sadcp'])
