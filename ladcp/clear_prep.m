function [] = clear_prep(stn)
% function [] = clear_prep(stn)
%
% clear the MAT files from prepare_cast
%
% input  :  stn     - station number
%
% version 0.1   last change 06.03.2006

% G.Krahmann, IFM-GEOMAR, March 2006

disp(' ')
disp(['clear_prep : removing mat-files for station ',int2str(stn)])

delete(['data/ctdprof/ctdprof',int2str0(stn,3),'.mat'])
delete(['data/ctdtime/ctdtime',int2str0(stn,3),'.mat'])
delete(['data/nav/nav',int2str0(stn,3),'.mat'])
delete(['data/sadcp/sadcp',int2str0(stn,3),'.mat'])
