function [] = clear_prep(stn_str)
% function [] = clear_prep(stn)
%
% clear the MAT files from prepare_cast
%
% input  :  stn     - station number
%
% version 0.1   last change 06.03.2006

% G.Krahmann, IFM-GEOMAR, March 2006

disp(' ')
disp(['clear_prep : removing mat-files for station ',stn_str])

delete(['data/ctdprof/ctdprof',stn_str,'.mat'])
delete(['data/ctdtime/ctdtime',stn_str,'.mat'])
delete(['data/nav/nav',stn_str,'.mat'])
delete(['data/sadcp/sadcp',stn_str,'.mat'])
