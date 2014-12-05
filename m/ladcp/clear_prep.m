function [] = clear_prep(stn,ndigits)
% function [] = clear_prep(stn)
%
% clear the MAT files from prepare_cast
%
% input  :  stn     - station number
%
% version 0.1   last change 06.03.2006

% G.Krahmann, IFM-GEOMAR, March 2006

if nargin<2
   ndigits = 3;
end

if ischar(stn)
   stn_str = stn;
else
   stn_str = int2str0(stn,ndigits);
end

disp(' ')
disp(['clear_prep : removing mat-files for station ',stn_str])

delete(['data/ctdprof/ctdprof',stn_str,'.mat'])
delete(['data/ctdtime/ctdtime',stn_str,'.mat'])
delete(['data/nav/nav',stn_str,'.mat'])
delete(['data/sadcp/sadcp',stn_str,'.mat'])
