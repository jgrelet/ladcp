function [err] = calc_rdi_std(info)
% function [err] = calc_rdi_std(info)
%
% calculate RDI's uncertainty estimate of a single bin
%
% input  : info           - LADCP data.down or data.up structure

% output : err            - RDI estimate of velocity standard deviation
%
% version 0.1  last change 02.06.2011

% G.Krahmann, IFM-GEOMAR   June 2011

%
% this error estimate is deduced from values given by PlanADCP version
% 2.06
% 150kHz, NB, 1 ping, 10m,  5.72cm/s
% 150kHz, NB, 2 ping, 10m,  4.04cm/s
% 150kHz, NB, 1 ping, 20m,  2.82cm/s
% 150kHz, NB, 1 ping,  5m, 11.19cm/s
% 150kHz, NB, 1 ping,  3m, 18.29cm/s
% 150kHz, BB, 1 ping, 10m,  2.84cm/s
% 150kHz, BB, 2 ping, 10m,  2.01cm/s
% 150kHz, BB, 1 ping, 20m,  1.36cm/s
% 150kHz, BB, 1 ping,  5m,  5.64cm/s
% 150kHz, BB, 1 ping,  3m,  9.25cm/s
% 300kHz, NB, 1 ping,  8m,  1.73cm/s
% 300kHz, NB, 2 ping,  8m,  1.22cm/s
% 300kHz, NB, 1 ping,  2m,  6.94cm/s
% 300kHz, NB, 1 ping,  4m,  3.55cm/s
% 300kHz, BB, 1 ping,  8m,  3.68cm/s
% 300kHz, BB, 2 ping,  8m,  2.60cm/s
% 300kHz, BB, 1 ping,  2m, 14.51cm/s
% 300kHz, BB, 1 ping,  4m,  7.44cm/s
%
% NB 150  err = 1/sqrt(ping) *( -9.7738 + 58.2524 * bl - 0.0467 * bl^2 )  
% BB 150  NB 150 / 2
% NB 300  err = 1/sqrt(ping) *( -2.6667 + 15.5600 * bl - 0.1733 * bl^2 )  
% BB 300  NB 300 / 2.1

% for testing purposes
if ~isstruct(info)
  dummy = info;
  clear info
  info.Cell_length = dummy(1);
  info.Frequency = dummy(2);
  info.bandwidth = dummy(3);
  info.Pings_per_Ensemble = dummy(4);
end

%
bl = 1/(info.Cell_length/100);
if info.Frequency==150
  err = 1/sqrt(info.Pings_per_Ensemble) * (-9.7738 * bl^2 + 58.2524 * bl - 0.0467);
  if info.bandwidth==0
    err = err / 2;
  end
elseif info.Frequency==300
  err = 1/sqrt(info.Pings_per_Ensemble) * (-2.6667 * bl^2 + 15.5600 * bl - 0.1733);
  if info.bandwidth==0
    err = err / 2;
  end
else
  disp('>   RDI standard deviation estimate not implemented for this frequency')
  err = nan;
end
err = err/100;