% script de traitement de rejeu des profils L-ADCP

% liste des stations
 st = ls('data/raw_ladcp/')
 
 for i=10: length(st)
   str = sprintf('process_cast(''%s'')', st(i, :));
   eval(str);
 end