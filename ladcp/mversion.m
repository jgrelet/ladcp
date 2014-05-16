function [id]=mversion()
% function [id]=mversion()
%
% matlab version number
%
% output :	id		- 61 for 6.1, 65 for 6.5 etc
%
% version 0.1	last change 02.06.2007

% G.Krahmann, IFM-GEOMAR June 2007

v = version;
id = str2num(v([1,3]));
