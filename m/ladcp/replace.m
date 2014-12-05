function [dataout,ind] = replace(datain,relop,val)
% function [dataout,ind] = replace(datain,relop,val)
%
% replace elements in datain which fullfill the relational
% operation relop (could also be a matrix of 0 and 1) by the
% the value val
%
% input  :	datain		- input matrix
%		relop		- relational operation (has to be the
%				  same size as datain
%		val		- scalar that will replace matches
%
% output :	dataout		- output matrix with replaced values
%		ind		- indices of replaced values
%
% version 0.1	last change 27.7.2005

% G.Krahmann, IFM-GEOMAR

ind = find(relop);
dataout = datain;
if ~isempty(ind)
  dataout(ind) = val;
end
