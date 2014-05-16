function [ind] = findany(data,vals)
% function [ind] = findany(data,vals)
% 
% find occurences of any of the elements of vals within data
%
% input  :	data		- data matrix
%		vals		- vector of values
%
% output :	ind		- indices
%
% version 0.1 last change 22.11.2004

% G.Krahmann, LDEO Nov 2004

% loop over vals
ind = repmat(0,size(data));
for n=1:length(vals(:))
  ind = ind + (data==vals(n));
end

% handle NaNs separately
if any(isnan(vals))
  ind = ind + isnan(data);
end

% find results
ind = find(ind>0);
