function d=setdefv(d,n,v)
% function d=setdefv(d,n,v)
%
% if field (n) does not exist in structure (d) set it to (v)

% G.K. modified to use standard Matlab functions

if ~isfield(d,n)
  d = setfield(d,n,v);
end
