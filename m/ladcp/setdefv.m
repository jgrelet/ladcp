function self = setdefv(self, prop, value)
% function self = setdefv(self, prop, value)
% replace the original  G.K. setdefv function
% J. Grelet IRD May 2014
% function d=setdefv(d,n,v)
%
% if field (n) does not exist in structure (d) set it to (v)
% G.K. modified to use standard Matlab functions
if isstruct(self)
  if ~isfield(self, prop)
    self.(prop)= value;
  end
  
elseif isobject(self)
  
  % create dynamics properties and initialize to empty structure
  % ------------------------------------------------------------
  if ~isprop(self, prop)
    addprop(self, prop);
    self.(prop) = value;
  end
else
  error('ladcp:setdefv', 'invalid parameter');
end
