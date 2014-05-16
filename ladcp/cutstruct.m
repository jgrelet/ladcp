function a=cutstruct(a,ii)
% function a=cutstruct(a,ii)
%
% reduce array size in structure
lz = length(ii);
iok = find(ii==1);
if isfield(a,'cutindx')
  a.cutindx = a.cutindx(1)-1+[iok(1) iok(end)];
else
  a.cutindx = [iok(1) iok(end)];
end
if isstruct(a)
  fnames = fieldnames(a);
  for n=1:size(fnames,1)
    dummy = getfield(a,fnames{n});
    [ly,lx]=size(dummy);
    if ly==lz
      a=setfield(a,fnames{n},dummy(iok,:));
    elseif lx==lz
      a=setfield(a,fnames{n},dummy(:,iok));
    end
  end
else
  error('first argument must be a structure')
end
