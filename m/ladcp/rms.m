function y = rms(x,dim)
%RMS    root-mean square. For vectors, RMS(x) returns the standard
%       deviation.  For matrices, RMS(X) is a row vector containing
%       the root-mean-square of each column. The difference to STD is
%       that here the mean is NOT removed.  
%       RMS(X,DIM) returns the root-mean-square of dimension DIM
%
%	See also STD,COV.
%

%       Uwe Send, IfM Kiel, Apr 1992
% added NaN handling   Gerd Krahmann, IfM Kiel, Oct 1993, Jun 1994
% removed bug in NaN handling   G.Krahmann, Aug 1994
% added compatibility to MATLAB 5	G.Krahmann, LODYC Paris, Jul 1997


  if nargin<2
    dim=min(find(size(x)>1));
  end

  if all(isnan(x))
    y=nan;
    return
  end    

  x = shiftdim(x,dim-1);
  s = size(x);
  so = s(1);
  s(1) = 1;

  for n = 1:prod(s)
    good = find(~isnan(x(:,n)));
    if ~isempty(good)
      y(1,n) = norm( x(good,n) ) / sqrt(length(good));
    else
      y(1,n) = NaN;
    end
  end
  y = reshape(y,s);
  y = shiftdim( y, ndims(x)-dim+1 );
