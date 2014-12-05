function [t] = num3str(x,n,m,b)
% function [t] = num3str(x,n,m,b)
%   NUM2STR Number to string conversion.
%
%   Input:     x           : real number        
%              n           : digits to the left of the decimal point
%              m           : digits to the right of the decimal point
%              b           : string of padding character, eg '0'
%
%	   T = NUM2STR(X)  converts the scalar number  X into a string
%  	   representation  T  with about  4  digits and an exponent if
%	   required.   This is useful for labeling plots with the
%	   TITLE, XLABEL, YLABEL, and TEXT commands.  See also INT2STR,
%	   SPRINTF, and FPRINTF.

%   Felix Tubiana

if nargin<1
   help num3str
   return
end
if nargin<2, n=4; end
if nargin<3, m=4; end
if nargin<4 b=' '; end
if ndims(x) == 2 & length(x) > 1
   [c,d] = size(x);
   for i = 1:c
      for j = 1:d         
         t{i,j} = num3str(x(i,j), n, m, b);
      end
   end
else
   if isstr(x)
      t = x;
   else
      t = sprintf(['%',int2str(n),'.',int2str(m),'f'],x);
      ii=find(t==' ');
      t(ii)=b;
   end
end
