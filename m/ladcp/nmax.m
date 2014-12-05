function [y,i] = nmax(x,dim)
%NMAX Column minimum with missing data.
%  Y = NMAX(X) returns the largest element in each column of X.
%  Missing data values must be encoded as NaNs. For vectors, NMAX(X)
%  returns the largest value of the elements in X.
%  [Y,I] = NMAX(X) stores the indices of the maximum values in vector I.
%  NMAX(X,DIM) returns maxima for dimension DIM.
%  If DIM is omitted, the first dimension of size larger than 1 is taken.
%
%  See also NMEAN, NMIN.

%  C. Mertens, IfM Kiel
%  added compatibility to MATLAB 5	G.Krahmann, LODYC Paris	

if nargin<2
  [y,i] = max(x);
else
  [y,i] = max(x,[],dim);
end
