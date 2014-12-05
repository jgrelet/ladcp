function [y,i] = nmin(x,dim)
%NMIN Column minimum with missing data.
%  Y = NMIN(X) returns the largest element in each column of X.
%  Missing data values must be encoded as NaNs. For vectors, NMIN(X)
%  returns the smallest value of the elements in X.
%  [Y,I] = NMIN(X) stores the indices of the minimum values in vector I.
%  NMIN(X,DIM) returns minima for dimension DIM.
%  If DIM is omitted, the first dimension of size larger than 1 is taken.
%
%  uses :	intvers.m
%
%  See also NMEAN, NMAX.

%  C. Mertens, IfM Kiel
%  added compatibility to MATLAB 5	G.Krahmann, LODYC Paris	


  if nargin<2
    [y,i] = min(x);
  else
    [y,i] = min(x,[],dim);
  end
