function [arr,indi]=nans(arr,newdummy,olddummy,mode);
% function [arr,ind]=nans(arr,newdummy,olddummy,mode);
%
% replace olddummy by newdummy
%
% input  :	arr		: data array
%		newdummy	: new dummy  [default nan]
%		olddummy	: old dummy  [default -99.0]
%				  if longer than 1 search for
%				  each of the values, --> mode=1
%		mode		: optional  1 : search for == olddummy
%					    2 : search for <= olddummy
%					    3 : search for >= olddummy
%				  [default 2] may be relational string
%
% output :	arr		: data array
%		ind		: indices of replaced values
%
% version 1.0.2		last change 18.08.2000

% G.Krahmann, IfM Kiel  3. Aug 1994
% added search for more than 1 olddummy   G.Krahmann, 28.11.1994
% check for empty input		G.K. Aug 2000		1.0.1-->1.0.2

if isempty(arr)
  return
end
if nargin<2
  newdummy=nan;
end
if nargin<3
  olddummy=-99;
end
if nargin<4
  mode=2;
end
if isstr(mode)
  cstr=mode;
else
  if length(olddummy)>1
    mode=1;
  end
  if mode==1
    cstr='==';
  elseif mode==2
    cstr='<=';
  elseif mode==3
    cstr='>=';
  end
end
indi=[];
for n=1:length(olddummy)
  if ~isnan(olddummy(n))
    eval(['ind=find(arr',cstr,'olddummy(n));'])
  else
    ind=find(isnan(arr));
  end
  if ~isempty(ind)
    arr(ind)=newdummy*ones(1,length(ind));
  end
  indi=[indi;ind];
end
