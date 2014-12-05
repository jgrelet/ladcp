function [mred,mind,sred]=meanreduce(data,n,dim);
% function [mred,mind,sred]=meanreduce(data,n,dim);
%
% average n values (NaN are allowed)
%
% input  :	data		: data array
%		n		: how many values shall be averaged
%				  (at the end all remaining values are used)
%		[dim]		: average which dimension
%
% output :	mred		: averaged array
%		mind		: averaged indices of averaged values
%		sred		: standard deviation of averaged values
%
% uses :	intvers.m
%
% version 1.1.0		last change 23.07.1997

% G.Krahmann, IfM Kiel, Sep 1995
% added 'sred'		G. Krahmann, IfM Kiel, Oct 1996 	1.0.0 --> 1.0.1
% removed bug		G. Krahmann, IfM Kiel, May 1997		1.0.1 --> 1.0.2
% comp. to MATLAB 5	G. Krahmann, LODYC Paris, Jul 1997	1.0.2 --> 1.1.0


  if nargin<3
    dim=min(find(size(data)>1));
  end
  data=shiftdim(data,dim-1);
  s=size(data);
  kok=floor((s(1)-0.001)/n);		% avoid even divisions !!
  s2=s;
  s2(1)=kok+1;
  mred=nan*ones(s2);
  mind=nan*ones(s2);
  sred=nan*ones(s2);
  for k=1:kok
    [d1,d2]=nmean(data((k-1)*n+[1:n],:),1);
    mred(k,:)=d1;
    mind(k,:)=d2;
    sred(k,:)=nstd(data((k-1)*n+[1:n],:));
  end
  if kok*n+1~=s(1)
    [d1,d2]=nmean(data(kok*n+1:s(1),:),1);
    sred(kok+1,:)=nstd(data(kok*n+1:s(1),:));
  else
    d1=data(s(1),:);
    d2=0*d1+s(1);
    sred(kok+1,:)=0*d1;
  end
  mred(kok+1,:)=d1;
  mind(kok+1,:)=d2;
  mred=reshape(mred,s2);
  mind=reshape(mind,s2);
  sred=reshape(sred,s2);
  mred=shiftdim(mred,ndims(data)-dim+1);
  mind=shiftdim(mind,ndims(data)-dim+1);
  sred=shiftdim(sred,ndims(data)-dim+1);
