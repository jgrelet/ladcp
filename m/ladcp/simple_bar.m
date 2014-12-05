function []=simple_bar(x,y,types)
% function []=simple_bar(x,y,types)
%
% a simple replacement for the regular stacked bar plot
% the regular Mathworks bar plot causes some difficulties with 
% legends and the storage as fig-file
%
% input  :	x		- row vector of x coordinates of bars
%		y		- array of lengths of bar elements
%		types		- character array of types of bar elements
%				  must contain one row per row of y
%
% version 0.2	last change 24.07.2008

% G.Krahmann, IFM-GEOMAR, June 2007

% colors repeat for more possibilities    GK, 24.07.2008   0.1-->0.2

%
% parse and check arguments
%
if length(x)~=size(y,2)
  error('x and y sizes are not as expected')
end


%
% 'integrate' y lengths
%
y = cumsum(y,1);


%
% set colors
%
cols = 'rgbkcm';
cols = [cols,cols];


%
% figure out linewidth
%
lw = ceil( 500/size(y,2) );


%
% plot bars
%
for n=1:size(y,1)
  if n==1
    dummy = plot([x;x],[0*y(1,:);y(1,:)],cols(n),'linewidth',lw);
    h(n) = dummy(1);
  else
    dummy = plot([x;x],[y(n-1,:);y(n,:)],cols(n),'linewidth',lw);
    h(n) = dummy(1);
  end
  hold on
end


%
% add legend, the MathWorks one causes trouble in 2006b
%
%legend(h,types)
ax = axis;
xx = linspace(ax(1),ax(1)+ax(2),length(h)+1);
for n=1:length(h)
  plot(xx(n),ax(4)+ax(4)/10,[cols(n),'.'],'markersize',20)
  text(xx(n),ax(4)+ax(4)/10,[' ',types(n,:)])
end
axis tight
axis(ax+[0,0,0,ax(4)/5])
