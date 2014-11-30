function handle = streamer(fig,TitleString)
%  STREAMER  Titles for an entire figure.
% 	     STREAMER('text') adds text at the top of the current figure,
% 	     going across subplots.
%        STREAMER(fig,'text') adds it to the specified figure.
%
% changed default to fontsize 14 and fontweight bold   GK
% 
% 	     See also XLABEL, YLABEL, ZLABEL, TEXT, TITLE.
%
%	Keith Rogers 11/30/93
%
% modified for LADCP processing, GK

% Copyright (c) by Keith Rogers 1995

% 
% Mods:
%	11/94 adapted to 4.2
%   06/95 clean up, added alternate figure option.

if(nargin<2)
	TitleString = fig;
	fig = gcf;
end
ax = gca;
sibs = get(fig, 'Children');
ind = findobj(sibs,'Type','axes','Tag','Streamer');
if isempty(ind)
   StreamerHand = axes('Parent',fig,...
                       'Units','normalized',...
                       'Position',[.1 .9 .8 .05],...
                       'Box','off',...
                       'Visible','off',...
                       'Tag','Streamer');
else
   StreamerHand = sibs(ind);
end
title(TitleString,'interpreter','none','Visible','on',...
                  'fontweight','bold','fontsize',16,...
                  'parent',StreamerHand);
axes(ax);
