function [] = plot_controls(fig)
% function [] = plot_controls(fig)
%
% control function for LADCP plotmenu
%
% reloads the stored figure into the display window

global mh

figload(['tmp/',int2str(fig),'.fig'],2)

figure(1)

for n=1:length(mh)
  if mh(n)~=0
    set(mh(n),'foregroundcolor',[0,0,0]);
  end
end
set(mh(fig),'foregroundcolor',[1,0,0]);
