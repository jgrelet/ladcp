% initialize the plot menu for LADCP processing
%
% version 0.3	last change 15.07.2008

% G.Krahmann, IFM-GEOMAR, 2004

% changed location of windows for OSX     GK, 15.07.2008   0.2--0.3

dd = dir('tmp/*.fig');
for n=1:length(dd)
  delete(['tmp/',dd(n).name])
end
global mh

% create 2 windows
% one for the control menu and one for the actual plots
figure(1)
clf
set(gcf,'position',[10,10+ismac*100,170,750],'numbertitle','off','menubar',...
	'none','name','LADCP 1');
figure(2)
clf
set(gcf,'position',[190,10+ismac*100,800,696],'numbertitle','off',...
	'name','LADCP 2');

% create the menu
figure(1)
fh(1) = uicontrol('style','frame','position',[10,10,150,600]);
fh(2) = uicontrol('style','frame','position',[10,620,150,120]);

th(1) = uicontrol('style','text','position',[15,705,140,30],...
	'horizontalalignment','center','string','LADCP',...
	'fontsize',18);
th(2) = uicontrol('style','text','position',[15,665,140,30],...
	'horizontalalignment','center','string','processing',...
	'fontsize',18);	
th(3) = uicontrol('style','text','position',[15,625,140,30],...
	'horizontalalignment','center','string','display',...
	'fontsize',18);	
	
mh(1) = uicontrol('style','push','position',[15,585,140,20],...
	'horizontalalignment','center','string','UV-profiles',...
	'callback','plot_controls(1)');
mh(2) = uicontrol('style','push','position',[15,560,140,20],...
	'horizontalalignment','center','string','W,Z,sensors',...
	'callback','plot_controls(2)');
mh(3) = uicontrol('style','push','position',[15,535,140,20],...
	'horizontalalignment','center','string','ERR-profiles',...
	'callback','plot_controls(3)');
mh(4) = uicontrol('style','push','position',[15,510,140,20],...
	'horizontalalignment','center','string','SURF/BOT recog.',...
	'callback','plot_controls(4)');
mh(5) = uicontrol('style','push','position',[15,485,140,20],...
	'horizontalalignment','center','string','Heading dual',...
	'callback','plot_controls(5)');
mh(6) = uicontrol('style','push','position',[15,460,140,20],...
	'horizontalalignment','center','string','HDG/PIT/ROL',...
	'callback','plot_controls(6)');
mh(7) = uicontrol('style','push','position',[15,435,140,20],...
	'horizontalalignment','center','string','CTD motion',...
	'callback','plot_controls(7)');
mh(8) = uicontrol('style','push','position',[15,410,140,20],...
	'horizontalalignment','center','string','CTD-LADCP merging',...
	'callback','plot_controls(8)');
mh(9) = uicontrol('style','push','position',[15,385,140,20],...
	'horizontalalignment','center','string','SADCP',...
	'callback','plot_controls(9)');
mh(10) = uicontrol('style','push','position',[15,360,140,20],...
	'horizontalalignment','center','string','Vel offset',...
	'callback','plot_controls(10)');
mh(11) = uicontrol('style','push','position',[15,335,140,20],...
	'horizontalalignment','center','string','Warnings',...
	'callback','plot_controls(11)');
mh(12) = uicontrol('style','push','position',[15,310,140,20],...
	'horizontalalignment','center','string','Inv. Weights',...
	'callback','plot_controls(12)');
mh(13) = uicontrol('style','push','position',[15,285,140,20],...
	'horizontalalignment','center','string','Bottom Track',...
	'callback','plot_controls(13)');
mh(14) = uicontrol('style','push','position',[15,260,140,20],...
	'horizontalalignment','center','string','TargetStr Edited',...
	'callback','plot_controls(14)');
mh(15) = uicontrol('style','push','position',[15,235,140,20],...
	'horizontalalignment','center','string','Correl. Edited',...
	'callback','plot_controls(15)');
	
