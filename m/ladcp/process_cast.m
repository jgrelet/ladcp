function [] = process_cast(stn,extraarg)
% function [] = process_cast(stn,extraarg)
%
% Process LADCP cast, including GPS, SADCP, and BT data.
%
% input  :    stn             - station number, all data needs to be
%                               numbered in this way
%                               Enter a negative number to clear loaded
%                               and saved ancillary data for that station.
%             extraarg  []    - add 'noplots' to not save plots
%                               useful to speed up testing
%
% version 0.8	last change 13.01.2013
 
% G.Krahmann, IFM-GEOMAR

% orient statements for figure saving     GK, 14.07.2008  0.1-->0.2
% orient gone, several drawnow's          GK, 15.07.2008  0.2-->0.3
% call clear_prep for negative stn        GK, 07.02.2009  0.3-->0.4
% variable p.print_formats                GK, 28.05.2011  0.4-->0.5
% noplots                                 GK, 04.11.2012  0.5-->0.6
% f as output of rdiload.m                GK, 08.11.2012  0.6-->0.7
% act on p.savemat                        GK, 13.01.2013  0.7-->0.8

%
% handle extra arguments
%
noplots = 0;
if nargin>1
  if strcmp(extraarg,'noplots')
    noplots = 1;
  end
end


%
% check current directory
%
if exist('logs')~=exist(pwd);
    disp('>   This directory is not prepared for the LADCP software ')
    disp('>   Sorry EXIT ')
    return
end


%
% clear already loaded and saved data for reloading and reprocessing
%
if stn<0
%  clear_prep(-stn);
%  return
   disp('please use clear_prep');
   return
end


%
% Initialize the processing by loading parameters
% and make sure that we have no leftovers from previous processings
%
default_params;
cruise_params;
cast_params;
files = misc_composefilenames(p,stn);

disp(p.software);


% 
% prepare the various data files for easy loading
%
[values] = prepare_cast(p,files);


%
% load RDI data
% and perform some simple processing steps
%
[data,values,messages,p,files] = rdiload(files,p,messages,values);
[data,messages] = misc_prepare_rdi(data,p,messages);

tic;					% start timer

% 
% plot the display menu
%
plot_menu
drawnow


%
% convolution of the loading routines for the 
% various data sets
%
[data,messages] = loading(files,data,messages,p);


%
% merge LADCP data with NAV and CTD
%
[data,p,messages,values] = mergedata(data,p,messages,values);


% 
% improve the data quality by removing spikes etc
%
[data,p,values,messages] = improve(data,p,values,messages);


%
% extract the bottom track
%
[data,p] = getbtrack(data,p,values);  


%
% Find the depth and bottom and surface using ADCP data 
%
% This needs to be done twice. First run without sound speed
% correction. Then apply the depth dependent sound speed 
% correction. And then recalculate the depths.
%
[data,p,values,messages] = calc_depth(data,p,values,messages);
[data,p,values,messages] = calc_soundsp(data,p,values,messages);
[data,p,values,messages] = calc_depth(data,p,values,messages);
drawnow


%
% cut off the parts before and after the main profile
% i.e. the surface soak and similar parts
%
[data,p,values,messages] = misc_cut_profile(data,p,values,messages);
drawnow

%
% for the real profile time period extract the SADCP data
% and average it
%
[data,messages] = calc_sadcp_av(data,p,values,messages);
drawnow


%
% Plot a summary plot of the raw data
%
plot_rawinfo( data, p, values );
drawnow


%
% apply some editing of the single bins
%
data = edit_data(data,p,values);
drawnow


%
% form super ensembles
%
%data1 = data;
[p,data,messages] = prepinv(messages,data,p,[],values);
%[p,data1,messages] = prepinv_with_old_rotation_options(messages,data1,p,[],values);
[di,p,data] = calc_ens_av(data,p,values);
%[di1,p,data1] = calc_ens_av(data1,p,values);
drawnow


%
% remove super ensemble outliers
%
if ps.outlier>0 | p.offsetup2down>0
  [messages,p,dr,de,der] = lanarrow(messages,values,di,p,ps);
%  [messages,p,dr1,de1,der1] = lanarrow(messages,values,di1,p,ps);
end


%
% once we have a first guess profile we recompute the super ensemble
%
if (p.offsetup2down>0 & length(data.izu)>0)
  [p,data,messages] = prepinv(messages,data,p,dr,values);
%  [p,data1,messages] = prepinv_with_old_rotation_options(messages,data1,p,dr1,values);
%  keyboard
  [di,p,data] = calc_ens_av(data,p,values);
end

 
% 
%  take advantage of presolve if it existed  ?? GK
%  call the main inversion routine
%
[messages,p,dr,ps,de] = getinv(messages,values,di,p,ps,dr,1);
drawnow


%
% check inversion constraints
% 
p = checkinv(dr,de,der,p,ps,values);
if isfield(de,'bvel') 
  p = checkbtrk(data,di,de,dr,p); 
end


%
% Compute 'old fashioned' shear based solution 
%  two choices, fisrt us all data
%  second use super ensemble data
%
if ps.shear>0
  if ps.shear==1
    [ds,dr,ps,p,messages] = calc_shear2(data,p,ps,dr,messages);
  else
    [ds,dr,ps,p,messages] = calc_shear2(di,p,ps,dr,messages);
  end
end


%
% Plot final results
%
plot_result(dr,data,p,ps,values)
drawnow


%  
% Convert p.warn to one line of text with newline characters
%
p.warnings = [];
for n = 1:size(messages.warnp,1)
  p.warnings = [p.warnings deblank(messages.warnp(n,:)) char(10)];
end
 
sfigure(2);
clf
% experimental diagnostic of battery voltage
%
%[p,messages] = calc_battery(p,values,messages);
  
%
% complete task by repeating the most important warnings
%
if size(messages.warn,1) + size(messages.warnp,1) > 2
  disp(' ')
  disp(messages.warn)
  disp(' ')
  disp(messages.warnp)
  for j=1:size(messages.warn,1)
    text(0,1.1-j/10,messages.warn(j,:),'color','r','fontsize',14,'fontweight','bold')
  end
  for j=1:size(messages.warnp,1)
    text(0,1.1-(size(messages.warn,1)+1+j)/10,messages.warnp(j,:),'color','r','fontsize',14,'fontweight','bold')
  end
  axis off
else
  text(0,1.1-1/10,'LADCP profile OK','color','g','fontsize',30,'fontweight','bold')
  axis off
end
  
streamer([p.name,' Figure 11']);
hgsave('tmp/11')


%----------------------------------------------------------------------
% STEP 18: SAVE OUTPUT
%----------------------------------------------------------------------

disp(' ')
disp('SAVING RESULTS')
if length(files.res)>1
  
  %
  % save results to ASCII, MATLAB and NETCDF files
  %
  saveres(data,dr,p,ps,files,values)
%  da = savearch(values,dr,data,p,ps,f);

  %
  % save plots
  %
  % handle newer matlab versions
  if version('-release')>=14
     imac = 0;
  else
    imac = ismac;
  end

  if noplots==0
    for n = 1:length(p.saveplot)
      j = p.saveplot(n);
      if exist(['tmp/',int2str(j),'.fig'],'file')
        figload(['tmp/',int2str(j),'.fig'],2)
      end
      warning off
      if imac
        if findstr(p.print_formats,'ps')
          eval(['print -depsc ',files.plots,'_' int2str(j) '.eps '])
        end
      else
        if findstr(p.print_formats,'ps')
          eval(['print -dpsc ',files.plots,'_' int2str(j) '.ps '])
        end
      end
      if findstr(p.print_formats,'jpg')
        eval(['print -djpeg ',files.plots,'_' int2str(j) '.jpg '])
      end
    if findstr(p.print_formats,'png')
       eval(['print -dpng ',files.plots,'_' int2str(j) '.png '])
      end
      warning on
    end
  end
  
  % save a protocol
  saveprot

  % save full information into mat file
  if p.savemat==1
    disp(['    Saving full information to ',files.res,'_full.mat'])
    save6([files.res,'_full.mat'])
  end
  
end

% switch to final result figure
plot_controls(1)
    

%----------------------------------------------------------------------
% FINAL STEP: CLEAN UP
%----------------------------------------------------------------------

fclose('all');				%  close all files just to make sure

disp(' ')				% final message
disp(['    Processing took ',int2str(toc),' seconds'])

save6 tmp/last_processed

diary off
