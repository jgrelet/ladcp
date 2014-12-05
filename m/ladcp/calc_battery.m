function [p,messages] = battery(p,values,messages)
% function [p,messages]=battery(p,values,messages)
% try to find know battery calibration and issue warning if low


%
% general function info
%
disp(' ')
disp('CHECK_BATTERY:  display the battery voltage')


%
% define battery levels: first low   second critical
%
batlevel=[40 37];


%
% which instrument to compare with
%
if ~isnan(values.xmv(1))
  instrid = values.instid(1);
  n_instr = 1;
elseif length(values.xmv)==2
  instrid = values.instid(2);
  n_instr = 2;
else
  instrid = values.instid(1);
  n_instr = 1;
end
 

%
% load battery info file
%
if exist('m/batt_info.dat','file')
  load m/batt_info.dat
else
  batt_info = [];
end


%
% look whether current instrument is in there
%
data = [];
p.battery = nan;
for n=1:size(batt_info,1)
  if instrid==batt_info(n,1)
    disp(['    Found CPU board of serial ',int2str(batt_info(n,2))])
    data = batt_info(n,:);
    p.battery = data(3)*values.xmv(n_instr);
  end
end

if isempty(data) 
  if ~isnan(values.xmv(n_instr))
    disp('>   Do not have calibration coefficient for this instrument.')
    disp('>     To get battery information in the future, please measure the')
    disp('>     battery voltage after the cast and divide it by')
    disp(['>     ',num2str(values.xmv(n_instr))])
    disp('>     Then enter the line:')
    disp(['>     ',int2str(instrid),'   Your-Serialnumber   Your-Value'])
    disp('>     in the file  m/batt_info.dat  .')
    disp('>     You can also send the info to gkrahmann@ifm-geomar.de')
    disp('>     so that it will be included in future software updates.')
  else
    disp('    No battery info available since this is an old instrument.')
  end
end

if p.battery>batlevel(1)
  bc = 'g';
  disp(['    Battery Voltage is ',num2str(round(p.battery*10)/10),' V'])
elseif p.battery>batlevel(2)
  bc = 'y';
  disp(['    Battery Voltage is ',num2str(round(p.battery*10)/10),' V'])
else
   if ~isnan(p.battery)
     warn = ['>   Battery voltage is low : ',...
	num2str(round(p.battery*10)/10),' V'];
     messages.warn = strvcat(messages.warn,warn);
     bc = 'r';
  else
     warn = ['>   Battery voltage is unknown. Need calibration coeff.'];
     messages.warn = strvcat(messages.warn,warn);
     bc = 'y';
  end
  disp(warn)
end
text(0,0,['Battery Voltage is ',num2str(round(p.battery*10)/10),' V'],...
	'color',bc,'fontsize',14,'fontweight','bold')
