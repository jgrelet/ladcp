function [data,messages] = misc_prepare_rdi(data,params,messages)
% function [data,messages] = misc_prepare_rdi(data,params,messages)
%
% do some simple data processing steps on the RDI raw data
% such as apply thresholds and cut the profile from raw data
%
% input  :	data		- LADCP data structure
%		params		- LADCP parameter structure
%		messages 	- LADCP message structure
%
% output :	data		- LADCP data structure
%		messages	- LADCP message structure
%
% version 0.1	last change 11.08.2005

% G.Krahmann, IFM-GEOMAR  Aug 2005

%
% apply any time offset to the LADCP
%
if params.timoff~=0
  data.tim = data.tim + params.timoff;
  disp(['  WARNING: adjusted ADCP time by ',num2str(params.timoff),' days']),
end


%
% we now have perfect Julian time
%
data.time_jul = data.tim(1,:);


%
% cut out the profile part from the raw data
% will add five minutes before and after time_start time_end
% so that the exact cutting is left for later
% Cutting here is only for the separation of multiple profiles
% within one raw data file
%
if ~isempty(params.time_start) & ~isempty(params.time_end)
  time_start = julian(params.time_start) - 5/1440;
  time_end = julian(params.time_end) + 5/1440;
  good = ( data.time_jul>=time_start & data.time_jul<=time_end );
  data = cutstruct(data,good);
end


%
% remove huge values in internal ADCP pressure data
%
if ~isempty(data.pres)
  data.pres = replace( data.pres, data.pres>10000, nan );
end


%----------------------------------------------------------------------
% Asynchronous ping interference detection
%----------------------------------------------------------------------
if params.detect_asynchronous==1
  disp('    Asynchronous ping interference detection')
  % determine average shape of combined TS profile
  dts_start = data.ts;
  for n=1:5
    dts = dts_start;
    dts = dts - ones(size(dts,1),1)*nmin(dts);
    [dummy,ind] = nmax(dts);
    dummy = hist(ind,[1:size(dts,1)]);
    [dummy,ind] = nmax(dummy);
    dtspeak = dts(ind,:);
    bad = find(isnan(dtspeak));
    good = find(~isnan(dtspeak));
    if ~isempty(bad)
      dtspeak(bad) = interp1(good,dtspeak(good),bad,'nearest');
    end
    dts = dts./(ones(size(dts,1),1)*dtspeak);
    dts = dts+1;
    dts = dts./(nmean(dts')'*ones(1,size(dts,2)));
    bad = find(dts>1.3);
    dts(bad) = nan;
    dts_start(bad) = nan;
  end
  bad = find(isnan(dts));
  data.ts(bad) = nan;
  data.weight(bad) = nan;
  data.ru(bad) = nan;
  data.rv(bad) = nan;
  data.rw(bad) = nan;
  data.re(bad) = nan;
  data.u(bad) = nan;
  data.v(bad) = nan;
  data.w(bad) = nan;
  data.e(bad) = nan;
end
