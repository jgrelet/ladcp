function jd = date2jd(varargin)
% Modified Julian day number from Gregorian date.
%
% check value 2451545.0 = January 1, 2000, 00:00:00
%
% this one differs from the correct Julian day by 0.5 days
%
%   JD = DATE2JD(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the Julian
%   day number of the given date (Gregorian calendar) plus a fractional part
%   depending on the time of day.
%
%   Any missing MONTH or DAY will be replaced by ones.  Any missing HOUR,
%   MINUTE or SECOND will be replaced by zeros.
%
%   If no date is specified, the current date and time is used.
%
%   Start of the JD (Julian day) count is from 0 at 12 noon 1 January -4712
%   (4713 BC), Julian proleptic calendar.  Note that this day count conforms
%   with the astronomical convention starting the day at noon, in contrast
%   with the civil practice where the day starts with midnight.
%
%   Astronomers have used the Julian period to assign a unique number to
%   every day since 1 January 4713 BC.  This is the so-called Julian Day
%   (JD).  JD 0 designates the 24 hours from noon UTC on 1 January 4713 BC
%   (Julian proleptic calendar) to noon UTC on 2 January 4713 BC.

%   Sources:  - http://tycho.usno.navy.mil/mjd.html
%             - The Calendar FAQ (http://www.faqs.org)

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-05-24 13:30:06 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   narginchk(0, 6);
   if nargsin
      argv = {1 1 1 0 0 0};
      argv(1:nargsin) = varargin;
   else
      argv = num2cell(clock);
   end
   [year, month, day, hour, minute, second] = deal(argv{:});

if nargin==1
  second = year(:,6);
  minute = year(:,5);
  hour = year(:,4);
  day = year(:,3);
  month = year(:,2);
  year = year(:,1);
end

   % The following algorithm is a modified version of the one found in the
   % Calendar FAQ.

   a = floor((14 - month)/12);
   y = year + 4800 - a;
   m = month + 12*a - 3;

   % For a date in the Gregorian calendar:
   jd = day + floor((153*m + 2)/5) ...
        + y*365 + floor(y/4) - floor(y/100) + floor(y/400) - 32045 ...
        + ( second + 60*minute + 3600*hour )/86400;
