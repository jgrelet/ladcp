THIS IS A BETA VERSION

It works for me and on a number of test settings.
Please try it and report any problems or bugs to

gkrahmann@ifm-geomar.de

Gerd Krahmann, Kiel July 2006


----------------------------------------------------------

First a DISCLAIMER:

This software is NOT A COMMERCIAL package. It is provided to
you at no cost, but also without any guarantees for correct
results. We do, however, use it ourselves ;-)

Please also note that this is work in progress. There are 
likely problems and bugs in the software and I will be changing
and improving the code in the foreseeable future.

Any suggestions, comments, and even problems you should 
have are welcome. To simplify the all-important feedback, I
have created a Yahoo ladcp group. You are most welcome to 
become an active member of the group.

http://groups.yahoo.com/search?query=ladcp&submit=Search

If you do not want to become a member, you can also send an
e-mail directly to me. Please note that I have other work
to do, and that an answer might take a while.

gkrahmann@ifm-geomar.de


Gerd Krahmann

-------------------------------------------------------------------


		LADCP Processing v10.0-v10.2 beta

This is a description of how to set up the IFM-GEOMAR/LDEO
Matlab LADCP-Processing system. During the last rewrite of 
the code we decided to separate more strictly the 
cruise/ship/institute dependent part of the processing from
the actual inversion routines. For this we create a rather
strict interface through which data is delivered to the
processing. In order for this to work we need you, the user,
to create conversion routines to get your data into the
proper format. Once that is done, the processing should in
many cases be able to figure out many parameters by itself
and should proceed to the final result without any interference
or further input by the user.


1. SETUP OF THE LADCP SYSTEM

1.1 PROPER TIME BASE

Most important for the proper processing is to establish a proper
timing of the cast. In this version we rely on a perfect time
base of the LADCP (the time base should always be UTC). You
will be able to correct for wrong time bases, should
you know that your LADCP's clocks were set wrong.
(Set the parameter   params.timoff   to the offset in decimal days,
time base is always the clock of the downlooking instrument in
cases of a dual head system. This should be done in
cruise_id/cast_params.m )

While navigational data, typically GPS data,
is in all but the rarest cases, on a proper time base, this
is not necessarily true for other data. So when you intend
to use this software and set up your LADCP system, please make
100% sure that your LADCP clocks are set correctly (at least to
a few seconds). Since the ADCP's clocks can drift noticably, you
will need to check and possibly set the clocks for every single
cast. A time server on the ship to which you synchronize your
LADCP control computer on a regular basis (before each cast
or once a day) is the safest way to ensure the proper time base.

CTD-time data which we, when available, use to 
infer the depth of the system should similarly be on the same
time base, though our routines are able to correct for a different
time base in the CTD data.


2. SETUP OF THE PROCESSING

2.0 WHICH DATA DO WE NEED FOR PROCESSING

Lowered ADCP processing can be done with just the ADCP data
and logged times and positions when the casts begin and end.
This very basic information will result in water velocity 
profiles. Shear profiles can even be determined with just the
ADCP data.
This is, however, not the best result one can obtain. 

If we add the processed CTD profile, that in most cases was collected 
together with the LADCP data, we can correct for sound speed variations
within the water column. This will lower the error of the final
profile.

If one has access to it, the raw CTD data (not yet deptp binned, and
not loop-edited with Seabird software) is another highly recommneded
data addition. It will be used for locating the LADCP in the vertical
axis. If this data is not available, the depth will be inferred from
the vertical velocity measured by the ADCP. Usually not the best way.

If a ship-ADCP was installed on the ship and running while collecting
the LADCP profile, this data can be added to further reduce the errors.

Navigational data, that is on most research vessels regularly collected
can also be fed into the processing and will save us the need to
enter start and end positions by hand. On some ships this data is
being fed into Seabird CTD systems. That is the best way for
LADCP processing.


2.1 CREATING THE DIRECTORY STRUCTURE

To maintain the separation of ship/cruise/institute dependent
and independent parts of the processing we introduce from this
version on a new directory structure. This will allow you to
easily install a new processing software by replacing just one
directory.
We suggest you start the installation in
a base directory named    'ladcp'   or similarly.

The directory structure thereafter will be as follows. A script 
(create_cruise.m) is provided to create this directory structure 
and copy some relevant files at the proper places.

ladcp				- base directory
    m				- ship/cruise/institute independent software	
	ladcp			- ladcp processing
	sw			- CSIRO's seawater routines
	initial_dir		- extra stuff to create new cruise directories
	netcdf			- this is the netcdf toolbox
	mexnc			- these are the platform independent
				  mex files for the netcdf toolbox
				  included here are only Linux mex files
				  all other operating system mex files
				  need to be installed by the user from
				  http://mexcdf.sourceforge.net/
    cruise_id			- ship/cruise/institude dependent files and data
				  also contains the startup.m script and
				  cruise_params.m and cast_params.m
	m			- user-modified loading routines
	logs			- processing logs
	profiles		- resulting profiles
	plots			- resulting plots
	tmp			- temporary files
	data			- input data
	    raw_ladcp		- LADCP raw data as from RDI
	    raw_sadcp		- SADCP raw data as from its processing
	    raw_nav		- navigational data as from its processing
	    raw_ctdprof		- CTD profiles as from its processing
	    raw_ctdtime		- CTD time data as from its processing
	    ladcp		- fixed format LADCP mat-files for fast loading
				  (not yet implemented)
	    sadcp		- fixed format SADCP mat-files for fast loading
	    nav			- fixed format NAV mat-files for fast loading
	    ctdprof		- fixed format CTD profile mat-files 
	    ctdtime		- fixed format CTD time mat-files

The base directory ( ladcp )will contain the m-file 'create_cruise.m' .
Execute this file to create a new directory structure under a different
'cruise_id' name. You will then have exactly one subdirectory 
named 'cruise_id' for each cruise.
'create_cruise.m' will also copy simple templates for the conversion
routines from your specific data formats to the fixed formats needed by
the LADCP processing. 

You MUST quit matlab before starting the processing.
This software relies on matlab being started in the cruise_id directory.


2.2 SET UP THE CRUISE SPECIFIC PART OF THE PROCESSING

You will need to edit EACH of the following files
cruise_id/m/prepladcp.m
cruise_id/m/prepsadcp.m
cruise_id/m/prepnav.m
cruise_id/m/prepctdtime.m
cruise_id/m/prepctdprof.m

These files contain instructions what you have to change within.
We have collected some examples in 'm/initial_dir/examples/'

prepladcp.m	will copy the raw ADCP data from a location determined
		by you to the location where the processing expects it.
		Here you can also rename the files to the convention
		used within the processing.

prepsadcp.m	will copy raw Ship-ADCP data, if it is available,
		to the location where the processing expects it.
		The routine will extract some data and store it
		in a mat-file with set variable names.

prepnav.m	will copy raw Navigation data, if it is available,
		to the location where the processing expects it.
		The routine will extract some data and store it
		in a mat-file with set variable names.

prepctdprof.m	will copy processed, depth binned, CTD data, if it is
		available, to the location where the processing expects it.
		The routine will extract some data and store it
		in a mat-file with set variable names.

prepctdtime.m	will copy not depth-binned CTD data, if it is
		available, to the location where the processing expects it.
		The routine will extract some data and store it
		in a mat-file with set variable names.

If you have done this once and you do another cruise with the same
setup, you can of course reuse your old prep*.m files with only
minor modifications.

In 'cruise_id/cruise_params.m' you should then change processing
parameters that are the same throughout your cruise.

In 'cruise_id/cast_params.m' you should change processing
parameters that are specific for a single cast. E.g. start and
end position, if you enter them manually.
This is the only file where you might have to enter information
on a per cast basis (usually not necessary if you provide correct 
time bases and navigational data).


2.3 PROCESSING A PROFILE

To process a cast you need to do the following steps
- copy the raw data files to the proper location
- start matlab in the 'cruise_id' directory
- if necessary, enter start and end position and time into
  'cruise_id/cast_params.m'
- call 'process_cast(stnno)'


2.4 REPROCESS A PROFILE

A profile can simply be reprocessed by calling 'process_cast(stnno)'
again. To speed up the processing it will load the previously 
prepared mat-files with the data. If you want to rerun the
preparation of data too, first execute 'clear_prep(stnno)' .
This second part NEEDS to be done, in case you are fixing problems
with the raw data. The reason is that some of the data preparation
depends on the previously prepared data.
