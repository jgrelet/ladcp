function    [f,v,vel,cm,ea,pg,bt] = whread_multi(fnames);
% function    [f,v,vel,cm,ea,pg,bt] = whread_multi(fnames);
%
% load multiple workhorse files and join the data so that
% it appears as one data set. This is for handling cases
% when the ADCP restarted during a cast because of power
% failures.
%
% input  :  fnames          - array of filenames
%
% output :  
%
% version 0.1   last change 7.12.2006

% G.Krahmann, IFM-GEOMAR, Dec 2006


%
% display info
%
disp(' ')
disp('WHREAD_MULTI:  load and concatenate multiple files')


%
% loop over input files
%
for n=1:size(fnames,1)

    %
    % open file
    %
    [fid,message] = fopen(fnames(n,:),'r','l');

    
    %
    % load data
    %
    [fd,vd,veld,cmd,ead,pgd,btd] = whread(fid);

    
    %
    % append structures and arrays
    %
    if n==1
        f = fd;
        v = vd;
        vel = veld;
        cm = cmd;
        ea = ead;
        pg = pgd;
        bt = btd;
    else
        v = cat(1,v,vd);
        vel = cat(1,vel,veld);
        cm = cat(1,cm,cmd);
        ea = cat(1,ea,ead);
        pg = cat(1,pg,pgd);
        bt = cat(1,bt,btd);
    end        
    
    
    %
    % close file
    %
    fclose(fid);
end 