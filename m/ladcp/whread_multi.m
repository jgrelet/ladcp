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
% version 0.3   last change 24.05.2011

% G.Krahmann, IFM-GEOMAR, Dec 2006

% pad gaps in the data with NaN           GK, 21.11.2009  0.1-->0.2
% bug in multiple gap padding             GK, 24.05.2011  0.2-->0.3


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
        mstep = nmedian(diff(v(:,1,1)));
        tmiss = vd(1,1,1)-v(end,1,1);
        nmiss = floor(tmiss/mstep)-1;
        if nmiss>1
          disp(['>   Found gap between files, padding the gap with ',...
            int2str(nmiss),' NaNs'])
          clear vn
          for m=1:12
            dummy = linspace(v(end,1,m),vd(1,1,m),nmiss+2);
            vn(:,1,m) = dummy(2:end-1);
          end
          veln = repmat(vel(1,:,:)*nan,[nmiss,1,1]);
          cmn = veln;
          ean = veln;
          pgn = veln;
          btn = repmat(bt(1,:,:)*nan,[nmiss,1,1]);
          v = cat(1,v,vn);
          vel = cat(1,vel,veln);
          cm = cat(1,cm,cmn);
          ea = cat(1,ea,ean);
          pg = cat(1,pg,pgn);
          bt = cat(1,bt,btn);
        end
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