function []=saveres(data,dr,p,ps,f,values)
% function []=saveres(data,dr,p,ps,f,values)
%
% store LADCP result in RODB format
% and store some other info as mat file
%
% version 0.7	 last change 16.11.2012

% changed values stored as lat/lon            G.K. May 2007   0.2-->0.3
% modified version id                         GK  Jun 2007    0.3-->0.4
% change save command                         GK  Aug 2007    0.4-->0.5
% added some header info                      GK  Jun 2008    0.5-->0.6
% write target strength info into mat file    GK, 16.11.2012  0.6-->0.7

%
% extract some target strength data and transmit voltage
%
ts.dn.ts = data.ts_edited(data.izd(3),:);
ts.dn.z = data.z-data.zd(3);
ts.dn.xmv = data.xmv(1,:);
if values.up==1
  ts.up.ts = data.ts_edited(data.izu(3),:);
  ts.up.z = data.z+data.zu(3);
  ts.up.xmv = data.xmv(2,:);
end


%
% store some results as a MAT file
%
save6([f.res,'.mat'],'ts','dr','p','ps','f')


%
% open file
%
fid = fopen([f.res,'.lad'],'wt');


%
% write file header
%
fprintf(fid,['Filename    = %s\n'],f.res);
fprintf(fid,['Date        = %s\n'],datestr(p.time_start,26));
fprintf(fid,['Start_Time  = %s\n'],datestr(p.time_start,13));
%[lats,lons] = pos2str([p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60]);
%[lats,lons] = pos2str([values.lat,values.lon]);
fprintf(fid,['Latitude    = %s\n'],num2str(values.lat));
fprintf(fid,['Longitude   = %s\n'],num2str(values.lon));
fprintf(fid,['Deviation   = %f\n'],values.magdev);
fprintf(fid,['Version     = %s\n'],p.software);
fprintf(fid,['Processed   = %s\n'],datestr(clock));
fprintf(fid,['Units       = m:m/s:m/s:m/s\n'],[]);
fprintf(fid,['Columns     = z:u:v:ev\n'],[]);
if ~isfield(dr,'uerr')
 dr.uerr=dr.u*NaN;
end


%
% write data
%
fprintf(fid,['%6.1f %6.3f %6.3f %6.3f \n'],[dr.z,dr.u,dr.v,dr.uerr]');


%
% close file
%
fclose(fid);


%
% in case we have bottom track data, store that in another file
%
if isfield(dr,'ubot')

  % save bottom track data
  % open file
  fid = fopen([f.res,'.bot'],'wt');

  fprintf(fid,['Filename    = %s\n'],f.res);
  fprintf(fid,['Date        = %s\n'],datestr(p.time_start,26));
  fprintf(fid,['Start_Time  = %s\n'],datestr(p.time_start,13));
  [lats,lons] = pos2str([p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60]);
  fprintf(fid,['Start_Lat   = %s\n'],lats);
  fprintf(fid,['Start_Lon   = %s\n'],lons);
  fprintf(fid,['Deviation   = %f\n'],values.magdev);
  fprintf(fid,['Bottom depth= %d\n'],fix(p.zbottom));
  fprintf(fid,['Version     = %s\n'],p.software);
  fprintf(fid,['Processed   = %s\n'],datestr(clock));
  fprintf(fid,['Units       = m:m/s:m/s:m/s\n'],[]);
  fprintf(fid,['Columns     = z:u:v:err\n'],[]);
  fprintf(fid,['%6.1f %6.3f %6.3f %6.3f\n'],...
            [dr.zbot,dr.ubot,dr.vbot,dr.uerrbot]');

  fclose(fid);

end
