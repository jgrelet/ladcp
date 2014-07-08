function [messages,p,dr,de,der] =...
	lanarrow(messages,values,di,p,ps)
%function [messages,p,dr,de,der] =...
%	lanarrow(messages,values,di,p,ps)
%
% remove super ensemble outliers
%
% input  :	messages	- structure containing messages
%		values		- structure containing derived values
%		di		- structure ?
%		p		- structure containing general parameters
%		ps		- structure containing solution parameters
% 
% output :	messages	- modified input structure
%		p		- modified input structure
%		dr		- ?
%		de		- ?
%		der		- ?
%
% version 0.1	last change 18.06.2008

% M.Visbeck, converted to function GK

disp(' ')
disp('LANARROW:  remove outliers in superensembles')


ps1 = ps;
ps1.down_up = 0;
ps1.solve = 0;
for n=1:ps.outlier
  if n>1
    [messages,p,dr,ps1,de,der] = getinv(messages,values,di,p,ps1,dr);
  else
    [messages,p,dr,ps1,de,der] = getinv(messages,values,di,p,ps1);
  end
  dif = (di.ru-der.ru_oce-der.ru_ctd).^2+...
         (di.rv-der.rv_oce-der.rv_ctd).^2;
  [es,ii] = sort(dif(:));
  iok = find(isfinite(es));
  ln = length(iok)*0.01;
  if ln>0
    di.weight(ii(iok(end-[0:ln]))) = NaN;
    disp(['    Giving low weight to the 1% of data scans deviating most'])
  end
end
