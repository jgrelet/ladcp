function a=whoami
% function a=whoami
if isunix
  [s,a] = system('whoami');
  b = double(a);
  a(b < 32) = []; % Get Rid of Control Characters
else
  a = 'unknown';
end
