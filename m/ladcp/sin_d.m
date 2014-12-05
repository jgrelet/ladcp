function [result]=sind(argu)
% function [result]=sind(argu)
%
% sine with degrees as argument

% Gerd Krahmann, IfM Kiel, Sep 1993

result=sin( argu / 180 * pi );
