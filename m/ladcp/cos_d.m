function [result]=cosd(argu)
% function [result]=cosd(argu)
%
% cosine with degrees as argument

% Gerd Krahmann, IfM Kiel, Sep 1993

result=cos( argu / 180 * pi );
