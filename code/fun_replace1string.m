function fun_replace1string(path,name,str1,str2)

% ######  comma2point  ######
% changes one string (e.g. ',') to another (e.g. '.')
%
% Inputs:
%    pfad: directory at which the data is located
%    name: name of data which needs to be loaded
%    str1: string to be replaced
%    str2: string which replaces the other
%
% Outputs: changed file
%
% Other m-files required: -
% functions: -
% MAT-files required: -
%
% Author: Sonja Tesschemacher
% email: sonja.teschemacher@tum.de
% August 2016; Last revision: 26-Oct-2019


cd(path)

file = memmapfile(name,'Writable',true);
comma = uint8(str1);
point = uint8(str2);
file.Data(( file.Data==comma)' ) = point;

clearvars