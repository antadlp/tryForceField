%A template for a single water molecule is created. The idea is
%to have the flexibility to position a variaty of water molecules
%emulating a thermal bath. In a real setting the water molecules
%have diferent orientations. The script attempt to create a
%single molecule with a specif orientation with a geometryc
%origin centered on the oxigen atom

%For the orientation it will be used the Euler Angles.

%The space units will be Ångström

% Input:
% 
% Or: Origin (xo, yo, zo)
% Or(1): must be xo
% Or(2): must be yo
% Or(3): must be zo
%
% {angle1, angle2, angle3}: Euler Angles
%
%


function [x y z] = singleWaterTemplate(Or, angle1, angle2, angle3)

d = 0.9584; % distance between Oxigen and a Hidrogen in
% equilibrium

O = [Or(1), Or(2), Or(3)];
H1 = [Or(1), Or(2)+d, Or(3)];


