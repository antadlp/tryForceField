%A template for a single water molecule is created. The idea is
%to have the flexibility to position a variaty of water molecules
%emulating a thermal bath. In a real setting the water molecules
%have diferent orientations. The script attempt to create a
%single molecule with a specif orientation with a geometryc
%origin centered on the oxigen atom

%For the orientation it will be used rotation matrix.

%The space units will be Ångström

% Input:
% 
% Or: Origin (xo, yo, zo)
% Or(1): must be xo
% Or(2): must be yo
% Or(3): must be zo
%
% {anglex, angley, anglez}: around an axis
% the angle units are degrees
%
%


function [O H1 H2] = singleWaterTemplate(Or, anglex, angley, anglez)

d = 0.9584; % distance between Oxigen and a Hidrogen in
% equilibrium

O = [Or(1), Or(2), Or(3)]';
H1 = [Or(1), Or(2)+d, Or(3)]';

H1p = H1-O;

theta = 104.45;

seno = sind(theta);
coseno = cosd(theta);

H2p = rotz(-104.45)*H1p;
H2 = O + H2p;

% rotations

% x-axis

H1 = O + rotx(anglex)*H1p;
H2 = O + rotx(anglex)*H2p;
%


H1p = H1-O;
H2p = H2-O;

% y-axis

H1 = O + roty(angley)*H1p;
H2 = O + roty(angley)*H2p;

H1p = H1-O;
H2p = H2-O;

% z-axis

H1 = O + rotz(anglez)*H1p;
H2 = O + rotz(anglez)*H2p;






