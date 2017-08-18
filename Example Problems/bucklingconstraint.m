function F=bucklingconstraint(x)
% F=bucklingconstraint(x)
%
% Calculates buckling constraint
%
% Global variables used:
%   BeamProperties
%
% Copyright 2007 A I J Forrester
%
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.

global BeamProperties
b=x(1)*0.045+0.005;
h=x(2)*0.23+0.02;
Iz=b^3*h/12;
It=(b^3*h+b*h^3)/12;
F=(BeamProperties.F*BeamProperties.SF)-(4/BeamProperties.L^2)*sqrt((BeamProperties.G*It*BeamProperties.E*Iz)/(1-BeamProperties.Nu^2));