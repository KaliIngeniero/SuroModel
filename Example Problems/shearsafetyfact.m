function c = shearsafetyfact(design)
% c = shearsafetyfact(design)
%
% Calculates shear safety factor
%
% Copyright 2007 A Sobester
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

% ********* DESIGN VARIABLES ********
% Wire diameter
d = 0.5 + 6.5*design(1);

% Spring index
i = 4 + 12*design(2);

% Inter-coil distance coefficient
% kdelta = 0.1 + design(3);


% ******** PROBLEM DEFINITION *******

% Maximum and minimum load
Fmin = 40; Fmax = 500;

% Wire material coefficients
C3 = 0.625e5;
A1 = -0.19;

% Stress factor (Wahl factor)
kw = (4*i-1)/(4*i+4) + 0.615/i;

% Maximum/minimum shear stress [MPa]
tautmax = 8*kw*i*Fmax/(pi*d^2);
tautmin = 8*kw*i*Fmin/(pi*d^2);

% Amplitude/mean of the loading cycle [MPa]
taua = (tautmax - tautmin)/2;
taum = (tautmax + tautmin)/2;

% Yield shear strength [MPa]
Sys = C3*(1/25.4^2)*4.448*(d/25.4)^A1;

% Yield safety factor
SFy = Sys/(taua+taum);

% ********* CONSTRAINT **************
c = 1.05 - SFy;