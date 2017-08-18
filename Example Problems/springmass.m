function mass = springmass(design)
% mass = springmass(design)
%
% Calculates the weight of a helical compression spring
% defined by the design variable vector design
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
kdelta = 0.1 + design(3);



% ******** PROBLEM DEFINITION *******

% Wire steel density [kg/mm^3]
rho = 7.87e-6;


% Rigidity modulus
G = 0.78e5;

% Maximum and minimum load
Fmin = 40; Fmax = 500;

% Stroke (working range) [mm]
h = 50;


% Stiffness
c = (Fmax-Fmin)/h;

% Maximum / minimum deflection [mm]
fmax = Fmax/c; 

% Spring mean diameter
Dm = i*d; 

% Number of active coils
n = (G*d^4)/(8*c*Dm^3);

% Number of end coils
if n<=7
    nr = 1.5;
else
    nr = 2;
end

% Total number of spring coils
nt = n + nr;

% Pitch of unloaded spring [mm]
t = d + fmax/n + kdelta*d;

% Helix angle of the unloaded spring [rad]
alpha0 = atan(t/(pi*Dm));

% Wire length [mm]
ls = pi*Dm*nt/cos(alpha0);

% ********* OBJECTIVE ****************
% Spring total mass
mass = ls*0.25*pi*d^2*rho;