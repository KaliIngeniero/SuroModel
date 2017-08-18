function c = buckling(design)
% c = buckling(design)
% Maximum deflection must be less than buckling deflection
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

% Young's Modulus
E = 2.06e5;

% Rigidity modulus
G = 0.78e5;

% Maximum and minimum load
Fmin = 40; Fmax = 500;

% End support parameter
v = 0.5;

% Stroke (working range) [mm]
h = 50;

% Stiffness
c = (Fmax-Fmin)/h;

% Maximum deflection [mm]
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

% Spring solid length
Hb = nt*d;

% Pitch of unloaded spring [mm]
t = d + fmax/n + kdelta*d;

% Spring free length [mm]
H0 = Hb + n*(t-d);

% Slenderness ratio
lambda = H0/Dm;

% Critial slenderness ratio coefficients
cf1 = E/(2*(E-G));
cf2 = 2*pi^2*(E-G)/(E+2*G);

% Critical slenderness ratio
lambda_critic = sqrt(cf2)/v;

% Theoretical buckling deflection
if lambda > lambda_critic
    ff = H0*cf1*(1-sqrt(1-cf2/(v^2*lambda^2)));
else
    ff = H0;
end

c = fmax-ff;