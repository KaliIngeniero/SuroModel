function W = liftsurfw4(design)
% W = liftsurfw(design)
%
% Four dimensional version of liftsurfw.
% The four most significant variables (as chosen via the Morris screening
% algorithm) are active, the rest are kept at the baseline value.
% 
% S_w    = design(1);
% W_fw   = design(2);
% A      = design(3);
% Lambda = design(4)*pi/180;
% q      = design(5);  
% lambda = design(6);
% tc     = design(7);   
% N_z    = design(8); 
% W_dg   = design(9);
% W_p    = design(10); 
%
% Symbol  Parameter      Typical light aircraft value (C172)
% S_w    -Wing area (ft^2)                              174
% W_fw   -Weight of fuel in the wing (lb)               252
% A      -aspect ratio                                 7.52
% Lambda -quarter-chord sweep (deg)                       0
% q      -dynamic pressure at cruise (lb/ft^2)           34
% lambda -taper ratio                                 0.672
% tc     -aerofoil thickness to chord ratio            0.12
% N_z    -ultimate load factor (1.5x limit load factor) 3.8
% W_dg   -flight design gross weight (lb)             ~2000
% W_p    -paint weight (lb/ft^2)                     ~0.064
% *********************************************************
% W      -wing weight lbs (est. by the eq./actual)  245/236
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

S_w    = design(1);
W_fw   = 252;
A      = 7.52;
Lambda = 0*pi/180;
q      = 34;  
lambda = 0.672;
tc     = design(2);   
N_z    = design(3); 
W_dg   = design(4);
W_p    = 0.064; 

W = 0.036*S_w^0.758*W_fw^0.0035*(A/cos(Lambda)^2)^0.6*q^0.006*lambda^0.04*...
    (100*tc/cos(Lambda))^-0.3*(N_z*W_dg)^0.49 + S_w*W_p;
