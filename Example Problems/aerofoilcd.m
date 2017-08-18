function cd=aerofoilcd(X)
% cd=aerofoilcd(X)
%
% This function takes a scalar input X and 
% returns a scalar output of drag of an aerofoil
% with shape parameter X (rounded to nearest 0.01).
%
% Inputs: X scalar in range 0 1
%
% Outputs: cd scalar
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

load cd.mat
cd=cd(round(X*100)+1);
