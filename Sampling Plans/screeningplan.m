function X = screeningplan(k, p, xi, r)
% Generates a Morris screening plan with a specified number of elementary
% effects for each variable.
%
% Inputs:
%       k - number of design variables
%       p - number of discreet levels along each dimension
%       xi- elementery effect step length factor
%       r - number of random orientations (=number of elementary effects
%           per variable).
%
% Output:
%       X - screening plan built within a [0,1]^k box
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


X = [];
for i=1:r
    X = [X; randorient(k,p,xi)];
end