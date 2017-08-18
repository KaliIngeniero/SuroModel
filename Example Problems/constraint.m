function [c]=constraint(x)
% c=constraint(x)
%
% This function takes a input vector x and 
% returns a scalar output of the product constraint.
%
% Inputs:
%   x -2 x 1 vector in range [0;0] [1;1]
%
% Outputs: 
%   f -scalar constraint value
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

c=prod(x);

