function f = rastrigin(x)
% f = rastrigin(x)
%
% This function takes a input vector x and 
% returns a scalar output of the Rastrigin function:
% a difficult search problem used to test GA
% Global optimum is 0 at x = [0 0 .... 0]
%
% Inputs:
%   x - vector of any length in range {-5.12 5.12}
%
% Outputs: 
%   f -scalar Rastrigin funtion value
%

%
% Copyright 2009  A Sobester
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

f =10*length(x)+sum(x.^2-10*cos(2*pi*x));