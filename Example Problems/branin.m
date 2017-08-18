function f=branin(x)
% f=branin(x)
%
% This function takes a input vector x and 
% returns a scalar output of Branin's function.
%
% Inputs:
%   x -2 x 1 vector in range [0;0] [1;1]
%
% Outputs: 
%   f -scalar branin funtion value
%
% Copyright 2007 A I J Forrester and A Sobester
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

if length(x)~=2
      error('Branin''s function is for two varaibles only');
end
if x(1)<0 || x(1)>1 || x(2)<0 || x(2)>1
      error('Variable outside of range - use x in {0,1}');
end

X1 = 15*x(1)-5;
X2 = 15*x(2);
a = 1;
b = 5.1/(4*pi^2);
c = 5/pi;
d = 6;
e = 10;
ff = 1/(8*pi); 
f = (a*( X2 - b*X1^2 + c*X1 - d )^2 + e*(1-ff)*cos(X1) + e)+5*x(1);

