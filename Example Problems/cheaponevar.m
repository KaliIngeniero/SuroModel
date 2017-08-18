function f=cheaponevar(x)
% f=cheaponevar(x)
%
% This function takes a scalar input x and 
% returns a scalar output of the cheap one-variable test function.
%
% Inputs: 
%   x - scalar in range 0 1
%
% Outputs: 
%   f - scalar funtion value
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

if length(x)~=1
      error('Too many variables - function is for one varaible only');
end
if x<0 || x>1
      error('Variable outside of range - use x in {0,1}');
end
A=0.5;
B=10;
C=-5;
D=0;
f =A.*(((x+D).*6-2).^2).*sin(((x+D).*6-2).*2)+((x+D)-0.5).*B+C;
