function f = basis(varargin)
% USAGE basis(Code,r,[Sigma])
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

if nargin==0
    error('Basis requires input arguments.');
elseif nargin==1
    % We default to cubic splines
    Code = 2;
    r = varargin{1};
elseif nargin==2
    % If only two arguments given, only non-parametric bases can be
    % calculated - if first argument greater, we default to 3.
    Code = min([varargin{1},3]);
    r = varargin{2};
else
    % Defaults to inverse multi-quadric if first argument >6
    Code = min([varargin{1},6]);
    r = varargin{2};
    Sigma = varargin{3};
end



switch Code
case 1
   % Linear function
   f = r;
case 2
   % Cubic
   f = r^3;
case 3
   % Thin plate spline
   if r < 1e-200
      f=0;
    else   
      f = r^2*log(r);
   end
case 4
   % Gaussian
   f = exp(-(r^2)/(2*Sigma^2));
case 5
   % Multi-quadric
   f = (r^2+Sigma^2)^0.5;
case 6
   % Inverse Multi-Quadric
   f = (r^2+Sigma^2)^(-0.5);
end