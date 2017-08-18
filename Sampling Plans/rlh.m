function X = rlh(n, k, Edges)
% Generates a random Latin hypercube within the [0,1]^k hypercube.
%
% Inputs:
%       n - desired number of points
%       k - number of design variables (dimensions)
%       Edges - if Edges=1 the extreme bins will have their centres on the
%               edges of the domain, otherwise the bins will be entirely 
%               contained within the domain (default setting). 
%
% Output:
%       X - Latin hypercube sampling plan of n points in k dimensions.
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


if nargin < 3
    Edges = 0;
end

% Pre-allocate memory
X = zeros(n,k);

for i=1:k
   X(:,i) = randperm(n)';
end

if Edges == 1
    X = (X-1)/(n-1);
else
    X = (X-0.5)/n;
end