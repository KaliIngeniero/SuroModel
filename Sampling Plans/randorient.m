function Bstar = randorient(k, p, xi)
% Generates a random orientation for a screening matrix
%
% Inputs:
%       k - number of design variables
%       p - number of discreet levels along each dimension
%       xi- elementery effect step length factor
%
% Output:
%       Bstar - random orientation matrix
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


% Step length
Delta = xi/(p-1);

m = k + 1;

% A truncated p-level grid in one dimension
xs = (0:1/(p-1):1-Delta);
xsl = length(xs);

% Basic sampling matrix
B = [zeros(1,k); tril(ones(k))];

% Randomization

% Matrix with +1s and -1s on the diagonal with equal probability
Dstar = diag(2*round(rand(1,k))-1);

% Random base value
xstar = xs(floor(rand(1,k)*xsl)+1);

% Permutation matrix
Pstar = zeros(k);
rp = randperm(k);
for i=1:k, Pstar(i,rp(i))=1; end

% A random orientation of the sampling matrix
Bstar = (ones(m,1)*xstar+(Delta/2)*...
    ((2*B-ones(m,k))*Dstar+ones(m,k)))*Pstar;