function f=pred(x)
% f=pred(x)
%
% Calculates a Kriging prediction at x
%
% Inputs:
%	x - 1 x k vetor of design variables
%
% Global variables used:
%	ModelInfo.X - n x k matrix of sample locations
%	ModelInfo.y - n x 1 vector of observed data
%   ModelInfo.Theta - 1 x k vector of log(theta)
%   ModelInfo.U - n x n Cholesky factorisation of Psi
%
% Outputs:
%	f - scalar kriging prediction
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

global ModelInfo
% extract variables from data structure 
% slower, but makes code easier to follow
X=ModelInfo.X;
y=ModelInfo.y;
theta=10.^ModelInfo.Theta;
p=2;  % added p definition (February 10)
U=ModelInfo.U;

% calculate number of sample points
n=size(X,1);

% vector of ones
one=ones(n,1);

% calculate mu
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));

% initialise psi to vector of ones
psi=ones(n,1);

% fill psi vector
for i=1:n
	psi(i)=exp(-sum(theta.*abs(X(i,:)-x).^p));
end

% calculate prediction
f=mu+psi'*(U\(U'\(y-one*mu)));
