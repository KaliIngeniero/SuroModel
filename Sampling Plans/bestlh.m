function X = bestlh(n,k, Population, Iterations)
% Generates an optimized Latin hypercube by optimizing the Morris-Mitchell
% criterion for a range of exponents and plots the first two dimensions of
% the current hypercube throughout the optimization process.
%
% Inputs:
%       n - number of points required
%       k - number of design variables
%       Population - number of individuals in the evolutionary operation
%                    optimizer
%       Iterations - number of generations the evolutionary operation
%                    optimizer is run for
%       Note: high values for the two inputs above will ensure high quality
%       hypercubes, but the search will take longer.
%
% Output:
%       X - optimized Latin hypercube
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

if k<2
	error('Latin hypercubes are not defined for k<2');
end

% List of qs to optimize Phi_q for
q = [1 2 5 10 20 50 100];

% Set the distance norm to rectangular for a faster search. This can be 
% changed to p=2 if the Euclidean norm is required.
p = 1;

% We start with a random Latin hypercube
XStart = rlh(n,k,1);

% For each q optimize Phi_q
for i=1:length(q)
    fprintf('Now optimizing for q=%d...\n',q(i));
    X3D(1:n,1:k,i) = mmlhs(XStart, Population, Iterations, q(i));
end

% Sort according to the Morris-Mitchell criterion
Index = mmsort(X3D,p);

fprintf('Best lh found using q=%d...\n',q(Index(1)));

% And the Latin hypercube with the best space-filling properties is...
X = X3D(:,:,Index(1));

% Plot the projections of the points onto the first two dimensions
plot(X(:,1),X(:,2),'ro');drawnow;

title(strcat('Morris-Mitchell optimum plan found using q=',...
    num2str(q(Index(1)))));

xlabel('x_1'); ylabel('x_2');