function X = fullfactorial(q,Edges)
% Generates a full factorial sampling plan in the unit cube
%
% Inputs:
%       q - k-vector containing the number of points along each dimension
%       Edges - if Edges=1 the points will be equally spaced from edge to
%               edge (default), otherwise they will be in the centres of 
%               n = q(1)*q(2)*...q(k) bins filling the unit cube.
%
% X - full factorial sampling plan
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

if nargin < 2, Edges=1; end

if min(q) < 2
	error('You must have at least two points per dimension.');
end

% Total number of points in the sampling plan
n = prod(q);

% Number of dimensions
k = length(q);

%Pre-allocate memory for the sampling plan
X = zeros(n,k);

%Additional phantom element
q(k+1)=1;

for j=1:k
		if Edges==1
			one_d_slice = (0:1/(q(j)-1):1);
		else
			one_d_slice = (1/q(j)/2:1/q(j):1);
		end
		
		column = [];
		
		while length(column) < n
			for l=1:q(j)
				column = [column; ones(prod(q(j+1:k)),1)*one_d_slice(l)];
			end
		end
		
		X(:,j) = column;
end
