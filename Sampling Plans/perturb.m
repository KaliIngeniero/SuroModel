function X = perturb(X, PertNum)
% Interchanges pairs of randomly chosen elements within randomly chosen
% columns of a sampling plan a number of times. If the plan is a Latin
% hypercube, the result of this operation will also be a Latin hypercube.
%
% Inputs:
%       X - sampling plan
%       PertNum - the number of changes (perturbations) to be made to X.
%
% Output:
%       X - perturbed sampling plan
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

if ~exist('PertNum','var')
	PertNum = 1;
end


[n,k] = size(X);

for  pert_count=1:PertNum
	col = floor(rand*k)+1;
	
	% Choosing two distinct random points
	el1 = 1; el2 = 1;	
	while el1==el2
		el1 = floor(rand*n)+1;
		el2 = floor(rand*n)+1;
	end

	%Swap the two chosen elements
	buffer = X(el1,col);
	X(el1,col) = X(el2,col);
	X(el2,col) = buffer;
end