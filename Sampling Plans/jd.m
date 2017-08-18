function [J,distinct_d] = jd(X,p)
% Computes the distances between all pairs of points in a sampling plan X
% using the p-norm, sorts  them  in  ascending order and removes multiple 
% occurences. 
%
% Inputs:
%       X - sampling plan being evaluated
%       p - distance norm (p=1 rectangular - default, p=2 Euclidean)
%
% Outputs:
%       J - multiplicity array (that is, the number of pairs separated by 
%           each distance value).
%       distinct_d - list of distinct distance values
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

if ~exist('p','var')
    p = 1;
end

% Number of points in the sampling plan
n = size(X,1);

% Compute the distances between all pairs of points
d = zeros(1,n*(n-1)/2);

for i = 1:n-1
    for j = i+1:n
            % Distance metric: p-norm
            d((i-1)*n-(i-1)*i/2+j-i) = norm(X(i,:)-X(j,:),p);	    
    end
end

% Remove multiple occurences
distinct_d = unique(d);


% Pre-allocate memory for J
J = zeros(size(distinct_d));

% Generate multiplicity array
for i = 1:length(distinct_d)
    % J(i) will contain the number of pairs separated
    % by the distance distinct_d(i)
    J(i) = sum(ismember(d,distinct_d(i)));
end