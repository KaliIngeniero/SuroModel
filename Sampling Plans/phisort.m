function Index = phisort(X3D,q,p)
% Ranks sampling plans according to the Morris-Mitchell numerical quality
% metric. Note: similar to mmsort, which uses the definition of `space-
% fillingness' as a basis for the ranking.
%
% Inputs:
%       X3D - three-dimensional array containing the sampling plans to be
%       ranked.
%       p - the distance metric to be used (p=1 rectangular - default, p=2
%       Euclidean)
%       q - exponent used in the calculation of the metric
%
% Output:
%       Index - index array containing the ranking
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

if ~exist('q','var')
    q = 2;
end

% Pre-allocate memory
Index = (1:size(X3D,3));

% Bubble-sort
swap_flag = 1;

while swap_flag==1
    swap_flag = 0;
    i=1;
    while i<=length(index)-1
        if mmphi(X3D(:,:,Index(i)),q,p) > mmphi(X3D(:,:,Index(i+1)),q,p)
            buffer = Index(i);
            Index(i) = Index(i+1);
            Index(i+1) = buffer;
            swap_flag = 1;
        end
        i = i + 1;
    end
end