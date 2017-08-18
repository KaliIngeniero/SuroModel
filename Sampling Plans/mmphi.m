function Phiq = mmphi(X, q, p)
% Calculates the sampling plan quality criterion of Morris and Mitchell.
%
% Inputs:
%       X - sampling plan
%       q - exponent used in the calculation of the metric
%       p - the distance metric to be used (p=1 rectangular - default, p=2
%           Euclidean)
%
% Output:
%       Phiq - sampling plan `space-fillingness' metric
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

% Assume defaults if arguments list incomplete
if ~exist('p','var')
    p = 1;
end

if ~exist('q','var')
    q = 2;
end

% Calculate the distances between all pairs of
% points (using the p-norm) and build multipli-
% city array J
[J,d] = jd(X,p);

% The sampling plan quality criterion
Phiq = sum(J.*(d.^(-q)))^(1/q);