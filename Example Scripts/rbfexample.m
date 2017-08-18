% This script demonstrates the process of fitting a radial basis function
% surrogate
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

addpath('../Sampling Plans');
addpath('../Example Problems');
addpath('../Constructing a Surrogate');

clear global
global ModelInfo

% Sampling plan
ModelInfo.X = bestlh(6,2,25,25);


% Compute objective function values - in this case using
% the dome.m test function
for i=1:size(ModelInfo.X,1)
    ModelInfo.y(i) = dome(ModelInfo.X(i,:));
end
% y must be a column vector
ModelInfo.y = ModelInfo.y';

% Basis function type:
ModelInfo.Code = 4;

% Estimate model parameters
rbf

% Plot the surrogate
x=(0:0.025:1);

for i=1:length(x)
    for j=1:length(x)
        M(j,i) = predrbf([x(i) x(j)]);
    end
end

pcolor(x,x,M)
shading interp