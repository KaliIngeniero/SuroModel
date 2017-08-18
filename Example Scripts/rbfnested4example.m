% This script demonstrates the process of fitting a radial basis function
% surrogate to four-dimensional data and using nested4.m to generate a
% nested plot of this surrogate
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

addpath('../Example Problems');
addpath('../Sampling Plans');
addpath('../Advanced Concepts');

clear global
global ModelInfo

% Sampling plan
ModelInfo.X = bestlh(50,4,10,10);

% Compute objective function values - in this case using the dome.m test
% function - you would insert your own objective function here
for i=1:size(ModelInfo.X,1)
    ModelInfo.y(i) = dome(ModelInfo.X(i,:));
end

% y must be a column vector
ModelInfo.y = ModelInfo.y';

% Select basis function type:
ModelInfo.Code = 4;

% Estimate model parameters
rbf

% Plot the surrogate if the model was successfuly fitted
if ModelInfo.Success == 1
    nested4([1 2 3 4], [10 10], [ 0 0 0 0; 1 1 1 1],...
        {'x_1','x_2','x_3','x_4'}, @predrbf,...
        30, 0, 1, 1)
else
    display('Could not fit model. Try a different basis function.')
end
    