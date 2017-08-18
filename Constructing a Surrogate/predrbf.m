function y = predrbf(x)
% Calculates the value of a Radial Basis Function surrogate model at x,
% where all model information (including the parameters, if needed for the
% type of basis function being used) are specified in the global variable
% ModelInfo.
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

global ModelInfo

d = zeros(1,size(ModelInfo.X,1));

for k=1:size(ModelInfo.X,1)
   d(k) = norm(x-ModelInfo.X(k,:),2);
end


phi = [];
for k=1:size(ModelInfo.X,1)
    if isfield(ModelInfo,'Sigma')
        phi(k) = basis(ModelInfo.Code,d(k),ModelInfo.Sigma);
    else
        phi(k) = basis(ModelInfo.Code,d(k));
    end
end
    

y = phi*ModelInfo.Weights;