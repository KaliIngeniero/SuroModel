function GAOptions = gaoptimset(field1,val1,field2,val2)
% GAOptions = gaoptimset('field1',val1,'field2',val2,'field3',val3)
%
% This is a very simple version of Matlab's gaoptimset.m to work with the
% "Engineering Design via Surrogate Modelling" GA
%
% Inputs:
%   field1 - string, must be either 'PopulationSize' or 'Generations'
%   val1 - scalar value for field1
%   field2 - string, must be either 'PopulationSize' or 'Generations'
%   val2 - scalar value for field2
%
% Outputs:
%	GAOptions - structure to be passed to ga.m

%
% Copyright 2008 A I J Forrester
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

disp('Calling "Engineering Design via Surrogate Modelling" gaoptimset.')

if nargin==2
    GAOptions=struct(field1,val1);
elseif nargin==4
    GAOptions=struct(field1,val1,field2,val2);
else
    error('Incorrect number of input arguments (type help gaoptimset)')
end

