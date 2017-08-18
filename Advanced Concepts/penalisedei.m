function NegLogExpImp=penalisedei(x)
% NegLogConExpImp=penalisedei(x)
%
% Calculates the negative of the log of the constrained 
% expected improvement at x using a one pass penalty 
% when a constraint is violated
%
% Inputs:
%	x - 1 x k vetor of of design variables
%
% Global variables used:
%   ObjectiveInfo - structure
%	ConstraintInfo - structured cell array
%   ConstraintInfo{i}.Penalty - scalar penalty
%
% Outputs:
%	NegLogConExpImp - scalar -log(E[I(x)])+P
%
% Calls:
%   predictor.m
%
% Copyright 2007 A I J Forrester
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
global ObjectiveInfo
global ConstraintInfo

% Calculate penalty for each constraint
ModelInfo=ConstraintInfo;
ModelInfo.Option='Pred';
Penalty=0;
for i=1:size(ConstraintInfo,2)
    if predictor(x)>ModelInfo{i}.ConstraintLimit;
        Penalty=Penalty+ModelInfo{i}.Penalty;
    end
end
% Penalise E[I(x)]
ModelInfo=ObjectiveInfo;
ModelInfo.Option='NegLogExpImp';
NegLogExpImp=predictor(x)+Penalty;
