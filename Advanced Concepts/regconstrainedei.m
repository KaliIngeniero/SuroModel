function NegLogConExpImp=regconstrainedei(x)
% NegLogConExpImp=constrainedei(x)
%
% Calculates the negative of the log of the 
% constrained expected improvement at x using Kriging regression
%
% Inputs:
%	x - 1 x k vetor of of design variables
%
% Global variables used:
%   ObjectiveInfo - structure
%	ConstraintInfo - structured cell array
%
% Outputs:
%	NegLogConExpImp - scalar -log(E[I(x)]P[F(x)])
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

% Calculate unconstrained E[I(x)]
ModelInfo=ObjectiveInfo;
ModelInfo.Option='NegLogExpImp';
NegLogExpImp=regpredictor(x);

% Calculate P[F(x)] for each constraint
for i=1:size(ConstraintInfo,2)
    ModelInfo=ConstraintInfo{i};
    ModelInfo.Option='NegProbImp';
    NegProbImp(i)=regpredictor(x);
end

% Calculate E[I(x)]P[F(x)] (add realmin before taking logs)
NegLogConExpImp=-(-NegLogExpImp+sum(log10(-NegProbImp+realmin))); % changed from "NegLogConExpImp=-(-NegLogExpImp+sum(log10(-NegProbImp+1e-50)));" (September 2009)

