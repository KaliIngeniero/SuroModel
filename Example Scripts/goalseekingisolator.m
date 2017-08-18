% Vibration isolator feasibility search using Kriging goal seeking
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

addpath('../Sampling Plans');
addpath('../Example Problems');
addpath('../Constructing a Surrogate');
addpath('../Exploring and Exploiting a Surrogate');

rand('state',0)
randn('state',0)
global ModelInfo
% Create sampling plan
k=18;
n=20;
ModelInfo.X=bestlh(n,k,20,10)

% Calculate observed data
for i=1:n
    ModelInfo.y(i,:)=intersections(ModelInfo.X(i,:));
end

% Search goal
ModelInfo.Goal=0;

% Iterate until goal is attained
while min(ModelInfo.y)>ModelInfo.Goal
    % Tune kriging model of objective
    options=gaoptimset('PopulationSize',50);
    [ModelInfo.Theta,MaxLikelihood]=ga(@likelihood,k,[],[],[],[],ones(1,k).*-3,ones(1,k).*3,[],options)
    % Put Cholesky factorisation of Psi and theta into ModelInfo
    [NegLnLike,ModelInfo.Psi,ModelInfo.U]=likelihood(ModelInfo.Theta);
    
    % Find location which maximises likelhood of goal
    options=gaoptimset('PopulationSize',100);
    [OptVar,NegCondLike]=ga(@condlikelihood,2*k,[],[],[],[],[ones(1,k).*-3 zeros(1,k)],[ones(1,k).*3 ones(1,k)],[],options)
    
    % Add infill point and calculate objective value 
    ModelInfo.X(end+1,:)=OptVar(k+1:2*k);
    ModelInfo.y(end+1)=intersections(ModelInfo.X(end,:));
end
