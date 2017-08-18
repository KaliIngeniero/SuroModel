% Constrained spring optimization using max{E[I(x)n F(x)]}
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
addpath('../Advanced Concepts');
addpath('../Exploring and Expoiting a Surrogate');

rand('state',0)
randn('state',0)
global ModelInfo
global ObjectiveInfo
global ConstraintInfo
% Number of variables
k=3;
% Number of sample points
n=10;
% Create sampling plan
ObjectiveInfo.X=bestlh(n,k,20,10);
ConstraintInfo{1}.X=ObjectiveInfo.X;
ConstraintInfo{2}.X=ObjectiveInfo.X;
% Run 'experiments' to get observed data
for i=1:n
	ObjectiveInfo.y(i,1)=-1*springcycles(ObjectiveInfo.X(i,:));
	ConstraintInfo{1}.y(i,1)=buckling(ConstraintInfo{1}.X(i,:));
	ConstraintInfo{2}.y(i,1)=shearsafetyfact(ConstraintInfo{2}.X(i,:));
end
% Constraint limits
ConstraintInfo{1}.ConstraintLimit=0;
ConstraintInfo{2}.ConstraintLimit=0;

% Figure for progress plot
figure
for I=1:20
    % Use only successful sample points in objective function model
    ObjectiveInfo.X=ObjectiveInfo.X(find(~isnan(ObjectiveInfo.y)),:);
    ObjectiveInfo.y=ObjectiveInfo.y(find(~isnan(ObjectiveInfo.y)));
    % Set upper and lower bounds for search of log theta
    UpperTheta=ones(1,k).*2;
    LowerTheta=ones(1,k).*-3;
    % Tune kriging model of objective
    ModelInfo=ObjectiveInfo;
    [ObjectiveInfo.Theta,MaxLikelihood] = ga(@likelihood,k,[],[],[],[],LowerTheta,UpperTheta);
    [NegLnLike,ObjectiveInfo.Psi,ObjectiveInfo.U]=likelihood(ObjectiveInfo.Theta);
    for i=1:size(ConstraintInfo,2)
        % Tune kriging model of constraint
        ModelInfo=ConstraintInfo{i};
        [ConstraintInfo{i}.Theta,MaxLikelihood] = ga(@likelihood,k,[],[],[],[],LowerTheta,UpperTheta);
        [NegLnLike,ConstraintInfo{i}.Psi,ConstraintInfo{i}.U]=likelihood(ConstraintInfo{i}.Theta);
    end
    % Search constrained expeccted improvment
    options=gaoptimset('PopulationSize',100);
    [OptVar,OptEI]=ga(@constrainedei,k,[],[],[],[],[0 0 0],[1 1 1],[],options);
    % Add infill point
    ObjectiveInfo.X(end+1,:)=OptVar;
    ObjectiveInfo.y(end+1)=-1*springcycles(OptVar);
    ConstraintInfo{1}.X(end+1,:)=OptVar;
    ConstraintInfo{2}.X(end+1,:)=OptVar;
    ConstraintInfo{1}.y(end+1)=buckling(OptVar);
    ConstraintInfo{2}.y(end+1)=shearsafetyfact(OptVar);
    % Plot progress
    close all
    figure
    hold on
    for i=1:length(ObjectiveInfo.y)
        plot(i,-ObjectiveInfo.y(i),'ko')
        if buckling(ObjectiveInfo.X(i,:))>0 || shearsafetyfact(ObjectiveInfo.X(i,:)) >0
            plot(i,-ObjectiveInfo.y(i),'kx')
        end
    end
    set(gca,'YScale','log')
    drawnow
end

xlabel('sample points','FontSize',14)
ylabel('cycles','FontSize',14)
set(gca,'FontSize',14)

