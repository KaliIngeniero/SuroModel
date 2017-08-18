% Constrained multi-objective spring optimization using
% max{E[I(x*)n F(x*)]}
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
addpath('../Exploring and Exploiting a Surrogate'); % spelling corected (September 2009)

rand('state',0)
randn('state',0)
global ModelInfo
global ObjectiveInfo
global ConstraintInfo

%Create sampling plan
n=15;
k=3;
ObjectiveInfo{1}.X=bestlh(n,k,20,10)
ObjectiveInfo{2}.X=ObjectiveInfo{1}.X;
nConstraints=2;
for i=1:nConstraints
    ConstraintInfo{i}.ConstraintLimit=0;
    ConstraintInfo{i}.X=ObjectiveInfo{1}.X;
end
for i=1:n
    ObjectiveInfo{1}.y(i,1)=springmass(ObjectiveInfo{1}.X(i,:));
    ObjectiveInfo{2}.y(i,1)=-1*(springcycles(ObjectiveInfo{2}.X(i,:)));
    ConstraintInfo{1}.y(i,1)=buckling(ConstraintInfo{1}.X(i,:));
    ConstraintInfo{2}.y(i,1)=shearsafetyfact(ConstraintInfo{2}.X(i,:));
end

figure
for I=1:60
    ConstraintInfo{1}.X=ConstraintInfo{1}.X(find(~isnan(ObjectiveInfo{2}.y)),:);
    ConstraintInfo{2}.X=ConstraintInfo{2}.X(find(~isnan(ObjectiveInfo{2}.y)),:);
    ConstraintInfo{1}.y=ConstraintInfo{1}.y(find(~isnan(ObjectiveInfo{2}.y)));
    ConstraintInfo{2}.y=ConstraintInfo{2}.y(find(~isnan(ObjectiveInfo{2}.y)));
    ObjectiveInfo{1}.X=ObjectiveInfo{1}.X(find(~isnan(ObjectiveInfo{2}.y)),:);
    ObjectiveInfo{2}.X=ObjectiveInfo{2}.X(find(~isnan(ObjectiveInfo{2}.y)),:);
    ObjectiveInfo{1}.y=ObjectiveInfo{1}.y(find(~isnan(ObjectiveInfo{2}.y)));
    ObjectiveInfo{2}.y=ObjectiveInfo{2}.y(find(~isnan(ObjectiveInfo{2}.y)));
    % tune kriging models of objectives
    options=gaoptimset('PopulationSize',20,'Generations',10);
    for i=1:2
        ModelInfo=ObjectiveInfo{i};
        ObjectiveInfo{i}.Theta = ga(@likelihood,k,[],[],[],[],ones(k,1).*-3,ones(k,1).*3,[],options);
        [NegLnLike,ObjectiveInfo{i}.Psi,ObjectiveInfo{i}.U]=likelihood(ObjectiveInfo{i}.Theta);
    end
    % tune kriging models of constraints
    for i=1:nConstraints
        ModelInfo=ConstraintInfo{i};
        ConstraintInfo{i}.Theta = ga(@likelihood,k,[],[],[],[],ones(k,1).*-3,ones(k,1).*3,[],options);
        [NegLnLike,ConstraintInfo{i}.Psi,ConstraintInfo{i}.U]=likelihood(ConstraintInfo{i}.Theta);
    end
    
    y1temp=ObjectiveInfo{1}.y;
    y2temp=ObjectiveInfo{2}.y;
    Xtemp=ObjectiveInfo{2}.X;
    for i=1:length(y1temp)
        for j=1:nConstraints
            if ConstraintInfo{j}.y(i)>ConstraintInfo{j}.ConstraintLimit
                y1temp(i)=nan;
                y2temp(i)=nan;
            end
        end    
    end
    Xtemp=Xtemp(find(~isnan(y2temp)),:);
    y1temp=y1temp(find(~isnan(y1temp)));
    y2temp=y2temp(find(~isnan(y2temp)));
    % find Pareto set
    [a,b]=sort(y1temp);
    PX(1,1:k)=Xtemp(b(1),1:k);
    Py1(1)=y1temp(b(1));
    Py2(1)=y2temp(b(1));
    Pnum=1;
    for i=2:length(y1temp)
        if y2temp(b(i))<=Py2(end)
            Pnum=Pnum+1;
            PX(Pnum,1:k)=Xtemp(b(i),1:k);
            Py1(Pnum)=y1temp(b(i));
            Py2(Pnum)=y2temp(b(i));
        end
    end    

    plot(Py1,(-Py2),'ko')
    set(gca,'YScale','log')
    title('Tradeoff')
    xlabel('Obj 1')
    ylabel('Obj 2')
    axis square
    drawnow
    pause(1)
    % search -log(E[I(x*)n F(x*)])
    options=gaoptimset('PopulationSize',50,'Generations',20);
    [VarOpt,EIOpt] = ga(@constrainedmultiei,k,[],[],[],[],zeros(k,1),ones(k,1),[],options)
    ObjectiveInfo{1}.X(end+1,:)=VarOpt;
    ObjectiveInfo{2}.X(end+1,:)=VarOpt;
    for i=1:nConstraints
        ConstraintInfo{i}.X(end+1,:)=VarOpt;
    end
    ObjectiveInfo{1}.y(end+1)=springmass(ObjectiveInfo{1}.X(end,:));
    ObjectiveInfo{2}.y(end+1,1)=-1*(springcycles(ObjectiveInfo{2}.X(end,:)));
    ConstraintInfo{1}.y(end+1,1)=buckling(ConstraintInfo{1}.X(end,:));
    ConstraintInfo{2}.y(end+1,1)=shearsafetyfact(ConstraintInfo{2}.X(end,:));
    clear PX Py1 Py2
end

