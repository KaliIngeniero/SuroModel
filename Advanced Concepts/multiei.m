function [metric,Py1,Py2,PX]=multiei(x)
% [metric,Py1,Py2,PX]=multiei(x)
%
% Calculates the expection of f_1(x) and f_2(x) improving 
% on the Pareto front defined by ObjectiveInfo{1:2}.y 
%
% Inputs:
%	x - 1 x k vetor of design variables
%
% Global variables used:
%	ObjectiveInfo{1:2}.X - n x k matrix of sample locations for each
%	objective
%	ObjectiveInfo{1:2}.y - n x 1 vector of observed data for each objective
%   ObjectiveInfo{1:2}.Theta - 1 x k vector of log(theta) for each
%   objective
%   ObjectiveInfo{1:2}.U - n x n Cholesky factorisation of Psi  for each
%   objective
%   ModelInfo.Option - string: 'NegLogExpImp' or 'NegProbImp'
%
% Outputs:
%	metric - either -log(E[I(x*)]) -P[I(x*)], determined by
%	ModelInfo.option
%   Py1,Py2 - non-dominated objective function values on the Pareto front
%   PX - locations of non-dominated solutions
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
%



global ObjectiveInfo
global ConstraintInfo
global ModelInfo
X=ObjectiveInfo{1}.X;
y1=ObjectiveInfo{1}.y;
y2=ObjectiveInfo{2}.y;
k=size(ObjectiveInfo{1}.X,2);
Option=ModelInfo.Option;

% find points which satisfy constraint (if present)
y1temp=y1;
y2temp=y2;
Xtemp=X;
if exist('ConstraintInfo','var')==1 % changed from "if isfield(ConstraintInfo{1},'ConstraintLimit')==1"  (September 2009)
    for i=1:length(y1)
        for j=1:size(ConstraintInfo,2)
            if ConstraintInfo{j}.y(i)>ConstraintInfo{j}.ConstraintLimit
                y1temp(i)=nan;
                y2temp(i)=nan;
            end
        end
    end
    Xtemp=Xtemp(find(~isnan(y2temp)),:);
    y1temp=y1temp(find(~isnan(y1temp)));
    y2temp=y2temp(find(~isnan(y2temp)));
    
end
%find Pareto set
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

% prediction of each objective
ModelInfo=ObjectiveInfo{1};
ModelInfo.Option='Pred';
pred1=predictor(x);
ModelInfo=ObjectiveInfo{2};
ModelInfo.Option='Pred';
pred2=predictor(x);
% RMSE of each objective
ModelInfo=ObjectiveInfo{1};
ModelInfo.Option='RMSE';
s1=predictor(x);
ModelInfo=ObjectiveInfo{2};
ModelInfo.Option='RMSE';
s2=predictor(x);

% probability of improvement calculation
PITerm1=(0.5+0.5*erf((1/(2^0.5))*((Py1(1)-pred1)/s1)));
PITerm3=(1-(0.5+0.5*erf((1/(2^0.5))*((Py1(end)-pred1)/s1)))) ...
* (0.5+0.5*erf((1/(2^0.5))*((Py2(end)-pred2)/s2)));
if Pnum>1
    for I=1:length(Py1)-1
        PITerm2calc(I)=((0.5+0.5*erf((1/(2^0.5))*((Py1(I+1)-pred1)/s1)))...
        -(0.5+0.5*erf((1/(2^0.5))*((Py1(I)-pred1)/s1))))...
        *(0.5+0.5*erf((1/(2^0.5))*((Py2(I+1)-pred2)/s2)));
    end
    PITerm2=sum(PITerm2calc);
	PI=(PITerm1+PITerm2+PITerm3);
else
    PI=(PITerm1+PITerm3);
end

if strcmp(Option,'NegLogExpImp')==1
	% Ybar1 calculation
	Ybar1Term1=pred1*(0.5+0.5*erf((1/(2^0.5))*((Py1(1)-pred1)/s1)))...
	-s1*(1/((2*pi)^0.5))*exp(-0.5*((Py1(1)-pred1)^2/s1^2));

    Ybar1Term3=(pred1*(0.5+0.5*erf((1/(2^0.5))*((pred1-Py1(end))/s1)))... % changed this term (March 2010)
	+s1*(1/((2*pi)^0.5))*exp(-0.5*((pred1-Py1(end))^2/s1^2))) ...  
	*(0.5+0.5*erf((1/(2^0.5))*((Py2(end)-pred2)/s2)));

    if Pnum>1
        for I=1:length(Py1)-1
            Ybar1Term2calc(I)=((pred1*(0.5+0.5*erf((1/(2^0.5))*((Py1(I+1)-pred1)/s1)))...
            -s1*(1/((2*pi)^0.5))*exp(-0.5*((Py1(I+1)-pred1)^2/s1^2))) ...
            -(pred1*(0.5+0.5*erf((1/(2^0.5))*((Py1(I)-pred1)/s1))) ...
            -s1*(1/((2*pi)^0.5))*exp(-0.5*((Py1(I)-pred1)^2/s1^2))))...
            *(0.5+0.5*erf((1/(2^0.5))*((Py2(I+1)-pred2)/s2)));
        end
        Ybar1Term2=sum(Ybar1Term2calc);
        Ybar1=(Ybar1Term1+Ybar1Term2+Ybar1Term3)/PI;
    else
        Ybar1=(Ybar1Term1+Ybar1Term3)/PI;
    end
    % Ybar2 calculation - now integrates from last point to first (March
    % 2010)
    
	Ybar2Term1=pred2*(0.5+0.5*erf((1/(2^0.5))*((Py2(end)-pred2)/s2)))...  % changed this term (March 2010)
	-s2*(1/((2*pi)^0.5))*exp(-0.5*((Py2(end)-pred2)^2/s2^2));
    
    Ybar2Term3=(pred2*(0.5+0.5*erf((1/(2^0.5))*((pred2-Py2(1))/s2)))... % changed this term (March 2010)
	+s2*(1/((2*pi)^0.5))*exp(-0.5*((pred2-Py2(1))^2/s2^2))) ...  
	*(0.5+0.5*erf((1/(2^0.5))*((Py1(1)-pred1)/s1)));

    if Pnum>1
        for I=length(Py2):-1:2
            Ybar2Term2calc(I)=((pred2*(0.5+0.5*erf((1/(2^0.5))*((Py2(I-1)-pred2)/s2)))...
            -s2*(1/((2*pi)^0.5))*exp(-0.5*((Py2(I-1)-pred2)^2/s2^2))) ...
            -(pred2*(0.5+0.5*erf((1/(2^0.5))*((Py2(I)-pred2)/s2))) ...
            -s2*(1/((2*pi)^0.5))*exp(-0.5*((Py2(I)-pred2)^2/s2^2))))...
            *(0.5+0.5*erf((1/(2^0.5))*((Py1(I-1)-pred1)/s1)));
        end
        Ybar2Term2=sum(Ybar2Term2calc);
        Ybar2=(Ybar2Term1+Ybar2Term2+Ybar2Term3)/PI;
    else
        Ybar2=(Ybar2Term1+Ybar2Term3)/PI;
    end
    % find closest point on front
    for i=1:Pnum
        dist(i)=sqrt((Ybar1-Py1(i))^2+(Ybar2-Py2(i))^2);
    end
    [a,b]=sort(dist);
    % expected improvement calculation 
    if PI==0
        EI=0;
    else
    EI=PI*a(1);
    end
end
if strcmp(Option,'NegLogExpImp')==1
metric=-log10(EI+realmin); % changed from "metric=-log10(EI+1e-50);" (September 2009)
else
metric=-PI;
end
ModelInfo.Option=Option; % added (September 2009)







