% One variable co-Kriging example (similar to Figure 8.1)
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

addpath('../Example Problems');
addpath('../Advanced Concepts');

rand('state',0)
randn('state',0)
global ModelInfo
% Expensive points
ModelInfo.Xe=[0; 0.4; 0.6;1];
% Cheap points
ModelInfo.Xc=[0.1;0.2;0.3;0.5;0.7;0.8;0.9;0;0.4;0.6;1];
k=1;
% Calulate expensive observations
for i=1:size(ModelInfo.Xe,1)
ModelInfo.ye(i,1)=onevar(ModelInfo.Xe(i));
end
% Calculate cheap observations
for i=1:size(ModelInfo.Xc,1)
ModelInfo.yc(i,1)=cheaponevar(ModelInfo.Xc(i));
end
% Optimize cheap model paramters
ModelInfo.Thetac=fminbnd(@likelihoodc,-3,3);
% Optimise difference model paramters
Params=ga(@likelihoodd,k+1,[],[],[],[],[-3 -5],[3 5]);
ModelInfo.Thetad=Params(1:k);
ModelInfo.rho=Params(k+1);
% Construct covariance matrix
buildcokriging
% Make predictions in range 0,1
Xplot=0:0.01:1;
ModelInfo.Option='Pred';
for i=1:101
    pred(i)=cokrigingpredictor(Xplot(i));
    truee(i)=onevar(Xplot(i));
    truec(i)=cheaponevar(Xplot(i));
end

plot(Xplot,truee,'k','LineWidth',2)
hold on
plot(Xplot,truec,'k--','LineWidth',2)
plot(ModelInfo.Xe,ModelInfo.ye,'ko')
plot(ModelInfo.Xc,ModelInfo.yc,'ko')
plot(Xplot,pred,'r')
legend('f_e','f_c','y_e','y_c','co-Kriging',2)
xlabel('x','FontSize',14)
ylabel('f','FontSize',14)
set(gca,'FontSize',14)

