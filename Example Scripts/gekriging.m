% Gradient enhanced Kriging example
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

addpath('../Advanced Concepts');

rand('state',0)
randn('state',0)
global ModelInfo
% Sampling plan
X=[0;0.33;0.66;1];
for i=1:4
    % One variable function values
    f(i,1) =((X(i).*6-2).^2).*sin(X(i).*12-4);
    % One variable function gradients
    dfdx(i,1)=12*(X(i).*6-2)*sin(X(i).*12-4)+((X(i).*6-2).^2)*cos(X(i).*12-4)*12;
end
% Put into ModelInfo
ModelInfo.X=X;
ModelInfo.y=[f;dfdx];
% Maximize likelihood
[ModelInfo.Theta,MaxLikelihood] = fminbnd(@gelikelihood,-3,3)
[NegLnLike,ModelInfo.Psi,ModelInfo.U]=gelikelihood(ModelInfo.Theta)

% Plot results
Xplot=[0:0.01:1]';
for i=1:101
    % GE Kriging prediction
    pred(i)=gepred(Xplot(i));
    true(i)=((Xplot(i).*6-2).^2).*sin((Xplot(i).*6-2).*2);
end
figure
plot(Xplot,true,'k')
hold on
plot(ModelInfo.X,f,'ko')
plot(Xplot,pred,'r')
legend('f(x)','sample points','GE Kriging',2)

