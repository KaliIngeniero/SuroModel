function [f]=hepred(x)
%[f]=hepred(x)
%
% Calculates Hessian enhanced kriging prediction
% Inputs:
%	x - k x 1 vector of location of prediction
%	ModelInfo.X - n x n matrix of sample locations
%	ModelInfo.y - 1 x (2k+1)n vector of observed data
%   ModelInfo.Theta - model-parameter
%	ModelInfo.U - Choleski factorisation of PsiDot
%
% Outputs:
%	f - gradient enhanced kriging prediction
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

X=ModelInfo.X;
y=ModelInfo.y;
if size(y,1)==1
y=y';
end
theta=ModelInfo.theta;
U=ModelInfo.U;
n=length(X);


one=[ones(n,1);zeros(2*n,1)];
mu=(one'*(U\(L\y)))/(one'*(U\(L\one)));
psi=ones(n+2*k*n,1);
for i=1:n
	psi(i)=exp(-sum(theta.*abs(X(i,:)-x).^2));
	psi(n+i)=2*sum(theta.*(x-X(i,:)))*exp(-sum(theta.*abs(X(i,:)-x).^2));
	psi(2*n+i)=(-2*sum(theta)+sum(2.*theta.*(x-X(i,:)))^2)*exp(-sum(theta.*abs(X(i,:)-x).^2));
end
f=mu+psi'*(U\(U'\(y-one*mu)));
