function f=gepred(x)
%f=gepred(x)
%
% Calculates gradient enhanced kriging prediction
% Inputs:
%	x - k x 1 vector of location of prediction
%	ModelInfo.X - n x n matrix of sample locations
%	ModelInfo.y - 1 x (2k+1)n vector of observed data
%   ModelInfo.Theta - hyper-parameter
%	ModelInfo.U - Cholesky factorisation of PsiDot
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

global ModelInfo

X=ModelInfo.X;
y=ModelInfo.y;
n=size(X,1);
k=size(X,2);
theta=10.^ModelInfo.Theta;
U=ModelInfo.U;

one=[ones(n,1);zeros(k*n,1)];
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
psi=ones(n+k*n,1);
 for i=1:n
 	psi(i)=exp(-sum(theta.*abs(X(i,:)-x).^2));
     for j=1:k
         psi(j*n+i)=2*theta(j).*(x(j)-X(i,j))*exp(-sum(theta.*(X(i,:)-x).^2));
     end
 end
f=mu+psi'*(U\(U'\(y-one*mu)));
