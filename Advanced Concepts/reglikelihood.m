function [NegLnLike,Psi,U]=reglikelihood(x)
%[NegLnLike,Psi,L,U]=likelihood(x)
%
% Calculates the negative of the concentrated log-likelihood
% Inputs:
%	x - k+1 vetor of log(theta) and lambda hyper-parameters
% Global variables used:
%	ModelInfo.X - k x n matrix of sample locations
%	ModelInfo.y - 1 x n vector of observed data
%
% Outputs:
%	NegLnLike - concentrated log-likelihood *-1 for minimising
%	Psi - correlation matrix
%	L,U - lu decomposition of correlation matrix
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
theta=10.^x(1:end-1);
p=2;  % added p definition (February 10)
lambda=10.^x(end);
n=size(X,1);
one=ones(n,1);
% build upper half of correlation matrix
Psi=zeros(n,n); % initialise to zeros (saves time)
for i=1:n
	for j=i+1:n
		Psi(i,j)=exp(-sum(theta.*abs(X(i,:)-X(j,:)).^p)); % build top right corner
	end
end
% add upper and lower halves and diagonal of ones plus lambda
Psi=Psi+Psi'+eye(n)+eye(n).*lambda; 
% concentrated log-likelihood calculation
% Cholesky factorisation
[U,p]=chol(Psi);
% sum logs of diagonal to find ln(abs(det(Psi)))
LnDetPsi=2*sum(log(abs(diag(U))));
% use back-substitution of LU decomposition instead of inverse
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
SigmaSqr=((y-one*mu)'*(U\(U'\(y-one*mu))))/n;
NegLnLike=-1*(-(n/2)*log(SigmaSqr) - 0.5*LnDetPsi);
