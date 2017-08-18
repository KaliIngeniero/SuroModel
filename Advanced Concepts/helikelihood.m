function [NegLnLike,Cov,U]=helikelihood(x)
%[NegLnLike,Cov,U]=helikelihood(x)
%
% Calculates the negative of the concentrated log-likelihood
% Inputs:
%	x - vetor of paramters, either 1 x k specifies theta only (p=2 by default),
%	1 x 2k specifies theta and p
%	1 x k+1 specifies theta and lambda
%	1 x 2k+1 specifies theta, p and lambda
%	ModelInfo.X - k x n matrix of sample locations
%	ModelInfo.y - 1 x n vector of observed data
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
theta=10.^x;
n=size(X,1);
k=size(X,2);
one=[ones(n,1); zeros(2*n*k,1)];

%%% build upper half of correlation matrix %%%
Psi=zeros(n,n); % initialise to zeros (saves time)
for i=1:n
	for j=i+1:n
		Psi(i,j)=exp(-sum(theta.*abs(X(i,:)-X(j,:)).^2)); % build top right corner
	end
end
Psi=Psi+Psi'+eye(n)+eye(n).*lambda;


CovXdY=zeros(n,n);
CovdYdY=zeros(n,n);
CovXdYdY=zeros(n,n);
CovX1dYdY=zeros(n,n);
CovX2dYdY=zeros(n,n);
CovXdYdYdY=zeros(n,n);
CovdYdYdYdY=zeros(n,n);
for i=1:n
	for j=i+1:n		
		CovXdY(i,j)=-2*sum(theta.*(X(i,:)-X(j,:)))*Psi(i,j);
		CovdYdY(i,j)=(2*sum(theta)-sum(2.*theta.*(X(i,:)-X(j,:)))^2)*Psi(i,j);
		
		CovX1dYdY(i,j)=(-2*sum(theta)+sum(2.*theta.*(X(i,:)-X(j,:)))^2)*Psi(i,j);
		CovX2dYdY(i,j)=(-2*sum(theta)+sum(2.*theta.*(X(i,:)-X(j,:)))^2)*Psi(i,j);
		
		CovXdYdYdY(i,j)=(-12*theta^2*(X(i,:)-X(j,:))+8*theta^3*(X(i,:)-X(j,:))^3)*Psi(i,j);
		
		CovdYdYdYdY(i,j)=(12*theta^2-48*theta^3*(X(i,:)-X(j,:))^2+16*theta^4*(X(i,:)-X(j,:))^4)*Psi(i,j);
	end
end
%%% add upper and lower halves %%%
Cov=[Psi (-CovXdY+CovXdY') (CovX1dYdY+CovX1dYdY'-eye(n).*(2*sum(theta)));
(CovXdY-CovXdY') (CovdYdY+CovdYdY'+eye(n).*(2*sum(theta))) (CovXdYdYdY-CovXdYdYdY');
(CovX2dYdY+CovX2dYdY'-eye(n).*(2*sum(theta))) (-CovXdYdYdY+CovXdYdYdY') (CovdYdYdYdY+CovdYdYdYdY'+eye(n).*12*theta^2)];

%Cov=[Psi CovXdY;
%CovXdY' CovdYdY];

%%% concentrated log-likelihood calculation %%%
%%% LU decomposition %%%
[L,U]=lu(Cov);
%%% sum logs of diagonal to find ln(abs(det(Psi))) %%%
lndetPsi=sum(log(abs(diag(U))));
%%% use back-substitution of LU decomposition instead of inverse 
mu=(one'*(U\(L\y)))/(one'*(U\(L\one)));
sigma_sqr=((y-one*mu)'*(U\(L\(y-one*mu))))/(2*n);
NegLnLike=-1*(-((2*n)/2)*log(sigma_sqr) - 0.5*lndetPsi);

