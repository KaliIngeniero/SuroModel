function [NegLnLike,PsiDot,U]=gelikelihood(x)
%[NegLnLike,PsiDot,U]=gelikelihood(x,ModelInfo)
%
% Calculates the negative of the concentrated log-likelihood
% Inputs:
%	x - vetor of paramters, either 1 x k specifies theta only or
%	1 x k+1 specifies theta and lambda
%	ModelInfo.X - n x n matrix of sample locations
%	ModelInfo.y - 1 x (2k+1)n vector of observed data
%
% Outputs:
%	NegLnLike - concentrated log-likelihood *-1 for minimising
%	PsiDot - differentiated correlation matrix
%	L,U - LU decomposition of PsiDot
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
theta=10.^x;
n=size(X,1);
k=size(X,2);
one=[ones(n,1); zeros(n*k,1)];

% build upper half of correlation matrix
% pre-allocate memory
Psi=zeros(n,n); 
for i=1:n
	for j=i+1:n
		Psi(i,j)=exp(-sum(theta.*(X(i,:)-X(j,:)).^2)); % build top right corner
	end
end
Psi=Psi+Psi'+eye(n)+eye(n).*eps;
% build differentiated matrices
% Pre-allocate memory
PsiDot=zeros((k+1)*n,(k+1)*n);
% Build upper half of PsiDot
for l=1:k
 for m=l:k
  if l==1
   % Build upper half of dPsidX
   for i=1:n
    for j=i+1:n
     PsiDot(i,m*n+j)=2*theta(m)*(X(i,m)-X(j,m))*Psi(i,j);
    end
   end
   % Add upper and lower halves
   PsiDot(1:n,m*n+1:(m+1)*n)=PsiDot(1:n,m*n+1:(m+1)*n)...
                             -PsiDot(1:n,m*n+1:(m+1)*n)';	
  end
  if m==l
   % Build upper half of d2PsidX^2
   for i=1:n
    for j=i+1:n
     PsiDot(l*n+i,m*n+j)=...
         (2*theta(l)-4*theta(l)^2*(X(i,l)-X(j,l))^2)*Psi(i,j);
    end
   end
   % Add half diagonal
   PsiDot(l*n+1:(l+1)*n,m*n+1:(m+1)*n)=PsiDot(l*n+1:(l+1)*n,m*n+1:(m+1)*n)...
                                       +eye(n).*theta(l);
  else
  % Build upper half of d2PsidXdX
   for i=1:n
    for j=i+1:n
     PsiDot(l*n+i,m*n+j)=-4*theta(l)*theta(m)*(X(i,l)-X(j,l))*(X(i,m)-X(j,m))*Psi(i,j);
    end
   end		
   % Add upper and lower halves
   PsiDot(l*n+1:(l+1)*n,m*n+1:(m+1)*n)=PsiDot(l*n+1:(l+1)*n,m*n+1:(m+1)*n)...
                                       +PsiDot(l*n+1:(l+1)*n,m*n+1:(m+1)*n)';
  end
 end
end
% Add upper and lower halves to Psi
PsiDot=[Psi zeros(n,k*n);zeros(k*n,(k+1)*n)]+PsiDot+PsiDot';

%%% concentrated log-likelihood calculation %%%
% Cholesky factorisation
[U,p]=chol(PsiDot);
if p>0
    NegLnLike=1e4;
else

% twice sum logs of diagonal to find ln(abs(det(Psi)))
lndetPsi=2*sum(log(abs(diag(U))));
%%% use back-substitution of LU decomposition instead of inverse 
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
sigma_sqr=((y-one*mu)'*(U\(U'\(y-one*mu))))/((k+1)*n);
NegLnLike=-1*(-(((k+1)*n)/2)*log(sigma_sqr) - 0.5*lndetPsi);
end
