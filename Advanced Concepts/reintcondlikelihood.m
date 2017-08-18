function NegCondLnLike=reintcondlikelihood(x)
% NegCondLnLike=reintcondlikelihood(x)
%
% Calculates the negative of the conditional 
% ln-likelihood at x(2:k+1) using Kriging reinterpolation
%
% Inputs:
%	x - 1 x 2k vetor of theta and hypothesised point
%
% Global variables used:
%	ModelInfo.X - n x k matrix of sample locations
%	ModelInfo.y - n x 1 vector of observed data
%   ModelInfo.Lambda - scalar regulrization parameter
%   ModelInfo.U - n x n Cholesky factorisation of Psi
%   ModelInfo.Goal - scalar goal
%
% Outputs:
%	NegCondLnLike - scalar negative ln-likelihood
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
k=size(X,2);
n=size(X,1);
% replace observed data with regressed points
ModelInfo.Option='Pred';
for i=1:n
    y(i,1)=regpredictor(ModelInfo.X(i));
end
theta=10.^x(1:k);
p=2;  % added p definition (February 10)
if length(x)==k*2
    xHyp=x(k+1:2*k);
else
    xHyp=ModelInfo.x;
end
% build upper half of correlation matrix
Psi=zeros(n,n);
for i=1:n
	for j=i+1:n
		Psi(i,j)=exp(-sum(theta.*abs(X(i,:)-X(j,:)).^p)); % build top right corner
	end
end
% add upper and lower halves and diagonal of ones 
% plus small number to reduce ill-conditioning
Psi=Psi+Psi'+eye(n)+eye(n).*1e-15; 
[U,p]=chol(Psi);
if p>0
    NegCondLnLike=100;
else
one=ones(n,1);
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
psi=ones(n,1);
for i=1:n
	psi(i)=exp(-sum(theta.*abs(X(i,:)-xHyp).^p));
end
m=one*mu+psi*(ModelInfo.Goal-mu);
C=Psi-psi*psi';
% concentrated log-likelihood calculation
[U,p]=chol(C);
if p>0
    NegCondLnLike=100;
else
% sum logs of diagonal to find ln(abs(det(Psi)))
LnDetC=2*sum(log(abs(diag(U))));
SigmaSqr=((y-m)'*(U\(U'\(y-m))))/n;
NegCondLnLike=-1*(-(n/2)*log(SigmaSqr) - 0.5*LnDetC);
end
end

