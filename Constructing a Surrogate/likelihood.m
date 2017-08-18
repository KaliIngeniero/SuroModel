function [NegLnLike,Psi,U]=likelihood(x)
% [NegLnLike,Psi,U]=likelihood(x)
%
% Calculates the negative of the concentrated ln-likelihood
%
% Inputs:
%	x - vector of log(theta) parameters
%
% Global variables used:
%	ModelInfo.X - n x k matrix of sample locations
%	ModelInfo.y - n x 1 vector of observed data
%
% Outputs:
%	NegLnLike - concentrated log-likelihood *-1 for minimising
%	Psi - correlation matrix
%	U - Choleski factorisation of correlation matrix
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
p=2;  % added p definition (February 10)
n=size(X,1);
one=ones(n,1);

% Pre-allocate memory
Psi=zeros(n,n);
% Build upper half of correlation matrix
for i=1:n
	for j=i+1:n
		Psi(i,j)=exp(-sum(theta.*abs(X(i,:)-X(j,:)).^p)); % abs added (February 10)
	end
end

% Add upper and lower halves and diagonal of ones plus 
% small number to reduce ill-conditioning
Psi=Psi+Psi'+eye(n)+eye(n).*eps; 

% Cholesky factorisation
[U,p]=chol(Psi);

% Use penalty if ill-conditioned
if p>0
    NegLnLike=1e4;
else
    
    % Sum lns of diagonal to find ln(abs(det(Psi)))
    LnDetPsi=2*sum(log(abs(diag(U))));

    % Use back-substitution of Cholesky instead of inverse
    mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
    SigmaSqr=((y-one*mu)'*(U\(U'\(y-one*mu))))/n;
    NegLnLike=-1*(-(n/2)*log(SigmaSqr) - 0.5*LnDetPsi);
end
