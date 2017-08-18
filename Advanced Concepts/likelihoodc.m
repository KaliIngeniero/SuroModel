function NegLnLikec=likelihoodc(x)
% NegLnLikec=likelihoodc(x)
%
% Calculates the negative of the concentrated ln-likelihood of the cheap data
%
% Inputs:
%	x - vector of log(theta) parameters
%
% Global variables used:
%	ModelInfo.Xc - n x k matrix of sample locations
%	ModelInfo.yc - n x 1 vector of observed data
%
% Outputs:
%	NegLnLike - concentrated log-likelihood *-1 for minimising
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
Xc=ModelInfo.Xc;
yc=ModelInfo.yc;
nc=size(Xc,1); 
thetac=10.^x;
p=2;  % added p definition (February 10)
one=ones(nc,1);
PsicXc=zeros(nc,nc);
for i=1:nc
	for j=i+1:nc
		PsicXc(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xc(j,:)).^p)); % abs added (February 10)
	end
end
PsicXc=PsicXc+PsicXc'+eye(nc)+eye(nc).*eps; 
[U,p]=chol(PsicXc);
if p>0
    NegLnLikec=100;
else
LnDetPsicXc=2*sum(log(abs(diag(U))));
muc=(one'*(U\(U'\yc)))/(one'*(U\(U'\one)));
SigmaSqrc=(yc-one.*muc)'*(U\(U'\(yc-one.*muc)))/nc;
NegLnLikec=-1*(-(nc/2)*log(SigmaSqrc)-0.5*LnDetPsicXc);
end
