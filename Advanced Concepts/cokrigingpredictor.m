function metric=cokrigingpredictor(x)
% metric=cokrigingpredictor(x)
%
% Calculates the co-Kriging prediction, RMSE, -log(E[I(x)]) or -log(P[I(x)])
%
% Inputs:
%	x - 1 x k vetor of design variables
%
% Global variables used:
%	ModelInfo.Xc - n x k matrix of cheap sample locations
%	ModelInfo.yc - n x 1 vector of cheap observed data
%	ModelInfo.Xe - n x k matrix of expensive sample locations
%	ModelInfo.ye - n x 1 vector of expensive observed data
%   ModelInfo.Thetac - 1 x k vector of log(theta) of cheap model
%   ModelInfo.Thetad - 1 x k vector of log(theta) of difference model
%   ModelInfo.SigmaSqrc - scalar variance of cheap model
%   ModelInfo.SigmaSqrd - scalar variance of expensive model
%   ModelInfo.rho - scalar scaling parameter
%   ModelInfo.UC - nc x nc Cholesky factorisation of C
%   ModelInfo.Option - string: 'Pred', 'RMSE', 'NegLogExpImp' or 'NegProbImp'
%
% Outputs:
%	metric - prediction, RMSE, -log(E[I(x)]) or -log(P[I(x)]), determined
%	by ModelInfo.option
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
Xe=ModelInfo.Xe;
Xc=ModelInfo.Xc;
ye=ModelInfo.ye;
yc=ModelInfo.yc;
ne=size(Xe,1); 
nc=size(Xc,1); 
thetad=10.^ModelInfo.Thetad;
thetac=10.^ModelInfo.Thetac;
p=2;  % added p definition (February 10)
rho=ModelInfo.rho;
one=ones(nc+ne,1);
y=[yc; ye];
cc=ones(nc,1);
for i=1:nc
	cc(i)=rho*ModelInfo.SigmaSqrc*exp(-sum(thetac.*abs(Xc(i,:)-x).^p));
end
cd=ones(ne,1);
for i=1:ne
	cd(i)=rho^2*ModelInfo.SigmaSqrc*exp(-sum(thetac.*abs(Xe(i,:)-x).^p))+ModelInfo.SigmaSqrd*exp(-sum(thetad.*abs(Xe(i,:)-x).^p));
end
c=[cc;cd];
f=ModelInfo.mu+c'*(ModelInfo.UC\(ModelInfo.UC'\(y-one.*ModelInfo.mu))); 
if strcmp(ModelInfo.Option,'Pred')==0
    SSqr=rho^2*ModelInfo.SigmaSqrc+ModelInfo.SigmaSqrd-c'*(ModelInfo.UC\(ModelInfo.UC'\c));
    s=abs(SSqr)^0.5;
    if strcmp(ModelInfo.Option,'RMSE')==0
            yBest=min(ye);
        if strcmp(ModelInfo.Option,'NegProbImp')==1
            ProbImp=0.5+0.5*erf((1/sqrt(2))*((yBest-f)/s));
        else
            EITermOne=(yBest-f)*(0.5+0.5*erf((1/sqrt(2))*((yBest-f)/s)));
            EITermTwo=s*(1/sqrt(2*pi))*exp(-(1/2)*((yBest-f)^2/SSqr));
            ExpImp=log10(EITermOne+EITermTwo+realmin); % changed from "ExpImp=log10(EITermOne+EITermTwo+eps);" (September 2009)
        end
    end
end
if strcmp(ModelInfo.Option,'Pred')==1
metric=f;
elseif strcmp(ModelInfo.Option,'RMSE')==1
metric=s;
elseif strcmp(ModelInfo.Option,'NegLogExpImp')==1
metric=-ExpImp;
elseif strcmp(ModelInfo.Option,'NegProbImp')==1
metric=-ProbImp;
end



