function buildcokriging
% Builds co-Kriging model, following parameter estimation. Run prior to
% using cokrigingpredictor.m
%
% Global variables used:
%	ModelInfo.Xc - n x k matrix of cheap sample locations
%	ModelInfo.yc - n x 1 vector of cheap observed data
%	ModelInfo.Xe - n x k matrix of expensive sample locations
%	ModelInfo.ye - n x 1 vector of expensive observed data
%   ModelInfo.Thetac - 1 x k vector of log(theta) of cheap model
%   ModelInfo.Thetad - 1 x k vector of log(theta) of difference model
%   ModelInfo.rho - scalar scaling parameter
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

one=ones(ne+nc,1);
y=[yc; ye];
PsicXc=zeros(nc,nc);
for i=1:nc
	for j=i+1:nc
		PsicXc(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xc(j,:)).^p));
	end
end
ModelInfo.PsicXc=PsicXc+PsicXc'+eye(nc)+eye(nc).*eps; 
ModelInfo.UPsicXc=chol(ModelInfo.PsicXc);

PsicXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsicXe(i,j)=exp(-sum(thetac.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo.PsicXe=PsicXe+PsicXe'+eye(ne)+eye(ne).*eps; 
ModelInfo.UPsicXe=chol(ModelInfo.PsicXe);


PsicXcXe=zeros(nc,ne);
for i=1:nc
	for j=1:ne
		PsicXcXe(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo.PsicXcXe=PsicXcXe; % Deleted "+[zeros(nc-ne,ne);eye(ne)].*eps" November 2010
ModelInfo.PsicXeXc=ModelInfo.PsicXcXe';


PsidXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsidXe(i,j)=exp(-sum(thetad.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo.PsidXe=PsidXe+PsidXe'+eye(ne)+eye(ne).*eps;
ModelInfo.UPsidXe=chol(ModelInfo.PsidXe);


ModelInfo.muc=(ones(nc,1)'*(ModelInfo.UPsicXc\(ModelInfo.UPsicXc'\yc)))/(ones(nc,1)'*(ModelInfo.UPsicXc\(ModelInfo.UPsicXc'\ones(nc,1))));
ModelInfo.d=ye-rho.*yc(end-ne+1:end);
ModelInfo.mud=(ones(ne,1)'*(ModelInfo.UPsidXe\(ModelInfo.UPsidXe'\ModelInfo.d)))/(ones(ne,1)'*(ModelInfo.UPsidXe\(ModelInfo.UPsidXe'\ones(ne,1))));

ModelInfo.SigmaSqrc=(yc-ones(nc,1).*ModelInfo.muc)'*(ModelInfo.UPsicXc\(ModelInfo.UPsicXc'\(yc-ones(nc,1).*ModelInfo.muc)))/nc; 
ModelInfo.SigmaSqrd=(ModelInfo.d-ones(ne,1).*ModelInfo.mud)'*(ModelInfo.UPsidXe\(ModelInfo.UPsidXe'\(ModelInfo.d-ones(ne,1).*ModelInfo.mud)))/ne; 

ModelInfo.C=[ModelInfo.SigmaSqrc*ModelInfo.PsicXc rho*ModelInfo.SigmaSqrc*ModelInfo.PsicXcXe;
rho*ModelInfo.SigmaSqrc*ModelInfo.PsicXeXc rho^2*ModelInfo.SigmaSqrc*ModelInfo.PsicXe+ModelInfo.SigmaSqrd*ModelInfo.PsidXe];
ModelInfo.UC=chol(ModelInfo.C);

ModelInfo.mu=(one'*(ModelInfo.UC\(ModelInfo.UC'\y)))/(one'*(ModelInfo.UC\(ModelInfo.UC'\one)));

