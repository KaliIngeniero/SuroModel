function NegLogConExpImp=reintconstrainedei(x)
% NegLogConExpImp=reintconstrainedei(x)
%
% Calculates the negative of the log of the 
% constrained expected improvement at x
%
% Inputs:
%	x - 1 x k vetor of of design variables
%
% Global variables used:
%   ObjectiveInfo - structure
%	ConstraintInfo - structured cell array
%
% Outputs:
%	NegLogConExpImp - scalar -log(E[I(x)]P[F(x)])
%
% [NegCondLnLike]=regcondlikelihood(x)
%
% Calculates the negative of the conditional 
% ln-likelihood at x(2:k+1) using Kriging regression
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
% A I J Forrester 13th December 2007
%
global ModelInfo
global ObjectiveInfo
global ConstraintInfo

% Calculate unconstrained E[I(x)]
ModelInfo=ObjectiveInfo;
ModelInfo.Option='NegLogExpImp';
NegLogExpImp=reintpredictor(x);

% Calculate P[F(x)] for each constraint
for i=1:size(ConstraintInfo,2)
    ModelInfo=ConstraintInfo{i};
    ModelInfo.Option='NegProbImp';
    NegProbImp(i)=reintpredictor(x);
end

% Calculate E[I(x)]P[F(x)] (add realmin before taking logs)
NegLogConExpImp=-(-NegLogExpImp+sum(log10(-NegProbImp+realmin))); % changed from "NegLogConExpImp=-(-NegLogExpImp+sum(log10(-NegProbImp+1e-50)));" (September 2009)

