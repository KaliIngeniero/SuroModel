function [C,Ceq]=likelihoodratiotest(x)
% [C,Ceq]=likelihoodratiotest(x)
%
% Performs a likelihood ratio test to evaluate whether a
% hypothesised point at x falls within a confidence limit
%
% Inputs:
%	x - 1 x k+1 or 2k+1 vetor of theta parameter and hypothesised
%       point 
%       x(1) is the hypothesised function value
%       x(2:k+1) are the theta paramters
%       x(k+2:2k+1) is the location of the hypothesised point
%                   (can be set using ModelInfo.x instead)
%
% Global variables used:
%	ModelInfo.X - n x k matrix of sample locations
%	ModelInfo.y - n x 1 vector of observed data
%   ModelInfo.Theta - 1 x k vector of log(theta)
%   ModelInfo.U - n x n Cholesky factorisation of Psi
%   ModelInfo.Confidence - scalar confidence limit
%
% Outputs:
%	C - negative if hypothesis within confidence limit
%   Ceq - empty equality constraint
%
% Calls:
%   predictor.m, condlikelihood.m
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
k=size(ModelInfo.X,2);
if length(x)<2*k+1
    x(k+2:2*k+1)=ModelInfo.x;
end

% Prediction at imputed point
ModelInfo.Option='Pred';
ModelInfo.Goal=predictor(x(k+2:2*k+1));
% L_O
[L0]=condlikelihood([ModelInfo.Theta x(k+2:2*k+1)]);

% L_HYP
ModelInfo.Goal=x(1);
[L_HYP]=condlikelihood(x(2:2*k+1));

% Value of C must be less than zero
C=2*(-L0+L_HYP)-chi2inv(ModelInfo.Confidence,k);
% If C is complex set to positive (violation)
if ~isreal(C)
    C=1;
end
% Empty equality constraint
Ceq=[];

